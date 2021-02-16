# ==== ML preprocess =====

timestamp <- Sys.time()
library(microbiome)
library(caret)
library(plyr)
library(recipes)
library(dplyr)
require(randomForest);require(caret);require(caretEnsemble);require(ggplot2);require(tidyverse)

rm(list=ls()[!ls() %in% c("dftops", "pstodf", "stool.core.ps", "ps.stool")]) # from DMM_fit5.R
seed=81
set.seed(seed)


# ==== custom ML function for association rule mining (ARM); pass to caretList ======
customARM <- list(label = "Random Forest Rule-Based Model",
                  library = c("randomForest", "inTrees", "plyr"),
                  type = c("Classification", "Regression"),
                  parameters = data.frame(parameter = c("mtry","maxdepth"),
                                          class = rep("numeric",2),
                                          label = c("#Randomly Selected Predictors","Maximum Rule Depth")),
                  grid = function(x, y, len = NULL, search = "grid"){
                    if(search == "grid") {
                      out <- data.frame(mtry = caret::var_seq(p = ncol(x),
                                                              classification = is.factor(y),
                                                              len = len),
                                        maxdepth = (1:len)+1)
                    } else {
                      out <- data.frame(mtry = sample(1:ncol(x), size = len, replace = TRUE),
                                        maxdepth = sample(1:15, size = len, replace = TRUE))
                    }
                  },
                  loop = function(grid) {
                    loop <- plyr::ddply(grid, c("mtry"),
                                        function(x) c(maxdepth = max(x$maxdepth)))
                    submodels <- vector(mode = "list", length = nrow(loop))
                    for(i in seq(along = loop$maxdepth)) {
                      index <- which(grid$mtry == loop$mtry[i])
                      trees <- grid[index, "maxdepth"]
                      submodels[[i]] <- data.frame(maxdepth = trees[trees != loop$maxdepth[i]])
                    }
                    list(loop = loop, submodels = submodels)
                  },
                  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                    if(!is.data.frame(x) | inherits(x, "tbl_df"))
                      x <- as.data.frame(x, stringsAsFactors = TRUE)
                    RFor <- randomForest::randomForest(x, y, mtry = param$mtry, ...)
                    treeList <- inTrees::RF2List(RFor)
                    exec <- inTrees::extractRules(treeList,x, maxdepth=param$maxdepth, ntree = RFor$ntree)
                    ruleMetric <- inTrees::getRuleMetric(exec,x,y)
                    ruleMetric <- inTrees::pruneRule(ruleMetric,x,y)
                    ruleMetric <- inTrees::selectRuleRRF(ruleMetric,x,y)
                    out <- list(model = inTrees::buildLearner(ruleMetric,x,y))
                    if(!last) {
                      out$rf <- treeList
                      out$x <- x
                      out$y <- y
                      out$trees <- RFor$ntree
                    }
                    out
                  },
                  predict = function(modelFit, newdata, submodels = NULL) {
                    if(!is.data.frame(newdata) | inherits(newdata, "tbl_df"))
                      newdata <- as.data.frame(newdata, stringsAsFactors = TRUE)
                    out <- inTrees::applyLearner(modelFit$model, newdata)
                    if(modelFit$problemType == "Regression") out <- as.numeric(out)
                    if(!is.null(submodels)) {
                      tmp <- vector(mode = "list", length = nrow(submodels) + 1)
                      tmp[[1]] <- if(is.matrix(out)) out[,1] else out
                      for(i in seq(along = submodels$maxdepth)) {
                        exec <- inTrees::extractRules(modelFit$rf,
                                                      modelFit$x,
                                                      maxdepth=submodels$maxdepth[i],
                                                      ntree = modelFit$trees)
                        ruleMetric <- inTrees::getRuleMetric(exec,modelFit$x,modelFit$y)
                        ruleMetric <- inTrees::pruneRule(ruleMetric,modelFit$x,modelFit$y)
                        ruleMetric <- inTrees::selectRuleRRF(ruleMetric,modelFit$x,modelFit$y)
                        mod <- inTrees::buildLearner(ruleMetric,modelFit$x,modelFit$y)
                        tmp[[i+1]] <- inTrees::applyLearner(mod, newdata)
                        if(modelFit$problemType == "Regression") tmp[[i+1]] <- as.numeric(tmp[[i+1]])
                      }
                      out <- tmp
                    }
                    out
                  },
                  prob = NULL,
                  predictors = function(x, ...) {
                    split_up <- strsplit(x$model[,"condition"], "&")
                    
                    isolate <- function(x) {
                      index <- gregexpr("]", x, fixed = TRUE)
                      out <- NULL
                      for(i in seq_along(index)) {
                        if(all(index[[i]] > 0)) {
                          tmp <- substring(x[i], 1, index[[i]][1])
                          tmp <- gsub("(X)|(\\[)|(\\])|(,)|( )", "", tmp)
                          tmp <- tmp[tmp != ""]
                          out <- c(out, as.numeric(tmp))
                        }
                      }
                      as.numeric(unique(out))
                    }
                    
                    var_index <- unique(unlist(lapply(split_up, isolate)))
                    if(length(var_index) > 0) x$xNames[var_index] else NULL
                  },
                  varImp = function(object, ...) {
                    split_up <- strsplit(object$model[,"condition"], "&")
                    
                    isolate <- function(x) {
                      index <- gregexpr("]", x, fixed = TRUE)
                      out <- NULL
                      for(i in seq_along(index)) {
                        if(all(index[[i]] > 0)) {
                          tmp <- substring(x[i], 1, index[[i]][1])
                          tmp <- gsub("(X)|(\\[)|(\\])|(,)|( )", "", tmp)
                          tmp <- tmp[tmp != ""]
                          out <- c(out, as.numeric(tmp))
                        }
                      }
                      as.numeric(unique(out))
                    }
                    
                    var_index <- lapply(split_up, isolate)
                    
                    vars_dat <- lapply(var_index,
                                       function(x, p) {
                                         out <- rep(0, p)
                                         if(length(x) > 0) out[x] <- 1
                                         out
                                       },
                                       p = length(object$xNames))
                    vars_dat <- do.call("rbind", vars_dat)
                    colnames(vars_dat) <- object$xNames
                    freqs <- as.numeric(object$model[,"freq"])
                    vars_dat <- vars_dat * freqs
                    var_imp <- apply(vars_dat, 2, sum)
                    out <- data.frame(Overall = as.vector(var_imp))
                    rownames(out) <- names(var_imp)
                    out
                  },
                  levels = function(x) x$obsLevels,
                  tags = c("Random Forest", "Ensemble Model", "Bagging",
                           "Implicit Feature Selection", "Rule-Based Model"),
                  sort = function(x) x[order(x[,"maxdepth"]),])


#remove 0 sum features so distance calculations downstream are allowed
stool.core.ps.new <- prune_taxa(taxa_sums(stool.core.ps) > 0, stool.core.ps) #  keep only taxa with positive sums
stool.core.ps.rel <- microbiome::transform(stool.core.ps, "compositional")
rowSums(stool.core.ps.rel@otu_table)
stool.core.ps@sam_data["550.L1S274.s.1.sequence",]
stool.core.ps@sam_data["550.L1S334.s.1.sequence",]
stool.core.ps@sam_data["550.L2S57.s.2.sequence",]
stool.core.ps<-subset_samples(stool.core.ps, sample_names(stool.core.ps) !=("550.L1S274.s.1.sequence"))
stool.core.ps<-subset_samples(stool.core.ps, sample_names(stool.core.ps) !=("550.L1S334.s.1.sequence"))
stool.core.ps<-subset_samples(stool.core.ps, sample_names(stool.core.ps) !=("550.L2S57.s.2.sequence"))
summarize_phyloseq(stool.core.ps)
stool.core.ps # main ps

stool.core.ps.rel <- microbiome::transform(stool.core.ps, "compositional")

stool.core.ml<-pstodf(stool.core.ps.rel) # append host_subject_id, use df for ml

#class balance check
rowSums(stool.core.ml[1:10])
stool.core.ml$met_var<-as.factor(stool.core.ml$met_var)
percentage <- prop.table(table(stool.core.ml$met_var)) * 100
cbind(freq=table(stool.core.ml$met_var), percentage=percentage)

K<-list(stool.core.ml)

#K[[1]] <- K[[1]][sample(nrow(K[[1]])),]

ind.stool <- sample(2,nrow(K[[1]]),replace=TRUE,prob=c(0.7,0.3))
train.stool <- K[[1]][ind.stool==1,]
test.stool <- K[[1]][ind.stool==2,]


seeds <- vector(mode = "list", length = nrow(stool.core.ml) + 1)
seeds <- lapply(seeds, function(x) 1:20)
#grid <- expand.grid(mtry = c(1, 7), maxdepth = c(3, 7))
grid <- expand.grid(mtry = 3, maxdepth = 7) # max depth = rule length
ctrl1 <- trainControl(method = "cv", number = 3, returnResamp = "all",
                      verboseIter = TRUE,
                      seeds = seeds)

set.seed(849)
tr.mod <- train(met_var ~ ., data = train.stool, 
                method = customARM,
                trControl = ctrl1,
                preProc = c("center", "scale"),
                tuneGrid = grid,
                ntree = 1003)

tr.mod[["results"]];tr.mod[["finalModel"]][["model"]];tr.mod[["finalModel"]][["xNames"]]

test.stool2 <- test.stool[,1:(ncol(test.stool)-1)]; Y <- test.stool[,"met_var"]

preds <- predict(tr.mod, test.stool2) # pred on hold out using train model

#accuracy /balance of test model TODO: create ROC
length(preds)
pred<-as.character(preds)
obs<-as.character(test.stool$met_var)
obs

logic.preds<-as.data.frame(cbind(pred, obs))
logic.preds$logic <- logic.preds$pred == logic.preds$obs

perc <- prop.table(table(logic.preds$logic)) * 100
perc<- cbind(freq=table(logic.preds$logic), percentage=perc)

perc #86% accuracy

#show class balance of test.
percentage <- prop.table(table(test.stool$met_var)) * 100
cbind(freq=table(test.stool$met_var), percentage=percentage)


# ===== ML for predicted probability distribution =====

library(dplyr) # data manipulation
library(caret) # base model-building
#library(DMwR) # potential class imbalance
library(purrr) # functional programming (map)
library(pROC) # AUC calculations
library(PRROC) # Precision-Recall curve calculations
require(randomForest);require(caret);require(caretEnsemble);require(ggplot2);require(tidyverse);require(xgboost)

#define custom RF so we can tune optimal hyperparameters
#no submodels here
rf.predProbs <- list(type = "Classification", library = "randomForest", loop = NULL)
rf.predProbs$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
rf.predProbs$grid <- function(x, y, len = NULL, search = "grid") {}
rf.predProbs$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
rf.predProbs$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
rf.predProbs$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
rf.predProbs$sort <- function(x) x[order(x[,1]),]
rf.predProbs$levels <- function(x) x$classes


# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3,
                        index = createFolds(stool.core.ml$met_var, 10), #use full dataset not train subset
                        summaryFunction = twoClassSummary,
                        classProbs=T,
                        savePredictions = T,
                        preProc=c("center","scale"))
metric<- "logLoss"
tunegrid <- expand.grid(.mtry=c(2:4), .ntree=c(100,203,501,1002,2341))
set.seed(seed)
custom <- train(met_var~., data=stool.core.ml, method=rf.predProbs, metric=metric, tuneGrid=tunegrid, trControl=control)
summary(custom)
plot(custom);custom 

custom$bestTune
#identify ~generalizable model

# ==== Final model =====
#retune with best fit
set.seed(782)
#extract probs for each obs rather than for each fold etc.
tunegrid <- expand.grid(.mtry=2)

#retune with best fit
final.fit <- caretList(
  met_var~., data=stool.core.ml,
  trControl=control,
  ntree=2341,
  tuneGrid=tunegrid,
  # methodList=c('rf', 'adaboost', 'earth', 'xgbDART', 'svmRadial')
  methodList=c('rf')
)

final.fit

options(digits=3)
final.fit.pp<-extractProb(final.fit) # ***** has pred probs of final fit?
final.fit.pp

# ==== VarImp custom (rules) ====
require(caret)
varImp(custom)
plot(varImp(custom))


final.fit[["rf"]][["finalModel"]][["importance"]] # MDA



#get caretFit scaled to 100 not in mean decrease accuracy
library(caret) 
library(dplyr) 
library(caretEnsemble) 
library(lattice) 
#library(gridExtra) 


#if using a list of models - perhaps TODO: text xgb RRF and RF ; works with base ensemble
plot_importance <- function(importance_list, imp, algo_names) {
  importance <- importance_list[[imp]]$importance
  model_title <- algo_names[[imp]]
  #dotplot v s heatmap
  if (ncol(importance) < 2) { 
    importance %>%
      as.matrix() %>%
      dotplot(main = model_title)
  } else { 
    importance %>%
      as.matrix() %>%
      levelplot(xlab = NULL, ylab = NULL, main = model_title, scales = list(x = list(rot = 45)))
  }
}

importance <- lapply(final.fit, varImp)
importance_plots <- list()
for (imp in seq_along(importance)) {
  # importance_plots[[imp]] <- plot(importance[[imp]])
  importance_plots[[imp]] <- plot_importance(importance_list = importance, imp = imp, algo_names = names(final.fit))
}
do.call("grid.arrange", c(importance_plots))

