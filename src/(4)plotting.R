
# ===== make predicted probability graphs ====

rm(list=ls()[!ls() %in% c("dftops", "pstodf", "stool.core.ps", "ps.stool", "final.fit.pp", "tr.mod", "final.fit")]) # from ML_stool.R

#check if OTU table is same order as predprobs index
if(!identical(row.names(final.fit.pp), row.names(stool.core.ps@otu_table))) stop("Observations not identical") # sanity checks out

#TODO: lets add a logical column to tell which samples fell below predicted probability threshold

#define colours for plotting
number.taxa<-length(tax_table(stool.core.ps))

microbiomeseq_cols <- function(){
  colours <- c("#019964", "#73DB9A", "#FFBE00","#F4E38F","#663084","#CCCCFF","#EF2B2B","#FF9999","#00CCFF","#CCFFFF",
               "#B50488","#FFCCFF",grey.colors(1000));
  return(colours)
}

colours <- microbiomeseq_cols()

#melt the ps object so we can plot long

library(reshape2);library(ggplot2)
#cbind the abundance information onto the predicted probability index
transform <- microbiome::transform
abund.rel.ps <- transform(stool.core.ps, "compositional")
abund.rel.ps
rowSums(abund.rel.ps@otu_table)

#append an ordered sample number vector so we can plot by sample
abund.rel.vec<-as.factor(1:259)
abund.rel.df<-as.data.frame(abund.rel.ps@otu_table)
abund.rel.df$sample.num<-abund.rel.vec

pp.abund<-cbind(final.fit.pp, abund.rel.df) # MERGE THE pp from ML and otu table from ps object

#sanity checks
if(!identical(row.names(pp.abund), row.names(stool.core.ps@otu_table))) stop("Observations not identical") # sanity checks out
final.fit.pp["550.L1S89.s.1.sequence",]
pp.abund["550.L1S89.s.1.sequence",]
abund.rel.df["550.L1S89.s.1.sequence",]
abund.rel.ps@otu_table["550.L1S89.s.1.sequence",]
stool.core.ps@otu_table["550.L1S89.s.1.sequence",]
stool.core.ps@sam_data["550.L1S89.s.1.sequence",] # should be 79, 3, 50, etc.
#sanity checks out

stool.core.ps.rel<-abund.rel.ps


# Male: M3
stool.core.ps.rel.M3<-subset_samples(stool.core.ps.rel, host_subject_id=="M3")
RA.M3 <- plot_composition(stool.core.ps.rel.M3, 
                          sample.sort = rownames(stool.core.ps.rel.M3@otu_table), 
                          otu.sort = NULL, 
                          x.label = "sample.num", 
                          plot.type = "barplot", 
                          verbose = FALSE) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "sample.num M3", 
       y = "Relative abundance",
       title = "RA M3", 
       subtitle = "stool",
       caption = "Caption text.") + 
  scale_fill_brewer("ASV") +scale_fill_manual(values=colours[1:(number.taxa+1)])


#* identify which obvervations the predprobs fall below threshold of .75

M3.pp<-pp.abund %>%
  select(F4, M3, obs, pred, sample.num) %>%
  filter(obs=="M3")
dim(M3.pp)
M3.pp.vec<-as.factor(1:128)
M3.pp$new.sample.n<-M3.pp.vec


#assign threshold cutoff .75
M3.pp.75<-M3.pp %>%
  filter(M3 < .75)
M3.pp.75  # samples falling below threshold: 5,30,117

#add intercept lines where pp <.threshold are.
M3.pp.p<-ggplot(M3.pp) + geom_line(aes(x=sample.num, y=M3, colour="white", group =obs), lwd =0.1) + theme_biome_utils() +
  geom_hline(yintercept = .75, 
             color = "black", size=0.3) + theme(
               # Hide panel borders and remove grid lines
               panel.border = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               # Change axis line
               axis.line = element_line(colour = "black")
             )

# add intercepts to the observations that fall below .75 threshold selected
RA.M3<- RA.M3 +
  geom_vline(xintercept = 5, 
             color = "black", size=0.3) +
  geom_vline(xintercept = 30,
             color = "black", size=0.3) +
  geom_vline(xintercept = 117, 
             color = "black", size=0.3) + theme_biome_utils()



#Female: F4
stool.core.ps.rel.F4<-subset_samples(stool.core.ps.rel, host_subject_id=="F4")
RA.F4 <- plot_composition(stool.core.ps.rel.F4, 
                          sample.sort = rownames(stool.core.ps.rel.F4@otu_table), 
                          otu.sort = NULL, 
                          x.label = "sample.num", 
                          plot.type = "barplot", 
                          verbose = FALSE) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "sample.num F4", 
       y = "Relative abundance",
       title = "RA F4", 
       subtitle = "stool",
       caption = "Caption text.") + 
  scale_fill_brewer("ASV") +scale_fill_manual(values=colours[1:(number.taxa+1)])



F4.pp<-pp.abund %>%
  select(F4, M3, obs, pred, sample.num) %>%
  filter(obs=="F4")
#append new sample.numbers for reduced dimensions upon filtering
dim(F4.pp)
F4.pp.vec<-as.factor(1:131)
F4.pp$new.sample.n<-F4.pp.vec
#assign threshold cutoff .75
F4.pp.75<-F4.pp %>%
  filter(F4 < .75)

F4.pp.75 #
F4.pp.p<-ggplot(F4.pp) + geom_line(aes(x=sample.num, y=F4, colour="white", group =obs), lwd =0.1) + theme_biome_utils() +
  geom_hline(yintercept = .75, 
             color = "black", size=0.3) + theme(
               # Hide panel borders and remove grid lines
               panel.border = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               # Change axis line
               axis.line = element_line(colour = "black")
             )

# add intercepts to the observations that fall below .75 threshold selected
RA.F4 <- RA.F4 +
  geom_vline(xintercept = 2, 
             color = "black", size=0.3) +
  geom_vline(xintercept = 4,
             color = "black", size=0.3) +
  geom_vline(xintercept = 12, 
             color = "black", size=0.3) +
  geom_vline(xintercept = 68, 
             color = "black", size=0.3) +
  geom_vline(xintercept = 90, 
             color = "black", size=0.3) + theme_biome_utils()

require(patchwork)
RA.F4 + RA.M3 + F4.pp.p + M3.pp.p


# ==== ordinate phyloseq object to see class distribution====
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling  

#rm(list=ls()[!ls() %in% c("dftops", "pstodf", "stool.core.ps", "ps.stool", "final.fit.pp", "tr.mod", 
#                         "F4.pp.75", "M3.pp.75")]) # from ML_stool.R


stool.core.ps.rel
#ps1.log <- transform_sample_counts(ps1, function(x) log(1 + x))

ordu.pcoa <- ordinate(stool.core.ps.rel, "PCoA", "bray")


CommunityTypeColors.use<- c(
  "#dd5f72", "#5f7be0")

pcoa.community.type <- plot_ordination(stool.core.ps.rel, 
                                       ordu.pcoa, color="host_subject_id") 
pcoa.community.type <- pcoa.community.type + ggtitle("PCoA Bray Community TYpe")
pcoa.community.type
pcoa.community.type <- pcoa.community.type + theme_biome_utils()
print(pcoa.community.type)
pcoa.community.type2<-pcoa.community.type + facet_grid(~body_site) + scale_color_manual(values = CommunityTypeColors.use) + 
  geom_point(shape = 21, colour = "black", fill = "white", size = 0.1, stroke = 0.2)
pcoa.community.type2
print(pcoa.community.type2 + stat_ellipse())



# ==== Inference the 'misclassified' samples =======

#M3
require(phyloseq)
M3.pp.75
rownames(M3.pp.75)
stool.core.ps.rel
stool.core.ps.rel.75.M3 <- subset(otu_table(stool.core.ps.rel), rownames(otu_table(stool.core.ps.rel)) %in% c("550.L1S137.s.1.sequence", 
                                                                                                              "550.L1S229.s.1.sequence", 
                                                                                                              "550.L6S329.s.6.sequence"))
stool.core.ps.rel.75.M3 <- merge_phyloseq(stool.core.ps.rel.75.M3, tax_table(stool.core.ps.rel), sample_data(stool.core.ps.rel))
#gives OTU table filtered.

library(microbiomeutilities)
library(microbiome)
stool.core.ps.rel.75.M3.p <- plot_composition(stool.core.ps.rel.75.M3, 
                                              sample.sort = rownames(stool.core.ps.rel.75.M3@otu_table),
                                              otu.sort = NULL, 
                                              x.label = "sample.num", 
                                              plot.type = "barplot", 
                                              verbose = FALSE) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "sample.num M3", 
       y = "Relative abundance",
       title = "RA M3 .75", 
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  scale_fill_brewer("Genus") +scale_fill_manual(values=colours[1:(number.taxa+1)])




#F4
require(phyloseq)
F4.pp.75
rownames(F4.pp.75)
stool.core.ps.rel
stool.core.ps.rel.75.F4 <- subset(otu_table(stool.core.ps.rel), rownames(otu_table(stool.core.ps.rel)) %in% c("550.L1S44.s.1.sequence", 
                                                                                                              "550.L1S36.s.1.sequence", 
                                                                                                              "550.L1S35.s.1.sequence",
                                                                                                              "550.L1S37.s.1.sequence",
                                                                                                              "550.L1S27.s.1.sequence"))
stool.core.ps.rel.75.F4 <- merge_phyloseq(stool.core.ps.rel.75.F4, tax_table(stool.core.ps.rel), sample_data(stool.core.ps.rel))
#gives OTU table filtered.

library(microbiomeutilities)
library(microbiome)
stool.core.ps.rel.75.F4.p <- plot_composition(stool.core.ps.rel.75.F4, 
                                              sample.sort = rownames(stool.core.ps.rel.75.F4@otu_table),
                                              otu.sort = NULL, 
                                              x.label = "sample.num", 
                                              plot.type = "barplot", 
                                              verbose = FALSE) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "sample.num F4", 
       y = "Relative abundance",
       title = "RA F4 .75", 
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  scale_fill_brewer("Genus") +scale_fill_manual(values=colours[1:(number.taxa+1)])



require(patchwork)
stool.core.ps.rel.75.M3.p + stool.core.ps.rel.75.F4.p





# ==== Heatmap of predicted probability quantiles ====

# create a quantile variable to add to ps metadata so we can plot using plot_heatmap_taxa
library(microbiomeutilities)
library(RColorBrewer)
library(ggpubr)

#* EVERYTHING IS IN TERMS OF M3
#need to create a quantile column so we can append to metadata
# we need to know the length of obs that < 0.25 pred probs, ...
#Quantile 0-.25
quantilesM3.1<-pp.abund %>%
  filter(M3 < 0.25) 

dim(quantilesM3.1) #126

#Quantile .25-.5
quantilesM3.2<-pp.abund %>%
  filter(M3 > .25) %>%
  filter(M3 < .50) 

dim(quantilesM3.2) #5

#Quantile .5-.75
quantilesM3.3<-pp.abund %>%
  filter(M3 > .50) %>%
  filter(M3 < .75) 

dim(quantilesM3.3) #3

#Quantile .75-1
quantilesM3.4<-pp.abund %>%
  filter(M3 > .75)  
dim(quantilesM3.4) #125

dim(quantilesM3.1);dim(quantilesM3.2);dim(quantilesM3.3);dim(quantilesM3.4)

quantiles.pp.M3<-rep(c("1","2","3", "4"), times=c(126, 5, 3, 125))

pp.abund.arr.M3<-pp.abund %>%
  arrange(M3)
pp.abund.arr.M3["quantiles.pp.M3"]<-quantiles.pp.M3
pp.abund.arr.M3

#arrange back by sample number now that we appended quantiles.
pp.abund.arr.M3<-pp.abund.arr.M3 %>%
  arrange(sample.num)
pp.abund.arr.M3


#add new set of sample data to ps object
stool.core.ps.pp.q<-abund.rel.ps # core ps RA'd
stool.core.ps.pp.q@otu_table["550.L1S89.s.1.sequence",] 
sample_data(stool.core.ps.pp.q)$ordered <- pp.abund.arr.M3$sample.num
sample_data(stool.core.ps.pp.q)$M3pp <- pp.abund.arr.M3$M3
sample_data(stool.core.ps.pp.q)$F4pp <- pp.abund.arr.M3$F4
sample_data(stool.core.ps.pp.q)$pred <- pp.abund.arr.M3$pred
sample_data(stool.core.ps.pp.q)$obs <- pp.abund.arr.M3$obs
sample_data(stool.core.ps.pp.q)$quantiles.pp <- pp.abund.arr.M3$quantiles.pp

#sanity checks
stool.core.ps.pp.q@otu_table["550.L1S89.s.1.sequence",] 
stool.core.ps.pp.q@sam_data["550.L1S89.s.1.sequence",] 

#sanity checks
if(!identical(rownames(stool.core.ps.pp.q@otu_table), rownames(pp.abund))) stop("string not identical")
pp.abund["550.L1S89.s.1.sequence",] # can remove (sanity)
stool.core.ps.pp.q@otu_table["550.L1S89.s.1.sequence",] # can remove (sanity)
stool.core.ps@otu_table["550.L1S89.s.1.sequence",] # can remove (sanity)
stool.core.ps.pp.q@sam_data["550.L1S89.s.1.sequence",] # can remove (sanity)

#sanity checks out 


#Subset samples we want so we can append later, arranges by sex.
M3_ps.rel <- subset_samples((stool.core.ps.pp.q), host_subject_id =="M3")
F4_ps.rel <- subset_samples((stool.core.ps.pp.q), host_subject_id =="F4")


rel.q1 <- subset_samples(stool.core.ps.pp.q, quantiles.pp =="1")
rel.q2 <- subset_samples(stool.core.ps.pp.q, quantiles.pp=="2")
rel.q3 <- subset_samples(stool.core.ps.pp.q, quantiles.pp=="3")
rel.q4 <- subset_samples(stool.core.ps.pp.q, quantiles.pp=="4")



library(microbiomeutilities)
library(pheatmap)
library(RColorBrewer)

#optional step if required to rename taxa_names
#taxa_names(ps0) <- gsub("d__denovo", "OTU:", taxa_names(ps0))

# create a gradient color palette for abundance

grad_ab <- colorRampPalette(c("#faf3dd","#f7d486" ,"#5e6472"))
grad_ab_pal <- grad_ab(10)
grad_ab_pal
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")
# create a color palette for varaibles of interest
meta_colors <- list(c("M3" = "#5F7BE0",
                      "F4" = "#DD5F72"), 
                    c("1" = "#1B9E77", 
                      "2" = "#D95F02",
                      "3" = "#7570B3",
                      "4" = "#E6AB02"))
# add labels for pheatmap to detect
names(meta_colors) <- c("host_subject_id", "quantiles.pp")

grad_ab <- colorRampPalette(c("#faf3dd","#f7d486" ,"#5e6472"))
grad_ab_pal <- grad_ab(10)

merge.m3.f4.rel<-merge_phyloseq(M3_ps.rel, F4_ps.rel) # merge them to separate F4 and M3
merge.qpp<-merge_phyloseq(rel.q1, rel.q2, rel.q3, rel.q4) # merge them to separate F4 and M3

merge.m3.f4.rel@sam_data$quantiles.pp
merge.m3.f4.rel@otu_table
merge.m3.f4.rel@sam_data$quantiles.pp<-as.factor(merge.m3.f4.rel@sam_data$quantiles.pp)
merge.qpp@sam_data$quantiles.pp<-as.factor(merge.qpp@sam_data$quantiles.pp)
levels(merge.m3.f4.rel@sam_data$quantiles.pp)

#sanity
merge.qpp@otu_table["550.L2S24.s.2.sequence"]
stool.core.ps@otu_table["550.L2S24.s.2.sequence"]
abund.rel.ps@otu_table["550.L2S24.s.2.sequence"]

p.pp.heat.Z <- plot_taxa_heatmap(merge.qpp,
                                 subset.top = 12,
                                 VariableA = c("host_subject_id", "quantiles.pp"),
                                 heatcolors =rev(brewer.pal(6, "Spectral")),
                                 transformation = "Z-OTU",
                                 cluster_rows = T,
                                 cluster_cols = F,
                                 show_colnames = T,
                                 annotation_colors=meta_colors)

p.pp.heat.log10 <- plot_taxa_heatmap(merge.qpp,
                                     subset.top = 12,
                                     VariableA = c("host_subject_id", "quantiles.pp"),
                                     heatcolors =rev(brewer.pal(6, "RdPu")),
                                     transformation = "log10",#x+1
                                     cluster_rows = T,
                                     cluster_cols = F,
                                     show_colnames = T,
                                     annotation_colors=meta_colors)


# ===== Plot mean RA =====
#Some of the same GENUS link back to different ASVs, so we need to make a rank as 'ASV'
library(Biostrings)
rownames(merge.qpp@tax_table)
merge.qpp@tax_table
rankasv<-taxa_names(merge.qpp)
taxa.df<-as.data.frame(merge.qpp@tax_table)
taxa.df
taxa.df$ASV<-rankasv
taxa.df<-as.matrix(taxa.df)
TAX = tax_table(taxa.df)
TAX
merge.qpp2 <- subset(otu_table(merge.qpp))
merge.qpp2
merge.qpp2 <- merge_phyloseq(merge.qpp2, tax_table(TAX), sample_data(merge.qpp))
merge.qpp2@tax_table
merge.qpp@tax_table
#sanity
merge.qpp2@otu_table
merge.qpp2@otu_table["550.L1S99.s.1.sequence",]
merge.qpp@otu_table["550.L1S99.s.1.sequence",]
abund.rel.ps@otu_table["550.L1S99.s.1.sequence",]



grp_abund <- get_group_abundances(merge.qpp2, 
                                  level = "ASV", 
                                  group="host_subject_id",
                                  transform = "compositional")

mycols <- c("brown3", "steelblue","grey50")

# clean names 
grp_abund$OTUID <- gsub("p__", "",grp_abund$OTUID)
grp_abund$OTUID <- ifelse(grp_abund$OTUID == "", 
                          "Unclassified", grp_abund$OTUID)


mean.RA <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=host_subject_id)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("host_subject_id", values=mycols) + # manually specify colors
  theme_biome_utils() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 

mean.RA

# extract OTUID and MRA  to find greatest difference in mean RA between each ASV and each class.
mean.RA.MA<-mean.RA[["data"]][["mean_abundance"]]
mean.RA.OTUID<-mean.RA[["data"]][["OTUID"]]
mean.RA.class.lab<-as.character(mean.RA[["data"]][["host_subject_id"]])
mean.RA.df.val<-cbind(mean.RA.MA)

mean.RA.df.desc<-cbind(mean.RA.MA, mean.RA.OTUID, mean.RA.class.lab)


library(matrixStats)
mean.RA.df.F4<-as.vector(mean.RA.df.val[1:10,])
mean.RA.df.M3<-as.vector(mean.RA.df.val[11:20,])

RA.M3.F4.val<-cbind(mean.RA.df.M3, mean.RA.df.F4)

diffs <- rowDiffs(RA.M3.F4.val)
diffs
mean.RA.OTUID


# ===== pot taxa boxplot ====

library(microbiomeutilities)
library(microbiome)
library(RColorBrewer)
library(ggpubr)

mycols <- c("brown3", "steelblue", "grey50")
taxa.box.ra <- plot_taxa_boxplot(merge.qpp2,
                                 taxonomic.level = "ASV",
                                 top.otu = 12, 
                                 group = "host_subject_id",
                                 add.violin= TRUE,
                                 title = "Top three family", 
                                 keep.other = FALSE,
                                 group.order = c("F4", "M3"),
                                 group.colors = mycols,
                                 dot.size = 1)

taxa.box.ra
#print(taxa.box.ra + theme_biome_utils())

# Plot selected taxa of largest difference in mean RA between each class
select.taxa.diff.MRA <- c("ASV96", "ASV1174", "ASV702", "ASV1319")

# Plot selected taxa of top 5 varIMP on which model? custom final fit or ARM?
select.taxa.top5.varIMP <- c("x", "x", "x", "x")

mycols <- c("brown3", "steelblue", "grey50")

taxa.box.ra <- plot_listed_taxa(merge.qpp2, select.taxa, 
                                group = "host_subject_id",
                                group.order = c("F4","M3"),
                                group.colors = mycols,
                                add.violin = TRUE,
                                violin.opacity = 0.3,
                                dot.opacity = 0.25,
                                box.opacity = 0.25,
                                panel.arrange= "grid")
#> An additonal column Sam_rep with sample names is created for reference purpose
print(taxa.box.ra + ylab("Relative abundance") + scale_y_continuous(labels = scales::percent))

comps <- make_pairs(sample_data(merge.qpp2)$host_subject_id)
print(comps)

taxa.box.ra <- taxa.box.ra + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")
taxa.box.ra + scale_y_continuous(labels = scales::percent)
#+ theme_biome_utils()


# ===== Plot Abundance-Prevalence Relationship - detection threshold =====
#TODO: use abs abund ps object, not compositional
set.seed(81)
prev.det.1 <- plot_abund_prev(merge.qpp2, 
                              label.core = FALSE,
                              color = "ASV", # NA or "blue"
                              mean.abund.thres = 0.01,
                              mean.prev.thres = 0.99,
                              dot.opacity = 0.7,
                              label.size = 3,
                              label.opacity = 1.0,
                              nudge.label=-0.5,
                              bs.iter=9, # increase for actual analysis e.g. 999
                              size = 536, # increase to match your nsamples(asv_ps)
                              replace = TRUE,
                              label.color="#5f0f40")  
prev.det.1 <- prev.det.1 + 
  geom_vline(xintercept = 0.95, lty="dashed", alpha=0.7) + 
  geom_hline(yintercept = 0.01,lty="dashed", alpha=0.7) +
  scale_color_brewer(palette = "Set3") + theme_biome_utils()

prev.det.1


# ===== Detection Threshold plot ====

#detections <- 10^seq(log10(1), log10(max(abundances(merge.qpp))/10), length = 10)
#prevalences <- seq(.05, 1, .05)

rowSums(merge.qpp2@otu_table)
prev.det.2 <- plot_core(merge.qpp2, plot.type = "heatmap", 
                        #prevalences = prevalences,
                        #detections = detections,
                        colours = rev(brewer.pal(5, "Spectral")),
                        min.prevalence = .2, horizontal = TRUE) 
print(prev.det.2 + theme_biome_utils()) # this is detectiom threshold (%) RA over prevalence as a % of total samples. (minimum population prevalence)
#Detection Threshold (Relative Abundance (%))


# above p and p_v go hand in hand could be patchworked.



# ==== Q1-4 RA plots =====
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
merge.qpp2.M3<-subset_samples(merge.qpp2, host_subject_id=="M3")
merge.qpp2.M3@sam_data
merge.qpp2.F4<-subset_samples(merge.qpp2, host_subject_id=="F4")
RA.quantiles.M3 <- plot_composition(merge.qpp2.M3,
                                    taxonomic.level = "ASV",
                                    VariableA = c("host_subject_id", "quantiles.pp"),
                                    sample.sort = "quantiles.pp",
                                    x.label = "quantiles.pp") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_brewer("ASV") +scale_fill_manual(values=colours[1:(number.taxa+1)]) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "RA Q1-4",
       subtitle = "M3",
       caption = "Caption text.")
#theme_ipsum(grid="Y")
RA.quantiles.M3


RA.quantiles.F4 <- plot_composition(merge.qpp2.F4,
                                    taxonomic.level = "ASV",
                                    VariableA = c("host_subject_id", "quantiles.pp"),
                                    sample.sort = "quantiles.pp",
                                    x.label = "quantiles.pp") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_brewer("ASV") +scale_fill_manual(values=colours[1:(number.taxa+1)]) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "RA Q1-4",
       subtitle = "F4",
       caption = "Caption text.")
#theme_ipsum(grid="Y")
merge.qpp2.F4@sam_data


library(patchwork)
(RA.quantiles.M3 + RA.quantiles.F4) # this is ordered, so technically we can manually color the top of the gray lines by quantile



# ==== Mean RA Quantiles ====

# * lets order based on the output of z-scored heatmap quantiles so we can align the 2 plots together.
#heatmap:
ASV1174
ASV2247
ASV171
ASV1611
ASV3615
ASV96
ASV2090
ASV702
ASV1319
ASV1922

# plot mean.plot sex

grp_abund <- get_group_abundances(merge.qpp2, 
                                  level = "ASV", 
                                  group=c("quantiles.pp", "host_subject_id"),
                                  transform = "compositional")
grp_abund
mycols <- c("#1b9e77", "#d95f02", "#7570b3", "#e6ab02")

# clean names 
grp_abund$OTUID <- gsub("p__", "",grp_abund$OTUID)
grp_abund$OTUID <- ifelse(grp_abund$OTUID == "", 
                          "Unclassified", grp_abund$OTUID)

re.ord.heatmap<-c("ASV1174",
                  "ASV2247",
                  "ASV171",
                  "ASV1611",
                  "ASV3615",
                  "ASV96",
                  "ASV2090",
                  "ASV702",
                  "ASV1319",
                  "ASV1922")


meanRA.plot.quant <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=quantiles.pp)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("host_subject_id", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 

meanRA.plot.quant


# we can make a plot like, the sames that were predicted as 
#e.g, quantile 1 for M3 is 0-.25 and then it's 100-whatever the M3 pred prob is for a F4 pred prob.