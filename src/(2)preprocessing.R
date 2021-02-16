# ==== initiate pipeline ====

# Call Data
rm(list=ls())
seed=81
set.seed(seed)
rm(list=ls()[!ls() %in% "seed"])

# ==== experimental design - load moving pictures dataset ====
library(jeevanuDB)
ps <- moving_pictures
ps@sam_data[1:5,1:5]
ps@otu_table[1:5,1:5]

table(meta(ps)$sample_type, meta(ps)$host_subject_id) # biome distribution M vs F

#find balance between experimental design and individuals M/F - we have to balance by min frequency of samples between each sample_type

table(meta(ps)$host_subject_id, meta(ps)$common_sample_site) # sample site distribution M vs F
# add ASV ID to tax_table
library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
ps@otu_table[1:5, 1:5]
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
taxa_names(ps)[1:3]
ps # main PS object


#balance PS
ps.stool<-subset_samples(ps, sample_type == "stool")
table(meta(ps.stool)$sample_type, meta(ps.stool)$sex)


pstodf<-function(physeq, met_var,...) {
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq) # keep this here sinde we're not balancing.
  physeq.otu = as(otu_table(physeq), "matrix")
  if(taxa_are_rows(physeq)){physeq.otu <- t(physeq.otu)}
  physeq.otu.df = as.data.frame(physeq.otu)
  head(physeq.otu.df)[1:5]
  
  met_var<-physeq@sam_data$host_subject_id # extract response from metadata
  physeq.otu.df$met_var<-met_var #append response to otu df so we can subset with sample_each func
  
  return(physeq.otu.df)
}


L<-list(ps.stool)

L <- list(L[[1]])
L
names(L) <- c("stool")
for(i in 1:length(L)){
  L[[i]]<-pstodf(L[[i]])
}

if(!identical(ps.stool@sam_data$host_subject_id, L[[1]]$met_var)) stop("string not identical")
if(!identical(rownames(ps.stool@sam_data), rownames(L[[1]]))) stop("string not identical")
#passes sanity checks

# randomly sample and balance classes
sample_each <- function(data, var, n = 1L) {
  lvl <- table(data[, var])
  n1 <- setNames(rep_len(n, length(lvl)), names(lvl))
  n0 <- lvl - n1
  idx <- ave(as.character(data[, var]), data[, var], FUN = function(x)
    sample(rep(0:1, c(n0[x[1]], n1[x[1]]))))
  data[!!(as.numeric(idx)), ]
}


L.bal<-list()
table(meta(ps.stool)$host_subject_id, meta(ps.stool)$body_site) 
L.bal[["stool"]]<-sample_each(L[["stool"]], 'met_var', n = c(131,131)) # minimum sample size is 268



#normalize, convert otu table to ps obect for subsequent analysis
dftops<-function(physeq.df, orig_physeq, ...) {
  physeq.df$met_var<-NULL
  OTU = otu_table(physeq.df, taxa_are_rows = FALSE)
  TAX = tax_table(orig_physeq@tax_table)
  SAM = sample_data(orig_physeq@sam_data)
  physeq.new = phyloseq(OTU, TAX, SAM)
  
  physeq.new <- prune_taxa(taxa_sums(physeq.new) > 0, physeq.new) # singletons/0sumfeatures
  #physeq.new = tax_glom(physeq.new, "Phylum", NArm = TRUE, bad_empty=c(NA, "", " ", "\t")) # otu_glom
  #physeq.new <- core(physeq.new, detection = 0, prevalence = .20) # no det/prev thresholding here
  #physeq.new <- microbiome::transform(physeq.new, "compositional") # convert compositional 'after' removing samples
  return(physeq.new)
}

ps.stool<-dftops(L.bal[["stool"]], ps.stool)
ps.stool
stool.df.ml.bal<-pstodf(ps.stool)
percentage <- prop.table(table(stool.df.ml.bal$met_var)) * 100
cbind(freq=table(stool.df.ml.bal$met_var), percentage=percentage)
ps.stool # balanced ps object @ genus level


# ==== core community shared between communbity types ====
library(microbiomeutilities)
library(eulerr)
library(microbiome)


pseq <- ps.stool
pseq
table(meta(pseq)$host_subject_id, useNA = "always")
pseq.rel <- microbiome::transform(pseq, "compositional") #RA


classes <- unique(as.character(meta(pseq.rel)$host_subject_id))
print(classes)


list_core <- c() 
for (n in classes){ 
  
  ps.sub <- subset_samples(pseq.rel, host_subject_id == n) # Choose sample from host_subject_id 'n'
  
  core_m <- core_members(ps.sub,
                         detection = 0.001, # min RA
                         prevalence = 0.75) # 75% of samples # prev
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each sebum samples
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

mycols <- c(nonCRC="#bc55ff", CRC="#ff9357", H="#57d58c") 
plot(venn(list_core),
     fills = mycols)

#keep the shared 10.
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
F4.core<-list_core[["F4"]]
M3.core<-list_core[["M3"]]
F4.core
M3.core

#Shared:
shared.core<-c("ASV96", "ASV171", "ASV702", "ASV1174", "ASV1319","ASV1611","ASV1922", "ASV2090", "ASV2247", "ASV3615")
shared.core

pseq

pseq.df<-pstodf(pseq)

extract.stool<-pseq.df[, shared.core]

dftops<-function(physeq.df, orig_physeq, ...) {
  physeq.df$met_var<-NULL
  OTU = otu_table(physeq.df, taxa_are_rows = FALSE)
  TAX = tax_table(orig_physeq@tax_table)
  SAM = sample_data(orig_physeq@sam_data)
  physeq.new = phyloseq(OTU, TAX, SAM)
  
  #physeq.new <- prune_taxa(taxa_sums(physeq.new) > 0, physeq.new) 
  #physeq.new = tax_glom(physeq.new, "Phylum", NArm = TRUE, bad_empty=c(NA, "", " ", "\t")) 
  #physeq.new <- core(physeq.new, detection = 0, prevalence = .20) # 
  #physeq.new <- microbiome::transform(physeq.new, "compositional") 
  return(physeq.new)
}

pseq
stool.core.ps<-dftops(extract.stool, pseq)
stool.core.ps # main ps


