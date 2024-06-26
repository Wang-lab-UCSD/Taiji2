#### import packages---------------------------------------------------------
maindir <- "~/Documents/allcombn/"
setwd(maindir)
set.seed(42)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(gtable)
library(grid)
library(ggforce)
library(dplyr)
library(data.table)
library(ggpubr)
library(VennDiagram)
library(ggrepel)
library(amap)
library(gridExtra)
library(png)
library(tibble)
library(writexl)
library(cowplot)
library(Hmisc)
library(corrplot)
library(formattable)
library(kableExtra)
library(readr)
library(purrr)
library(ggfortify)
library(reshape2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(xlsx)

fl.sources <- list.files(paste0(maindir,"scripts/functions"), full.names = T)
tmp <- sapply(fl.sources,source)

theme_set(theme_get() + theme(text = element_text(family = 'Arial')))


#### import Data-------------------------------------------------------------
# pagerank scores
Data <- read.table("GeneRanks.tsv", check.names = F)
rownames(Data) <- firstup(tolower(rownames(Data)))
Data <- Data[!rownames(Data) %in% c("Dnmt1"),] # remove DNMT1


# group info
group_sorted <- read.table('group_file.txt',header = T, row.names = 1)
group_sorted$oldName <- rownames(group_sorted)
rownames(group_sorted) <- group_sorted$newName
names(Data) <- group_sorted$newName[match(names(Data),group_sorted$oldName)]

# rna data
rna_seq <- read.table("expression_profile.tsv", check.names = F)
rna_seq <- rna_seq[complete.cases(rna_seq),]
names(rna_seq) <- group_sorted$newName[match(names(rna_seq),group_sorted$oldName)]
rna_seq <- rna_seq[,rownames(group_sorted)]

# update new TFnames
TFnames <- intersect(rownames(rna_seq),rownames(Data))
rna_seq <- rna_seq[TFnames,]
rna_seq <- rna_seq[rowSums(rna_seq)>0,]
TFnames <- intersect(rownames(rna_seq),rownames(Data))
Data <- Data[TFnames,]
#### default scale
Data_normed <- zscore(Data)
#### scale to [0,1]
Data_normed_2 <- scaleData(Data)


#### kmeans to cluster TFs -------------------------------------------------
## pca plot shows cumulative variance explained vs # of components
setwd(maindir)
plotPCA(Data_normed)
# manual check the plot and determine PC number
PCNo <- 10
data_reduced <- as.data.frame(prcomp(Data_normed, rank. = PCNo)$x)

## kmeans 
findOptimal(data_reduced)
clusterNo <- 7
Cluster <- Kmeans(data_reduced, centers = clusterNo, nstart = 25,iter.max = 50, method = "pearson")

# extract tfs from same cluster
df2 <- as.data.frame(Cluster$cluster)
names(df2) <- "cluster"

##### 0. differentiation wave ---------------------------------------------------
# wavedf <- data.frame(x = c(1,2,2,3,3,4,4,4,5), 
#                  y = c(3,5,1,3,1,4,2,0,3),
#                  samplename = c("Naive","TE","TexProg","MP","TexInt","TEM","TRM","TexTerm","TCM"),
#                  labelposx = c(1,2,2,3,3.4,4,4,4,5),
#                  labelposy = c(3,5,1,3,1.4,4,2,0,3)-0.4)
wavedf <- data.frame(x = c(1,2,2,3,3,4,4,5,5), 
                     y = c(3,5,1,3,1,3,0,4,2),
                     samplename = c("Naive","TE","TexProg","MP","TexInt","TRM","TexTerm","TEM","TCM"),
                     labelposx = c(1,2,2,3,3,4,4,5,5),
                     labelposy = c(3,5,1,3,1.4,3,0,4,2)-0.4)

##### 1. transcriptional waves for clusters------------------
subdir <- paste0("differentiationWave/wave",clusterNo)
createDir(subdir)

df2 <- do.call("rbind", lapply(1:7, function(x){
  data.frame(TF=readLines(paste0("c",x,".txt")),cluster=x)
})) %>% tibble::column_to_rownames("TF")

write.csv(df2, "cluster_result.csv")
write.csv(Data_normed, "pagerank.csv")
df2 <- read.csv("cluster_result.csv", row.names = 1)
Data_normed <- read.csv("pagerank.csv", row.names = 1)


lapply(c(1:clusterNo), 
       function(x) outputWave(x,wavedf,df2 = df2))

##### 2. over-representation GO and KEGG terms-----------
L <- list.files(path="./", pattern = "c[0-9]+.txt")
lapply(L, gsea)

##### 3. get ranked list for each TF cluster-------------
get_ranked_list <- function(x, foreground){
  geneL <- readLines(paste0("c",x,".txt"))
  df <- Data_normed[geneL,] %>% tibble::rownames_to_column("id") %>% reshape2::melt(id.vars="id") %>% 
    dplyr::mutate(group=gsub("[0-9\\.].*","",sapply(strsplit(as.character(variable),"_"),`[`,1))) %>%
    dplyr::mutate(group2=ifelse(group%in%foreground, "foreground","background")) %>%
    aggregate(value~group2+id, data=., mean) %>% tidyr::pivot_wider(names_from = group2, values_from = value)%>%
    dplyr::mutate(diff=foreground-background) %>% arrange(-diff)
  write.csv(df,paste0("c",x,"_pagerank.csv"))
  df %>% pull(id) %>% writeLines(paste0("c",x,"_ranked.txt"))
}

get_ranked_list(x = 1, foreground = c("Naive","TE","MP","MP","TRM","TEM","TCM"))
get_ranked_list(x = 2, foreground = c("TexInt","TexProg","TexTerm"))
get_ranked_list(x = 3, foreground = c("MP","TE"))
get_ranked_list(x = 4, foreground = c("TEM","TCM"))
get_ranked_list(x = 5, foreground = c("Naive"))
get_ranked_list(x = 6, foreground = c("TRM"))
get_ranked_list(x = 7, foreground = c("Naive","TE","MP","MP","TRM","TEM","TCM"))
###### 3.1 write ranked list to excel -------------
setwd("~/Documents/allcombn/differentiationWave/wave7/")
L <- list.files(pattern = ".*ranked.txt")
sheets <- lapply(L, readLines %>% as.data.frame)
names(sheets) <- lapply(L, function(x) {sub("_.*","",x)})
write_xlsx(sheets, "differentiation_wave_ranked_TFs.xlsx")


##### 4. GO analysis for each cluster regulatees-------------
