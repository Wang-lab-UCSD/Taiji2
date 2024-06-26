#### differentiation wave ---------------------------------------------------
library(ggplot2)
source("utils/plotPCA.R")
source("utils/utils.R")
source("utils/findOptimal.R")
source("utils/outputWave.R")
source("utils/gsea.R")

#### input lists ---------
#### 1. Data: Taiji PageRank score, can be dataframe or file, required
#### 2. wavedf: pre-defined differentiation path, required for visualization
#### 3. cv: cut-off of cumulative variance explained used to select number of princial components, default is 0.75
#### 4. 
# ##### 0. give pre-defined differentiation path
# wavedf <- data.frame(x = c(1,2,2,3,3,4,4,5,5), 
#                      y = c(3,5,1,3,1,3,0,4,2),
#                      samplename = c("Naive","TE","TexProg","MP","TexInt","TRM","TexTerm","TEM","TCM"),
#                      labelposx = c(1,2,2,3,3,4,4,5,5),
#                      labelposy = c(3,5,1,3,1.4,3,0,4,2)-0.4)

#### 1. kmeans clustering of TFs based on PageRank -------------------------------------------------
## pca plot shows cumulative variance explained vs # of components
setwd(maindir)
Data_normed <- zscore(Data)
plotPCA(Data_normed)

# manual check the plot and determine PC number (needs to change)
PCNo <- 10
data_reduced <- as.data.frame(prcomp(Data_normed, rank. = PCNo)$x)

## kmeans 
findOptimal(data_reduced)
# manually check the plot and determine number of k and distance metric (needs to change)
clusterNo <- 7
Cluster <- Kmeans(data_reduced, centers = clusterNo, nstart = 25,iter.max = 50, method = "pearson")

# extract tfs from same cluster
df2 <- as.data.frame(Cluster$cluster)
names(df2) <- "cluster"
write.csv(df2, "cluster_result.csv")

##### 2. transcriptional waves visualization------------------
subdir <- paste0("differentiationWave/wave",clusterNo)
createDir(subdir)

df2 <- read.csv("cluster_result.csv", row.names = 1)
lapply(c(1:clusterNo), function(x) outputWave(x,wavedf,df2 = df2))

##### 3. over-representation GO and KEGG terms-----------
L <- list.files(path="./", pattern = "c[0-9]+.txt")
lapply(L, gsea)

# ##### 4. get ranked list for each TF cluster-------------
# get_ranked_list <- function(x, foreground){
#   geneL <- readLines(paste0("c",x,".txt"))
#   df <- Data_normed[geneL,] %>% tibble::rownames_to_column("id") %>% reshape2::melt(id.vars="id") %>% 
#     dplyr::mutate(group=gsub("[0-9\\.].*","",sapply(strsplit(as.character(variable),"_"),`[`,1))) %>%
#     dplyr::mutate(group2=ifelse(group%in%foreground, "foreground","background")) %>%
#     aggregate(value~group2+id, data=., mean) %>% tidyr::pivot_wider(names_from = group2, values_from = value)%>%
#     dplyr::mutate(diff=foreground-background) %>% arrange(-diff)
#   write.csv(df,paste0("c",x,"_pagerank.csv"))
#   df %>% pull(id) %>% writeLines(paste0("c",x,"_ranked.txt"))
# }
# 
# get_ranked_list(x = 1, foreground = c("Naive","TE","MP","MP","TRM","TEM","TCM"))
# get_ranked_list(x = 2, foreground = c("TexInt","TexProg","TexTerm"))
# get_ranked_list(x = 3, foreground = c("MP","TE"))
# get_ranked_list(x = 4, foreground = c("TEM","TCM"))
# get_ranked_list(x = 5, foreground = c("Naive"))
# get_ranked_list(x = 6, foreground = c("TRM"))
# get_ranked_list(x = 7, foreground = c("Naive","TE","MP","MP","TRM","TEM","TCM"))
# 
# ##### 5. summarize GO results-------------
# L <- list.files(pattern = "c.*GO_KEGG.xlsx")
# sheets <- lapply(L, function(x) {read.xlsx(x, sheetIndex = 1)})
# names(sheets) <- lapply(L, function(x) {sub("_GO.*","",x)})
# write_xlsx(sheets, "differentiation_wave_TFs_GO_terms.xlsx")
