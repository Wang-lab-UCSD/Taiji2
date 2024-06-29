# import packages and data
set.seed(42)
fl.sources <- list.files("../scripts/utils/", full.names = T)
tmp <- sapply(fl.sources,source)
df2 <- read.csv("cluster_result.csv", row.names = 1)
Data_normed <- read.csv("pagerank.csv", row.names = 1)
clusterNo <- 7
wavedf <- data.frame(x = c(1,2,2,3,3,4,4,5,5), 
                     y = c(3,5,1,3,1,3,0,4,2),
                     samplename = c("Naive","TE","TexProg","MP","TexInt","TRM","TexTerm","TEM","TCM"),
                     labelposx = c(1,2,2,3,3,4,4,5,5),
                     labelposy = c(3,5,1,3,1.4,3,0,4,2)-0.4)

# transcriptional waves for clusters------------------
lapply(c(1:clusterNo), function(x) outputWave(x, wavedf, df2 = df2))


