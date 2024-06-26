# set up and import libraries--------------------------------
library(dplyr)
library(igraph)
library(ggplot2)
library(data.table)

#### TRM -------------
sample <- "TRM"
maindir <- paste0('/new-stg/home/cong/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/network_analysis/',sample)
setwd(maindir)
getwd()
source('/new-stg/home/cong/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/scripts/function/gsea.R')
source('/new-stg/home/cong/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/scripts/function/utils.R')
func <- 'shrinkage'
res = 1.8
sparsity <- 0.05

##### 1. import graph-------------
g <- readRDS(paste0("graph_",func,"_res",res,"_sparsity_", sparsity,".rds"))
set.seed(42)
clustering = cluster_leiden(g, objective_function="modularity",resolution_parameter=res)

##### 2. customize graph------------
V(g)$cluster <- clustering$membership
# generate colors based on cluster
colors <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(length(clustering))
V(g)$color <- colors[V(g)$cluster]
V(g)$label.color <- ifelse(V(g)$group=="other",colors[V(g)$cluster],"black")
V(g)$frame.color <- "NA"

# size is non-linearly proportional to node degree
V(g)$size <- log2(degree(g))
V(g)$label.cex <- ifelse(V(g)$group == "other",NA,2)
V(g)$label <- ifelse(V(g)$group == "other",NA,V(g)$name)

# set edge color
E(g)$color <- "gray80"
E(g)$size <- 0.1
E(g)$lty <- 2

# only show edges that involve texterm-specific and common TFs
# lty=0 means that the edge won't be shown, default is 1
# it takes 30s
E(g)$lty <- unlist(lapply(E(g), function(x) {
    ifelse(V(g)$group[V(g)$name == ends(g, x)[1]] == "other" | V(g)$group[V(g)$name == ends(g, x)[2]] == "other", 0, 1)
}))

## set the within-module edges to some large weight, and the between module edges to some small weight
## and then choose 'layout_with_fr' to make the grouped layout 
edge.weights <- function(community, network, weight.within = 20, weight.between = 1) {
    bridges <- crossing(communities = community, graph = network)
    weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
    return(weights) 
}
E(g)$weight <- edge.weights(clustering, g)

##### 3. plot graph----------------
set.seed(100)
pdf(paste0("graph_copula_",func,"_res",res,"_s", sparsity,"_degree.pdf"), width = 15, height = 15)
# pdf(paste0("graph_copula_",func,"_res",res,"_s", sparsity,"_degree_no_label.pdf"), width = 15, height = 15)
plot(g, layout=layout_with_fr,
     mark.groups = lapply(unique(V(g)$cluster), function(x) V(g)[V(g)$cluster == x]), 
     mark.col = paste0(colors,'7f'), mark.border = colors)
dev.off()
print(paste0("clustering quality: ",clustering$quality))


#### TexTerm -------------
sample <- "TexTerm"
maindir <- paste0('/new-stg/home/cong/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/network_analysis/',sample)
setwd(maindir)
getwd()
source('/new-stg/home/cong/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/scripts/function/gsea.R')
source('/new-stg/home/cong/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/scripts/function/utils.R')
func <- 'shrinkage'
res = 1.8
sparsity <- 0.05

##### 1. import graph-------------
g <- readRDS("g2.rds")
set.seed(123)
clustering = cluster_leiden(g, objective_function="modularity",resolution_parameter=res)

##### 2. customize graph------------
V(g)$cluster <- clustering$membership
# generate colors based on cluster
colors <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(length(clustering))
V(g)$color <- colors[V(g)$cluster]
V(g)$label.color <- ifelse(V(g)$group=="other",colors[V(g)$cluster],"black")
V(g)$frame.color <- "NA"

# size is non-linearly proportional to node degree
V(g)$size <- log2(degree(g))
V(g)$label.cex <- ifelse(V(g)$group == "other",NA,2)
V(g)$label <- ifelse(V(g)$group == "other",NA,V(g)$name)
# V(g)$label <- NA


# set edge color
E(g)$color <- "gray80"
E(g)$size <- 0.1
E(g)$lty <- 2

# only show edges that involve texterm-specific and common TFs
# lty=0 means that the edge won't be shown, default is 1
# it takes 30s
E(g)$lty <- unlist(lapply(E(g), function(x) {
    ifelse(V(g)$group[V(g)$name == ends(g, x)[1]] == "other" | V(g)$group[V(g)$name == ends(g, x)[2]] == "other", 0, 1)
}))

## set the within-module edges to some large weight, and the between module edges to some small weight
## and then choose 'layout_with_fr' to make the grouped layout 
edge.weights <- function(community, network, weight.within = 50, weight.between = 1) {
    bridges <- crossing(communities = community, graph = network)
    weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
    return(weights) 
}
E(g)$weight <- edge.weights(clustering, g)

##### 3. plot graph----------------
set.seed(42)
pdf(paste0("graph_copula_",func,"_res",res,"_s", sparsity,"_degree.pdf"), width = 15, height = 15)
# pdf(paste0("graph_copula_",func,"_res",res,"_s", sparsity,"_degree_no_label.pdf"), width = 15, height = 15)
plot(g, layout=layout_with_fr,
     mark.groups = lapply(unique(V(g)$cluster), function(x) V(g)[V(g)$cluster == x]), 
     mark.col = paste0(colors,'7f'), mark.border = colors)
dev.off()
print(paste0("clustering quality: ",clustering$quality))

#### write network to excel ----------------
maindir <- "~/Documents/allcombn/network/community/TexTerm/"
setwd(maindir)

df <- read.csv("network.txt", sep = " ", header = F) %>% dplyr::select(V1,V2)
write_xlsx(df, "network_TexTerm.xlsx")

df <- read.csv("network.txt", sep = " ", header = F) %>% dplyr::select(V1,V2) %>% dplyr::filter(V1 %in% c("Nfil3","Prdm1","Zscan20"))
write_xlsx(df, "network_Nfil3_Prdm1_Zscan20.xlsx")

maindir <- "~/Documents/allcombn/network/community/TRM/"
setwd(maindir)

df <- read.csv("network.txt", sep = " ", header = F) %>% dplyr::select(V1,V2)
write_xlsx(df, "network_TRM.xlsx")

df <- read.csv("network.txt", sep = " ", header = F) %>% dplyr::select(V1,V2) %>% dplyr::filter(V1 %in% c("Nfil3","Prdm1"))
write_xlsx(df, "network_Nfil3_Prdm1.xlsx")

