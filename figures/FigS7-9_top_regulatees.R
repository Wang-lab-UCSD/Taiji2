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
sapply(fl.sources,source)

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

# remove two TRM.lung samples from Hayward dataset -------
group_sorted <- group_sorted[!grepl("TRM.*Hayward",rownames(group_sorted)),]
Data <- Data[,rownames(group_sorted)]

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
Data_normed <- as.data.frame(t(scale(t(as.matrix(Data)))))


n1 = length(unique(group_sorted$Group)) 
mycolors <- list(subsets = c("#A6A6A6", "#FF0000", "#0007F4", "#FFA75C", "#00AAFE", 
                             "#009051", "#FEB5B5", "#E6CFFF", "#8AAA75"))
names(mycolors$subsets) <- unique(group_sorted$Group)

annotation_col = data.frame(
  subsets = group_sorted$Group
)
rownames(annotation_col) = colnames(Data)

#### analysis of regulatees----------------------------------
#### 0. set up dir -------------------
subdir <- "network/"
dir.create(file.path(maindir, subdir), showWarnings = FALSE)
setwd(file.path(maindir, subdir))
samples <- c("texterm", "trm")
key <- c("CD74","NT5E","CXCR3","CXCR6","CX3CR1","TGFBR")
TF <- readLines("TFs.txt")

#### 1. reimport rna_seq ---------------
rna_seq <- read.table("../expression_profile.tsv", check.names = F)
# row.names(rna_seq) <- toupper(rownames(rna_seq))
rna_seq <- rna_seq[complete.cases(rna_seq),]
names(rna_seq) <- group_sorted$newName[match(names(rna_seq),group_sorted$oldName)]
rna_seq <- rna_seq[,rownames(group_sorted)]
rna_seq <- rna_seq[rowSums(rna_seq)>0,]

rna_seq_normed <- as.data.frame(t(scale(t(as.matrix(rna_seq)))))
rna_seq_log <- log2(rna_seq+1)

#### 6. expression plot of regulatees --------------------------
get_top_regulatee_heatmap_by_expression <- function(TF,top=100){
  L <- list.files(path = "result_20220518", pattern = paste0(TF,".*texterm.txt"), full.names = T, ignore.case = T)
  x <- L[1]
  r <- readLines(x)
  df <- rna_seq_normed[r,] %>% tibble::rownames_to_column("id") %>% reshape2::melt(id.vars="id") %>% 
    dplyr::mutate(group=gsub("[0-9\\.].*","",sapply(strsplit(as.character(variable),"_"),`[`,1))) %>%
    aggregate(value~group+id, data=., mean) 
  
  selected_genes <- df %>% reshape(idvar = "id", timevar = "group", direction = "wide") %>%
    mutate(diff.texterm.te = value.TexTerm-value.TE) %>% 
    # dplyr::arrange(-value.TexTerm, value.TE, -ratio.texterm.te)
    top_n(top, diff.texterm.te) %>% arrange(-diff.texterm.te) %>%
    column_to_rownames("id") %>%
    dplyr::select(1:(ncol(.)-1)) %>% dplyr::rename_with(function(x) gsub("value.","",x))
  
  ## heatmap of expression
  p3 <- pheatmap(selected_genes, fontsize = 7, show_rownames = T,
                 angle_col = 45, show_colnames = T,
                 cellwidth = 10, cellheight = 7,
                 # clustering_distance_cols = 'correlation', 
                 # clustering_distance_rows = 'correlation', 
                 # clustering_method = 'average',
                 # treeheight_row = 0, treeheight_col = 10,
                 # annotation_col = annotation_col, annotation_colors = mycolors, 
                 cluster_rows = F, cluster_cols = F,
                 border_color = NA,
                 filename = paste0("hp_rna_normed_top_",top,"_",TF,"_regulatees.pdf"))
  
  
}
lapply(c("Nfil3","Zscan20","Jdp2"), get_top_regulatee_heatmap_by_expression)


