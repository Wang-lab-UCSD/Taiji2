## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----import packages and functions--------------------------------------------------------------------------------------------------------------------------------------------------------------------
suppressMessages(library(igraph))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(huge))
suppressMessages(library(xlsx))

set.seed(42)
fl.sources <- list.files("../../scripts/utils/", full.names = T)
tmp <- sapply(fl.sources,source)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# sample <- 'TRM' # 'TRM' or 'TexTerm'
# tag <- '_subset_TFs_v2'
# 
# L <- readLines(paste0(sample,'_samples.txt'))
# TFs <- readLines(paste0(sample,'_TFs_for_network_v3.txt')) 
# TFs <- lapply(TFs, toupper)
#  
# 
# edge_weight <- do.call("rbind", lapply(L, function(x) read.table(paste0(x,'/edges_combined.csv'), sep = ",", header = T) %>% select(c(1,2,3)) %>% setNames(c("TF","regulatee","weight")))) %>% dplyr::filter(TF %in% TFs) %>% 
#                group_by(TF, regulatee) %>% dplyr::summarize(weight = mean(weight, na.rm=TRUE)) %>% 
#                reshape2::dcast(regulatee ~ TF, value.var = "weight", fill=0) %>% 
#                tibble::column_to_rownames(var = "regulatee") 
#                                        
# df <- edge_weight %>% rowwise() %>% mutate(sd=sd(c_across(everything()))) %>% ungroup() %>% 
#                 filter(sd>1) %>% select(!matches("MOUSE$|[0-9\\.]{6}$")) %>% select(-sd) %>%
#                 rename_with(firstup)
# write.csv(df,file = paste0("mean_",sample,"_edge_weight",tag,".csv")) # write to file


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sample <- 'TexTerm' ### 'TRM' or 'TexTerm'
tag <- '_subset_TFs_v2'
df <- read.csv(paste0("mean_",sample,"_edge_weight",tag,".csv"), row.names = 1)
colnames(df) <- firstup(colnames(df))
knitr::kable(df[10000:10006,1:6], caption = 'TF-gene edge weight') |> kableExtra::kable_styling(latex_options = 'scale_down')
print(dim(df))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
npn.func <- 'shrinkage'
df.npn <- huge.npn(df, npn.func = npn.func) # default                      
out.npn <- huge(df.npn, method = "glasso", nlambda = 30, lambda.min.ratio = 0.01) # fit glasso model to data
saveRDS(out.npn,paste0("result_copula_",sample,'_',npn.func,tag,".rds"))


## ---- fig.width=8, fig.height=5-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# default visualization provided by huge package
plot(out.npn)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sparsity <- 0.15
idx = which(out.npn$sparsity>=sparsity & out.npn$sparsity<=sparsity+0.01)[1]
print(paste0("lambda: ", out.npn$lambda[idx]))
print(paste0("sparsity: ", out.npn$sparsity[idx]))                       
m = out.npn$icov[[idx]]
rownames(m) <- colnames(df)
colnames(m) <- colnames(df)                       
knitr::kable(head(m[,1:6]), caption = 'sparse TF-TF correlation matrix') |> kableExtra::kable_styling(latex_options = 'scale_down')

## write to file
write.csv(m, paste0("precision_copula_",sample,"_",npn.func,"_sparsity_",sparsity,tag,".csv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- data.frame(row=rownames(m)[row(m)[upper.tri(m)]], 
                 col=colnames(m)[col(m)[upper.tri(m)]], 
                 corr=m[upper.tri(m)])
### remove zero entries
df <- df[df$corr!=0,]
knitr::kable(head(df), caption = 'long form of TF correlation') |> kableExtra::kable_styling(latex_options = 'scale_down')

print(dim(df))
print(paste('percent of correlation score >','0',':',sum(df$corr > 0) *100 / nrow(df),"%"))
print(paste('percent of correlation score <','0',':',sum(df$corr < 0) *100 / nrow(df),"%"))

df$corr = abs(df$corr)
write.csv(df, paste0("correlation_intensity_table_copula_",sample,"_",npn.func,"_sparsity_",sparsity,tag,".csv"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nodes <- read.xlsx('TRM_Tex_TF_list_20240121.xlsx', sheetName = sample) %>% tidyr::drop_na(TF) %>% dplyr::mutate(group=ifelse(is.na(Specificity), Important, Specificity)) %>% select(TF, group)
names(nodes) <- c('name','group')
knitr::kable(head(nodes), caption = 'node meta') |> kableExtra::kable_styling(latex_options = 'scale_down')
print(dim(nodes))   
  


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
g = graph_from_data_frame(df, directed=FALSE, vertices = nodes)
print(g, e=TRUE, v=TRUE)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(42)
res <- 1.1
clustering = cluster_leiden(g, objective_function="modularity",resolution_parameter=res)
print(clustering)
print(sizes(clustering))
print(paste0('modularity:',modularity(g, membership(clustering))))
saveRDS(clustering,paste0('clustering_copula_',sample,'_',npn.func,'_res',res,'_s',sparsity,tag,'.rds'))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## write clustering result
clust <- data.frame(name = names(membership(clustering)),cluster = as.numeric(membership(clustering))) %>%
        dplyr::inner_join(nodes, by = "name") %>%
        dplyr::arrange(group, cluster)
knitr::kable(head(clust), caption = 'community clustering') |> kableExtra::kable_styling(latex_options = 'scale_down')
write.csv(clust,paste0('clustering_membership_copula_',sample,'_',npn.func,'_res',res,'_s',sparsity,tag,'.csv'), row.names=F)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

