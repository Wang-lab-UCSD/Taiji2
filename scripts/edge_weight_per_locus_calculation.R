library(dplyr)
library(plyranges)
library(valr)
setwd("~/Taiji/mouse_immune/mouse11_allcombn/output/ATACSeq/BigWig_norm/batch_corrected/")

# read locus info
# genes <- c('RBPJ','NR4A3', 'JUNB')
genes <- c('RBPJ','IRF4','CD86','NR4A3','EGR2','KIF11','MYC','JUNB','PHLDA1','CREB3L2','SERPINB1A')
df <- read.csv('~/Taiji/mouse_immune/mouse11_allcombn/output/Network/edges_binding_subset_v4.csv') |> filter(X.END_ID %in% genes) |> select(-X.TYPE)

# change sample name
mapping_df <- read.table("~/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/group_file.txt", sep=' ') %>% tibble::rownames_to_column("oldName") |> select(oldName,newName)
df <- df |> left_join(mapping_df, by = c("sample"="oldName")) |> select(-sample)
names(df)[8] <- 'sample'


tfs <- unique(df$X.START_ID)


# prepare for overlap
df2 <- df |> select(chr, start.int, end.int) %>% distinct()
names(df2)[1:3] <- c("chrom","start","end")

# read expression value
rna <- read.csv("~/Taiji/mouse_immune/mouse11_allcombn/output/post-analysis/rna_all.csv", header = T, row.names = 1) 
rownames(rna) <- toupper(rownames(rna))
rna <- rna[,grepl('TRM|TexTerm', names(rna))] %>% tibble::rownames_to_column('gene') |> tidyr::pivot_longer(cols=!gene, names_to = 'sample', values_to = 'rna_expr')

# read bigwig files
bw_files = list.files(pattern = '^(TexTerm|TRM_liver|TRM_IEL).*bw')

# main function
# formula for calculating edge weight per locus: sqrt(TF_expression * peak_intensity * binding_affinity)
get_ew <- function(file_bw,file_ew,file_rna){
    # file_bw: bigwig file contains peak intensity. File path: /path_to_taiji_output/ATACSeq/BigWig/sample_name_ATAC_rep0.bw
    # file_ew: edges_binding file contains binding affinity. File path: /path_to_taiji_output/Network/sample_name/edges_binding.csv
    # file_rna: gene expression file. File path: /path_to_taiji_output/RNASeq/expression_profile.tsv
    gr <- as.data.frame(read_bigwig(file_bw)) |> dplyr::select(chrom, start, end, score)
    # gr1 <- as.data.frame(gr) %>% dplyr::select(seqnames, start, end, score)
    # names(gr1)[1] <- "chrom"

    df_peak <- valr::bed_intersect(df2,gr, suffix = c("_locus", "_score"))
    df_peak_s <- df_peak |> filter(.overlap>=0.5*(end_locus-start_locus)) |> 
            select(chrom, start_locus, end_locus, score_score) |> group_by(chrom, start_locus, end_locus) %>% summarize(mean_score = mean(score_score))
    dff <- df |> left_join(df_peak_s, by = join_by(chr==chrom, start.int==start_locus, end.int==end_locus)) |> 
        left_join(rna, by = join_by(sample, X.START_ID==gene)) |> mutate(edge.weight = sqrt(affinity*mean_score*rna_expr))
    write.csv(dff, paste0('edge_weight_per_locus_',sub('bw','csv',file)), row.names=F, quote=F)
}

tmp <- lapply(bw_files, function(x) get_ew(file = x, df = df, df2 = df2, rna=rna))

# # get average value
# fs <- list.files(pattern = 'edge_weight_per_locus_T.*csv')
# tmp <- do.call('rbind',lapply(fs, function(x) {na.omit(read.csv(x)) |> mutate(file=x,group=strsplit(x,'_')[[1]][5]) |> filter(grepl(strsplit(x,'_')[[1]][5],sample))}))
# output <- tmp |> group_by(X.START_ID, X.END_ID, chr, start.int, end.int, annotation, group) %>% summarize(mean_edge_weight=mean(edge.weight))
# write.csv(output,"edge_weight_per_locus_all.csv",quote=F,row.names=F)