df <- lapply(unique(group_sorted$Group), function(x) read.xlsx('cell_state_specific_255TFs.xlsx', sheetName = x) %>% pull(TF))
names(df) <- unique(group_sorted$Group)

mt <- readLines('MultiTasker_119TFs.txt')
df2 <- lapply(df, function(x) setdiff(x, mt))
df3 <- lapply(df, function(x) intersect(x, mt))

library(UpSetR)
pdf("figS6B_multitasker_upsetR.pdf", width = 10)    
p1 <- UpSetR::upset(fromList(df3), nsets = length(df3), nintersects = NA, 
                    order.by = "freq", text.scale = 1.8)
print(p1)
dev.off()
