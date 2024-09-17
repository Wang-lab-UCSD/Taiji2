gsea <- function(x,key=NULL,output_file=NA,is.plot=F,show.item=30, go.width=10, pvalue=0.1, qvalue=0.5){
  
  # get genes and gene id -------------------
  prefix <- gsub(".txt","",x)
  prefix2 <- sub(".*/_","",prefix)
  gene <- readLines(x)
  gene.id <- AnnotationDbi::select(org.Mm.eg.db, gene, "ENTREZID", "SYMBOL")[,2]
  
  if (is.na(output_file)){
    output_file <- paste0(prefix,"_GO_KEGG.xlsx")
  }
  
  # GO term over-representation test-------------------
  ego <- enrichGO(gene = gene, 
                  OrgDb = org.Mm.eg.db, ont = "BP",keyType = "SYMBOL",
                  pvalueCutoff = pvalue,qvalueCutoff = qvalue
  )
  ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  
  # select specific GO terms (optional)
  if (!is.null(key)){
    key.term <- which(ego$ID%in%key)
    if (length(key.term)>=1){
      o <- ego[key.term]
      write.table(o, paste0(prefix,"_specificGOterms.txt"),quote=F, row.names = F)
    }
  }

  # KEGG over-representation test ---------------------------
  kk <- enrichKEGG(gene = gene.id, organism = "mmu", 
                   minGSSize = 2, pvalueCutoff = pvalue, qvalueCutoff = qvalue)
  if (!is.null(kk) && nrow(kk)>=1){
    kk2 <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    common_cols <- intersect(names(kk2@result),names(ego2@result))
    o <- rbind(kk2@result[,common_cols],ego2@result[,common_cols])
  }else{
    o <- ego2@result
  }
  
  # write result to table
  if (!file.exists(output_file)){
    xlsx::write.xlsx(o,output_file,sheetName = prefix2)
  }else{
    xlsx::write.xlsx(o,output_file,sheetName = prefix2, row.names = F, append = T)
  }

  # (optional) generate dot plot
  if (is.plot){
    ## GO dot plot
    p1 <- dotplot(ego2, showCategory=show.item)+ggtitle("GO")
    pdf(file = paste0(prefix,"_go.pdf"),width = go.width, height = 10)
    print(p1)
    dev.off()
    
    ## KEGG dot plot
    if (!is.null(kk) && nrow(kk)>2){
      p2 <- dotplot(kk, showCategory=show.item)+ggtitle("KEGG")
      pdf(file = paste0(prefix,"_kegg.pdf"), width = 6, height = 8)
      print(p2)
      dev.off()
    }
  }
}
