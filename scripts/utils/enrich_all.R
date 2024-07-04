enrich_all <- function(x,key=NULL,output_file=NA,is.plot=F,
                 show.item=30, go.width=10, pvalue=0.1){
  
  # get genes and gene id -------------------
  prefix <- gsub(".*/","",gsub(".txt","",x))
  gene <- readLines(x)
  tmp <- AnnotationDbi::select(org.Mm.eg.db, gene, "ENTREZID", "SYMBOL") %>% tidyr::drop_na()
  gene.id <- tmp[,2]
  if (is.na(output_file)){
    output_file <- paste0(prefix,"_enrich_all.xlsx")
  }
  
  # GO term over-representation test-------------------
  ego <- enrichGO(gene = gene, 
                  OrgDb = org.Mm.eg.db, ont = "BP",keyType = "SYMBOL",
                  pvalueCutoff = pvalue, qvalueCutoff = 0.5
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
                   minGSSize = 2, pvalueCutoff = pvalue, qvalueCutoff = 0.5)
  # # Reactome over-representation analysis ----------------------------
  # rr <- enrichPathway(gene=gene.id, organism = "mouse", pvalueCutoff = pvalue, readable=TRUE)
  
  # MSigDb-Hallmark dataset over-representation analysis ----------------------------
  m_t2g <- msigdbr::msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, entrez_gene)
  hh <- enricher(gene=gene.id, pvalueCutoff = pvalue, TERM2GENE=m_t2g)
  
  # MSigDb-C7 immunology dataset over-representation analysis ----------------------------
  m_t2g <- msigdbr::msigdbr(species = "Mus musculus", category = "C7") %>% dplyr::select(gs_name, entrez_gene)
  c7 <- enricher(gene=gene.id, pvalueCutoff = pvalue, TERM2GENE=m_t2g)
  
  # write result to table
  # o <- rbind(ego2@result, kk2@result, rr@result, hh@result, c7@result)
  if (!is.null(kk)){
    kk2 <- setReadable(kk, 'org.Mm.eg.db', 'ENTREZID')
    o <- rbind(ego2@result, kk2@result, hh@result, c7@result)
  }
  o <- rbind(ego2@result, hh@result, c7@result)
  
  # if (!file.exists(output_file)){
  xlsx::write.xlsx(o,output_file,row.names = F)
  # }else{
  #   xlsx::write.xlsx(o,output_file,sheetName = prefix2, row.names = F, append = T)
  # }

  # (optional) generate dot plot
  if (is.plot){
    ## GO dot plot
    p1 <- dotplot(ego2, showCategory=show.item)+ggtitle("GO")
    pdf(file = paste0(prefix,"_go.pdf"),width = go.width, height = 10)
    print(p1)
    dev.off()
    
    ## KEGG dot plot
    if (nrow(kk)>2){
      p2 <- dotplot(kk, showCategory=show.item)+ggtitle("KEGG")
      pdf(file = paste0(prefix,"_kegg.pdf"), width = 6, height = 8)
      print(p2)
      dev.off()
    }
  }
}
