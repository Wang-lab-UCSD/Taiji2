# unpairedtTest uses unpaired t test select tissue-specific TFs -----------------------

# -------------- arguments -------------------------------------------------------
#   Data: rows are genes, columns are samples, could be original or normalized data
#   n:    indices of all target cells
#   p:    p-value cut-off
#   f:    "single" means the identifed TFs have higher PageRank in target group;
#         "double" means that using abs(log2FoldChange)
#   d:    log2 fold change cut-off
#   ps:   pseudocount added to Data, to avoid the "Inf" when log2

# -------------- main ------------------------------------------------------------

unpairedtTest <- function(Data,n,p=0.001,d=1,f="single",ps=0,type="original",test.type="two.sided"){
  
  Data <- Data[rowSums(Data[, -1])>0, ] # remove the row if rowSum is 0
  pvalue <- rep(0,nrow(Data))
  avg <- rep(0,nrow(Data))
  
  if (type=="original"){
    Data <- log2(Data+ps)
  }
  
  ## check if one-sample or two-sample t.test----
  if (length(n)==1){
    ## one-sample t.test
    for (i in 1:nrow(Data)){
      tryCatch({
        pvalue[i] <- t.test(as.numeric(as.character(unlist(Data[i,-n]))), 
                            mu = Data[i,n], alternative = "less")$p.value
        avg[i] <- Data[i,n]-mean(as.numeric(as.character(unlist(Data[i,-n]))))
      },error=function(e){})
    }
    
  }else if (length(n)>=1){
    ## two-sample t.test
    for (i in 1:nrow(Data)){
      tryCatch({
        pvalue[i] <- t.test(as.numeric(as.character(unlist(Data[i,n]))),
                            as.numeric(as.character(unlist(Data[i,-n]))),
                            alternative = test.type)$p.value
        avg[i] <- mean(as.numeric(as.character(unlist(Data[i,n]))))-
          mean(as.numeric(as.character(unlist(Data[i,-n]))))
        
      },error=function(e){})
    }
  }
  ## order and filter----
  x <- cbind(Data, 'P-value'=pvalue, 'log2 fold change' = avg)
  
  if (f == "single"){
    x <- x[which(x$'P-value' <= p & x$'P-value' > 0 & x$'log2 fold change' >= d), ]
  }else{
    x <- x[which(x$'P-value' <= p & x$'P-value' > 0 &  abs(x$'log2 fold change') >= d), ]
  }
  
  # results are ordered first by increasing p-value and second by decreasing log2FC
  x <- x[order(x$'P-value',x$'log2 fold change',decreasing = c(FALSE, TRUE)),]
  return (x)
  
}


