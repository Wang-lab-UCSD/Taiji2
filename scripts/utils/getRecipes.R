getRecipes <- function(soN, taN, recipe_size=3,
                       rank_cut=20, p=0.001, 
                       method="prod",m="ratio.abs",
                       g=group_sorted,D=Data, R=rna_seq){
  
  sources = rownames(g)[g$Group==soN]
  targets = rownames(g)[g$Group==taN]
  
  getRecipeType <- function(o,s,t,c="TF"){
    o$exp_S <- rowMeans(R[match(o[[c]], rownames(R)),s,drop=F])
    o$exp_T <- rowMeans(R[match(o[[c]], rownames(R)),t,drop=F])
    with(o, ifelse(exp_T>exp_S, "oe", "ko"))
  }
  
  getCandidates <- function(idx_s,idx_t,m){
    name_s <- sources[idx_s]
    name_t <- targets[idx_t]
    message(name_s, " -> ", name_t)
    df = data.frame(diff=D[,name_t]-D[,name_s],
                    diff.abs=abs(D[,name_t]-D[,name_s]),
                    ratio=D[,name_t]/D[,name_s])
    if (sum(df$diff)==0){
      message("pagerank of source and target is exactly the same!")
      writeLines(paste0("pagerank of ",name_s," and ", name_t, "is the same!"),
                 paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"))
      return()
    }
    
    df$ratio.abs = unlist(lapply(df$ratio, function(x) if(x<1){1/x}else{x}))
    rownames(df) = rownames(D)
    df_2 = head(df[order(df[[m]],decreasing = TRUE),], n = rank_cut)
    
    if (recipe_size == 1){
      df_2$rank <- rep(1:nrow(df_2))
      return(df_2)
    }else{
      tmp <- as.data.frame(t(combn(rownames(df_2),recipe_size)))
      if (m=='ratio.abs'){
        tmp2 <- scale(apply(tmp,1,function(x) prod(df_2[x,m])))
      }else if (m=='diff.abs'){
        tmp2 <- scale(apply(tmp,1,function(x) sum(df_2[x,m])))
      }
      tmp$agg <- tmp2
      tmp$count <- 1
      tmp <- tmp[tmp$agg>=qnorm(1-p),]
      if(nrow(tmp)>=1){
        output <- as.data.frame(table(unlist(tmp[,1:recipe_size])))
        output <- output[order(output$Freq, decreasing = T),]
        names(output) <- c("TF","freq")
        output$type <- getRecipeType(output, name_s, name_t)
        write.table(output, paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"),
                    row.names = F, quote = F)
        
      }else{print("No recipes pass the p-value cutoff")}
      return(tmp)
    }
  }
  
  if (recipe_size==1){
    candidatesL <- as.data.frame(table(mapply(function(x,y) getCandidates(idx_s = x, idx_t = y, m=m),
                                              rep(1:length(sources),each=length(targets)),
                                              rep(1:length(targets),length(sources)))))
    output <- head(candidatesL[order(candidatesL$Freq, decreasing = T),], n = 10)
    names(output) <- c("TF","freq")
    output$type <- getRecipeType(output, sources, targets)
    
  } else{
    candidatesL <- do.call(rbind,lapply(1:length(sources), function(x) do.call(rbind,lapply(1:length(targets), function(y) getCandidates(x,y,m)))))
    
    # get the TF freq summary
    output <- as.data.frame(table(unlist(candidatesL[,1:recipe_size])))
    output <- output[order(output$Freq, decreasing = T),]
    names(output) <- c("TF","freq")
    output$type <- getRecipeType(output, sources, targets)
    
    # keep those having non-zero expression values in more than 20% of the samples
    st <- c(sources, targets)
    output <- output[unlist(lapply(output$TF, function(x) sum(R[match(x,rownames(R)), st]>0) >= ceiling(0.2*length(st)))),]
    
    write.table(output, paste0("TF_freq_summary_",soN, "_to_",taN,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"),
                row.names = F, quote = F)
    
    # still needs manual modification
    if (recipe_size==2){
      output <- aggregate(cbind(agg,count) ~ V1 + V2, candidatesL, sum)
      names(output) <- c("TF1","TF2","z-score.sum","freq")
      output$TF1_type <- getRecipeType(output,sources,targets,c="TF1")
      output$TF2_type <- getRecipeType(output,sources,targets,c="TF2")
    }else if (recipe_size==3){
      output <- aggregate(cbind(agg,count) ~ V1 + V2 + V3, candidatesL, sum)
      names(output) <- c("TF1","TF2","TF3","z-score.sum","freq")
      output$TF1_type <- getRecipeType(output,sources,targets,c="TF1")
      output$TF2_type <- getRecipeType(output,sources,targets,c="TF2")
      output$TF3_type <- getRecipeType(output,sources,targets,c="TF3")
    }

    output$`z-score.sum` <- round(output$`z-score.sum`)
    output <- head(output[order(output$freq, output$`z-score.sum`, decreasing = T),], n = rank_cut)
  }
  write.table(output, paste0("recipe_summary_",soN, "_to_",taN,"_",method,"_",m,"_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"),
              row.names = F, quote = F)
}
