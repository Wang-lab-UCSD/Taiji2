# callSpecificTFs function is to generate specific TF list along with the bubble plot------------------
# bubble plot: p-value is increasing from top to bottom

# --------- arguments: ---------------------------------------------------
#   t: target group name;
#   p: pValue cut-off
#   d: log2 fold change cut-off
#   cv: most conserved TFs in target group with covariance smaller than cv
#   per: TFs whose prs are all ranked top per percent in target group
#   f: "single" or "double", whether the log2FC is single or double ended
#   sort.col: logical, whether to separate target and background visually
#   name.prefix: add prefix in output *txt file
#   normalize.method: c("default","0-1","-2-2")

# ---------main body------------------------------------------------------
library(grid)
callSpecificTFs <- function(t,pr=Data,rna=rna_seq,
                            g=group_sorted,subset="Group",
                            p=0.001,d=1,cv=0.5,per=0.5,
                            f="single",normalize.method="default",
                            is.uniform=F, is.highly=F,
                            sort.col=F, only.target=F,show.legend=T,
                            name.prefix = "",max.row = 30, angle=90,
                            label.color = brewer.pal(10, "Set3"),
                            split.width=c(1.5,10,1), split.height=c(1.5,10,1),
                            plot.margin=c(4,4.5,0.5,1), bubble.type=1, fl.type="pdf"){
  
  #### bubble plot function ####
  getBubblePlot <- function(Tcell){
    
    # check the rows of bubble plot
    Tcell <- Tcell[1:min(length(Tcell), max.row)]
    
    # reverse the order of Tcell to make the top rows as most significant
    Tcell <- rev(Tcell)
    
    # check if putting T cell columns together
    if (sort.col){
      m <- setdiff(c(1:length(pr)),n)
      pr = pr[c(Tcell),c(n,m)]
      log_exp = log2(rna+1)[c(Tcell),c(n,m)]
    }else{
      pr = pr[c(Tcell),]
      log_exp = log2(rna+1)[c(Tcell),]
    }
    
    # normalize pagerank, different strategies
    if (normalize.method=="default"){
      rank <- as.data.frame(t(scale(t(as.matrix(pr)))))
    }else if (normalize.method=="0-1"){
      rank <- as.data.frame(t(apply(as.matrix(pr), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
    }else if (normalize.method=="-2-2"){
      tmp <- as.data.frame(t(scale(t(as.matrix(pr)))))
      tmp <- apply(tmp,1:2, function(x) min(x, 2))
      rank <- as.data.frame(apply(tmp,1:2, function(x) max(x, -2)))
    }


    
    # # switch the element position
    # switch <- function(x,i,j) {x[c(i,j)] <- x[c(j,i)]; x} 
    # Tcell <- switch(Tcell, 11, 13)
    
    if (only.target){
      plot <- bubble(Tcell, log_exp[,n], rank[,n], 
                     plot.margin = plot.margin, 
                     angle = angle)
      w <- ifelse(length(n)<=5, 3.5, max(5, 8 / 25 * length(n)))
      h <- max(3.5, 7.5 / 25 * length(Tcell))
    }else{
      if (bubble.type==1){
        plot <- bubble(Tcell, log_exp, rank, 
                       plot.margin = plot.margin, 
                       angle = angle)
        w <- ifelse(ncol(pr)<=5, 3.5, max(5, 8 / 25 * ncol(pr)))
        h <- max(3.5, 7.5 / 25 * length(Tcell))
      }else if (bubble.type==2){
        plot <- bubble_v2(Tcell,log_exp,rank,g,label.color,subset,split.width,split.height,plot.margin)
        w <- max(3.5, 8 / 25 * ncol(pr))
        h <- max(3.5, 9 / 25 * length(Tcell))
      }else if (bubble.type==3){
        plot <- bubble_v3(Tcell,log_exp,rank,g,label.color,subset,plot.margin)
        w <- max(3.5, 8 / 25 * ncol(pr))
        h <- max(3.5, 9 / 25 * length(Tcell))
      }
      
    }
    
    if (fl.type=="png"){
      png(gsub(" ", "", paste(name.prefix, "_",gsub("/",".",paste0(t, collapse = "_")),"_P",gsub("\\.","",p),"_F",gsub("\\.","",d),"_L",length(Tcell),"_",f, ".png")),
          units="in", width=w, height=h, res=300)
      grid.draw(plot)
      dev.off()
    }else if (fl.type=="pdf"){
      pdf(gsub(" ", "", paste(name.prefix, "_",gsub("/",".",paste0(t, collapse = "_")),"_P",gsub("\\.","",p),"_F",gsub("\\.","",d),"_L",length(Tcell),"_",f, ".pdf")),
          width=w, height=h)
      grid.draw(plot)
      dev.off()
    }

  }
  
  #### check if file exists ####
  fl <- list.files(path = "./", 
                  pattern = gsub(" ","",paste(name.prefix, "_",gsub("/",".",paste0(t, collapse = "_")),"_P",gsub("\\.","",p),"_F",gsub("\\.","",d),"_L.*_",f,".txt")))
  if (length(fl)==1){
    message("File exists!")
    Tcell <- readLines(fl)
    getBubblePlot(Tcell)
    
  }else if (length(fl)==0){
    
    #### unpaired t-test -----------------------
    if (length(t)==1){n <- which(g[[subset]]==t)}else{n <- which(g[[subset]]%in%t)}
    
    x <- unpairedtTest(pr,n,p=p,d=d,f=f)
     
    # pr in target group are high (top "per") and uniformly distributed (<cv)-----
    # every sample in the target group is rank top 10% (per=0.1)
    # if that's the case, maybe we don't need the "cv" parameter
    if (is.uniform && is.highly){
      
      x2 <- mostConserveTFs(pr[,n],cv)
      Tcell <- intersect(rownames(x),rownames(x2))
      x3 <- as.data.frame(t(apply(as.matrix(pr[Tcell,]), MARGIN = 1, 
                                    FUN = function(X) (X - min(X))/diff(range(X)))))
      # per.abs <- quantile(as.numeric(unlist(x3)),1-per)
      per.abs <- apply(x3, quantile, probs=1-per, MARGIN = 1)
      x4 <- x3[,n]
      x4 <- cbind(x4,per.abs)
      x4[] <- lapply(x4, function(x) ifelse(x>=per.abs, x, NA))
      Tcell <- rownames(x4[complete.cases(x4),])
    }else if (is.uniform && !is.highly){
      x2 <- mostConserveTFs(pr[,n],cv)
      Tcell <- intersect(rownames(x),rownames(x2))
    }else if (!is.uniform && is.highly){
      x3 <- as.data.frame(t(apply(as.matrix(pr), MARGIN = 1, 
                                  FUN = function(X) (X - min(X))/diff(range(X)))))
      per.abs <- quantile(as.numeric(unlist(x3)),1-per)
      x4 <- x3[,n]
      x4[] <- lapply(x4, function(x) ifelse(x>per.abs, x, NA))
      Tcell <- rownames(x4[complete.cases(x4),])
    }else{
      Tcell <- rownames(x)
    }
    
    if (length(Tcell) == 0){
      return()
    }else{
      
      # keep those having non-zero expression values in more than 20% of the cells in target group
      Tcell <- Tcell[unlist(lapply(Tcell, function(x) sum(rna[x, n]>0) >= ceiling(0.2*length(n))))]
      
      # if log2FC is double ended, bubble plot and txt file is arranged in increasing log2FC
      if (f != "single"){
        x <- x[order(x$'log2 fold change'),]
        Tcell <- rownames(x)
        Tcell <- Tcell[unlist(lapply(Tcell, function(x) sum(rna[x, n]>0) >= ceiling(0.2*length(n))))]
      }
      writeLines(Tcell, gsub(" ","",paste(name.prefix, "_",gsub("/",".",paste0(t, collapse = "_")),"_P",gsub("\\.","",p),"_F",gsub("\\.","",d),"_L",length(Tcell),"_",f,".txt")))
      write.table(x[rownames(x)%in%Tcell,],gsub(" ","",paste(name.prefix, "_",gsub("/",".",paste0(t, collapse = "_")),"_P",gsub("\\.","",p),"_F",gsub("\\.","",d),"_L",length(Tcell),"_",f,"_table.txt")),
                  quote = F, sep = "\t")
      
    }
    
    getBubblePlot(Tcell)
    
  }else{
    message("Multiple files exist! Will update all")
    Tcells <- lapply(fl, readLines)
    lapply(Tcells, getBubblePlot)
    
  }
  


}





