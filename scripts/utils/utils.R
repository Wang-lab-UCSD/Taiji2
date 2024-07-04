# make the first letter uppercase and other letters as lower case
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# create new sub-dir and set the new sub-dir as the working dir
createDir <- function(subdir){
    dir.create(file.path(maindir, subdir), showWarnings = FALSE)
    setwd(file.path(maindir, subdir))
    print(getwd())
}

# normalize data across rows to make sure the values is in range [0,1] after norm
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaleData <- function(df){
  return(as.data.frame(t(apply(as.matrix(df), MARGIN = 1, FUN = range01))))
}
                               
# get z-scores of data 
zscore <- function(df){
  return(as.data.frame(t(scale(t(as.matrix(df))))))
}                      
                           
## calculate covariance, x is the matrix
co.var <- function(x) { matrixStats::rowSds(x) / matrixStats::rowMeans2(x) }

