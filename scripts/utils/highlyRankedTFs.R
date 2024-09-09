# =======================================================================================
#  highlyRankedTFs is a function to select TFs which rank high, average is used
#  Data: row is TF, column is sample
#  p is the threshold for percentage, default is 10%
#  return the Data with trimmed rows, row now is most  TFs
# =======================================================================================
highlyRankedTFs <- function(Data,p){
  rsum <- apply(Data, MARGIN=1, FUN=mean)
  Data <- cbind(Data, 'Mean'=rsum)
  Data <- Data[order(Data$Mean,decreasing = 'True'),]
  active <- Data[1:as.integer(dim(Data)[1]*p),-dim(Data)[2]]
  return(active)
  
}
