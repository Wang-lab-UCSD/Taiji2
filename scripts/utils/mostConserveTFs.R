# =======================================================================================
#  mostConserveTFs is a function to select most invariable TFs based on CV
#  Data: row is TF, column is sample
#  cv is the threshold, default is 0.5
#  return the Data with trimmed rows, row now is most invariable TFs
# =======================================================================================
mostConserveTFs <- function(Data,cv){
  co.var <- function(x) ( sd(x)/mean(x) )
  rcv <- apply(Data, MARGIN=1, FUN=co.var)
  Data <- cbind(Data, 'Co_Variance'=rcv)
  Data <- Data[order(Data$Co_Variance),]
  Data <- Data[Data$Co_Variance <= cv, ] 
  Data <- Data[,1:length(Data)-1]
  return(Data)
}
