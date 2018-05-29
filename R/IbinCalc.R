#'To evaluate the adequacy of the set of parameters the information content per bin (Ibin)
#'was used (Arrigo et al., 2009)

#Calculates Ibin
IbinCalc<-function(matrix)
{
  Ibins<-numeric(0)
  for(i in 1:nrow(matrix))
  {
    Ibins<-append(Ibins,(mean(colSums(apply(matrix[-i,], 1,function(x) abs(x-matrix[i,]) )==1)))/ncol(matrix))
  }#for i
  return(round(median(as.matrix(Ibins)),digits = 3))
}#IbinCalc
