#'ParamGeno
#'
#'Exact search for "best" parameters for AFLP binary matrix composition with
#'package RawGeno V 2.0  \insertCite{arrigo2012automated}{CookedGeno}.
#'
#'Recomendations from \insertCite{arrigo2012automated}{CookedGeno} has been followed with an arbitrary discretization of the ranges as follow:
#'
#'* Maximun bin with (MaxBin): from 1.5 to 2.0 by 0.1 increase
#'
#'* Minimun bin with (MinBin): from 1.0 to 1.5 by 0.1 increase
#'
#'* Relative fluorescens units treshold (cutRFU) for the normalized intensity average:
#' from 50 to 200 by 25 increase.
#'
#'* Reproducibility filter (TRESH): 80
#'
#'It requiers RawGeno to be running and the AFLP data already loaded.
#'
#' @author Sergio E. Sclovich
#' @param AFLP: An AFLP object from Package RawGeno
#' @param case: A vector with the data section to be evaluated.
#'              i.e. if you want to evaluate parameters from 100 to 250 bp
#'              the vector should be:
#'
#'              case = c(100,250).
#'
#' @examples #Load data into RawGeno
#'           ParamGeno(AFLP,c(100,400))
#'
#' @importFrom Rdpack reprompt
#' @references \insertAllCited{}

#Individual parameter evaluation
ParamGeno<-function(AFLP,case, who = who, keep = keep, thresh = thresh)
{
  MaxW=seq(from = 1.5, to = 2, length.out = 6)
  MinW=seq(from = 1, to = 1.5, length.out = 6)
  RFU=seq(from = 50, to = 200, length.out = 7)
  output<-matrix(0,ncol = 7)
  colnames(output)<-c("MaxBinWith","MinBinWith","from","to","RFU","Ibin","#bins")

  limitinf<-case[1]
  limitsup<-case[2]

  for(i in MaxW){
    for(j in MinW){
      for(k in RFU){
        EXTRACTAFLP(all.dat=AFLP$all.dat,samples.names=AFLP$samples.names, MAXBIN=i,
                    MINBIN=j,RMIN=limitinf,RMAX=limitsup, cutRFU=k, who = who, thresh = thresh, keep = keep)
        output<-rbind(output,c(i,j,limitinf,limitsup,k,data.binary$table.stats[7,1],data.binary$table.stats[2,1]))
      }#for RFU
    }#for MinW
  }#for MaxW
  return(as.table(output[-1,]))
}#ParamGeno
