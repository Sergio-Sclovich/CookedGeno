#'BestComb
#'
#'Search for a "Best" combination in a partitioned analysis of AFLP crhomatograms.
#'
#'There are several ovservations that peaks which represents bigger DNA fragments
#'in a AFLP crhomatograms had lower intensities \insertCite{questiau1999amplified,papa2005genome,trybush2006getting,whitlock2008objective,herrmann2010selection}{CookedGeno}
#'Therefore, the function BestComb() allow to test partitioned ranges. This way it is possible to search a local better fit of
#'parameters in order to retain a maximun of loci with minimum genotyping error.
#'
#'
#' @author Sergio E. Sclovich
#' @param AFLP: an AFLP profile loaded into RawGeno
#' @param mini: from where to start the analysis of the chromatogram in bp (i.e. 100)
#' @param maxi: to where analyse the cromatogram in bp (i.e. 400)
#' @param who: characters defining the replicates as define in RawGeno (i.e. "_R")
#' @param keep: whether to keep (T) or not (F, defaut) non replicable bins
#' @param thresh: threshold to evaluate reproducivitlity as implemented in RawGeno (i.e. 80, default)
#' @examples #Load data into RawGeno
#'           BestComb(AFLP,100,400)
#'
#' @importFrom Rdpack reprompt
#' @references \insertAllCited{}

#Check for the maximun Ibin
BestComb<-function(AFLP, mini, maxi, who ="_R",keep = F, thresh = 80)
{
  table<-TestComb(AFLP,mini,maxi, who = who, keep = keep, thresh = thresh)
  Total<-table[6,1]
  TotalIbin<-numeric(0)
  TotalNBin <- numeric(0)

  for(case in c(2,4,6))
  {
    for(half in 0:1)
    {
      EXTRACTAFLP(all.dat=AFLP$all.dat,samples.names=AFLP$samples.names, MAXBIN=table[1,case+half],
                  MINBIN=table[2,case+half],RMIN=table[3,case+half],RMAX=table[4,case+half], cutRFU=table[5,case+half], who= who, thresh= thresh, keep = keep)
      if(half==0)#first half saved
      {
        mxbinary<-t(data.binary$data.binary)
      }#if half==0
      else#second half saved for those individuals that has their range long enough
      {
        mxbinary<-cbind(mxbinary[intersect(rownames(t(data.binary$data.binary)),rownames(mxbinary)),],t(data.binary$data.binary)[intersect(rownames(t(data.binary$data.binary)),rownames(mxbinary)),])
      }#else half==0
    }#for half
    TotalIbin<- append(TotalIbin, IbinCalc(mxbinary))
  }#for case
  aux<-2
  TotalNBin <- table[7,1]

  for(i in 2:7)
  {
    if(i%%2 == 1)#if uneven
    {
      Total[i]<- TotalIbin[i-aux]
      TotalNBin[i] <- table[7,i-1]+table[7,i]
      aux<-aux+1
    }#if even
    else
    {
      Total[i]<- "-"
      TotalNBin[i]<- "-"
    }#else
  }#for disposition
  table <- rbind(table,TotalNBin)
  table <- rbind(table,Total)
  rownames(table)[8]<-"Total#bins"
  rownames(table)[9]<-"TotalIbin"
  return(table)
}#BestComb

generate.partitions <- function(mini, maxi)
{
  ranges = matrix(c(mini, maxi), nrow = 2)
  part = round((maxi-mini)/4) #define lenght of the partition

  ####declarations####

  for(i in 1:3)
  {
    ranges <- cbind(ranges,c(mini,mini+(part*i)))
    ranges <- cbind(ranges,c(mini+(part*i),maxi))
  }#for
  return(ranges)
}#generate.cases

TestComb<-function(AFLP,mini,maxi, who = who, keep = keep, thresh = thresh) #dif between maxi and mini can not be less than 100
{
  ##Declarations###
  output1<-matrix(0, nrow = 7)
  ranges = generate.partitions(mini,maxi)
  #cases =dim(ranges)[2]

  ####evaluate cases#######
  for(i in 1:dim(ranges)[2])
  {
    aux<-ParamGeno(AFLP,ranges[,i], who, keep, thresh)
    print(i)

    if(is.null(dim(aux[which(as.numeric(aux[,6])==max(as.numeric(aux[,6]))),])))#check if there is more than one case with the same maximum
    {
      output1<-cbind(output1,aux[which(as.numeric(aux[,6])==max(as.numeric(aux[,6]))),])
    }#if
    else
    {
      aux1<-aux[which(as.numeric(aux[,6])==max(as.numeric(aux[,6]))),]
      if(all(aux1[,7]==aux1[1,7])==FALSE)
      {
        output1<-cbind(output1,aux1[which(aux1[,7]==max(aux1[,7])),])#check which has more bins
      }#if
      else#Random choice
      {
        output1<-cbind(output1,aux1[sample(nrow(aux1),1),])
      }#else
    }#else
  }#for cases

  ####organizing the data####
  output1<-output1[,-1]
  colnames(output1)<-c("total","1.1","1.2","2.1","2.2","3.1","3.2")
  rownames(output1)<-c("MaxBinWith","MinBinWith","from","to","RFU","Ibin","#bins")
  return(output1)
}#TestComb




