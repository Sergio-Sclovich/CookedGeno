#'MDS
#'
#'It simply returns a pcoa object based on the Gower \insertCite{gower1971general}{CookedGeno} or Jaccard distance matrix.
#'
#'It dependens in packages APE \insertCite{paradis2004ape}{CookedGeno}, FD \insertCite{laliberte2014package}{CookedGeno} and vegan \insertCite{vegano}{CookedGeno}.
#'
#'option ordinal sets option ord of gowdis() function from package FD:
#'
#' * "classic" (default) treats ordinal variables as continuous variables.
#'
#' * "podani" refers to Eqs. 2a-b of \insertCite{podani1999extending;textual}{CookedGeno}
#'
#' * "metric" refers to Eq. 3 ibid.
#'
#' @param filename: data matrix
#' @param var_rang: which columns from the data matrix to be used (i.e. 1:20)
#' @param dist: type of distance to be used (options Gower "G" or Jaccard "J")
#' @param ordinal: method to be used for ordinal variables (see details).
#' @importFrom Rdpack reprompt
#' @references \insertAllCited{}



MDS <- function(filename, var_rang, dist, ordinal="classic") {
  require(ape)
  require(FD)
  #Choose which distance to be used
  switch(dist,
         G={
           mx_dist<-gowdis(filename[,var_rang], ord = ordinal)#calculation of Gower's distance matrix
         },
         J={
           mx_dist<-vegdist(filename[,var_rang], method = "jaccard",na.rm = T)#calculation of Jaccard's distance matrix
         }
  )
  return(pcoa(mx_dist))
}#MDS

