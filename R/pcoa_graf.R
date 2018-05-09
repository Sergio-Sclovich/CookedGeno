#'pcoa_graf
#'
#'make a 3D plot based on scatterplot3d \insertCite{ligges2002scatterplot3d}{CookedGeno}
#'of a pcoa object from package APE  \insertCite{paradis2004ape}{CookedGeno}
#'
#'Colors and simbols for the objects are chosen randomly
#'
#' @param filename: a pcoa object
#' @param angle: angle for the graph to be plot
#' @param x,y,z: Which eigenvector will be plot in each axis
#' @param nomb: Vector of names
#' @examples pcoa_graf(pcoa_object, 120, 1, 2, 3, names(Data))
#' @importFrom Rdpack reprompt
#' @references \insertAllCited{}

pcoa_graf<-function(filename, angle, x,y,z,nomb){
  require(scatterplot3d)
  pch <- round(runif(nrow(table(nomb)),1,18))
  scatterplot3d(filename$vectors[,x],
                filename$vectors[,y],
                filename$vectors[,z],
                type = "p",
                xlab = sprintf("axis %s",x),
                ylab = sprintf("axis %s",y),
                zlab = sprintf("axis %s",z),
                main = sprintf("PCoA_%s",deparse(substitute(filename))),
                color = rainbow(nrow(table(nomb)))[as.factor(nomb)],
                pch = pch[as.factor(nomb)],
                angle = angle,
                cex.axis = .5
  )
  lista_nombres<-aggregate(nomb, list(nomb), FUN=list)
  legend("topright",legend = sort(lista_nombres$Group.1), col = rainbow(nrow(table(nomb))), pch = pch, cex = .5)
}#pcoa_graf
