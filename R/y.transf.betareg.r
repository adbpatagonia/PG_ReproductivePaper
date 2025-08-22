#' Transformation function for beta regression data.
#' If data includes 0 or 1, this function transforms the data to the open standard unit interval (0, 1)
#' 
#' @param y Response variable as a vector
#' @return Response variable as a vector, transformed to the open standard unit interval (0, 1)
#' @examples
#' y <- c(runif(8, .01, .99), 0, 1)
#' y.transf.betareg(y)

y.transf.betareg <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

#' Smithson and Verkuilen 2006
#' https://stackoverflow.com/questions/26385617/proportion-modeling-betareg-errors
#' 