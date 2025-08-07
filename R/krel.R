## ADB
## relcondition function
## this function estimates a growth curve from a data set of lengths and weights and then calculates relative condition,
## based on Le Cren (1951) Kr = W/Wˆ, where Wˆ is the predicted body weight from the length-weight relationship
## Le Cren ED (1951) The length-weight relationship and seasonal cycle in gonad weight and condition in the perch (Perca fluviatilis). J Anim Ecol 20:201-219

krel <- function(length, weight) {
  lmlenweight <- lm(log(weight) ~ log(length))
  krel <- weight / (exp(coef(lmlenweight)[1]) * ((length)^coef(lmlenweight)[2]))
  return(krel)
}


Eweight <- function(length, weight) {
  len <- floor(min(length, na.rm = TRUE)):ceiling(max(length, na.rm = TRUE))
  lmlenweight <- lm(log(weight) ~ log(length))
  eweight <- (exp(coef(lmlenweight)[1]) * ((len)^coef(lmlenweight)[2]))
  return(eweight)
}

bpar <- function(length, weight) {
  len <- floor(min(length, na.rm = TRUE)):ceiling(max(length, na.rm = TRUE))
  lmlenweight <- lm(log(weight) ~ log(length))
  bpar <- coef(lmlenweight)[2]
  return(bpar)
}


Expweight <- function(length, weight) {
  # len <- floor(min(length, na.rm = TRUE)):ceiling(max(length, na.rm = TRUE))
  lmlenweight <- lm(log(weight) ~ log(length))
  eweight <- (exp(coef(lmlenweight)[1]) * ((length)^coef(lmlenweight)[2]))
  return(eweight)
}
