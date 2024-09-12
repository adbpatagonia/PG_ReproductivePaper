#Standardize continuous covariates


standardize_x <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}