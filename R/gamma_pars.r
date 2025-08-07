# https://math.stackexchange.com/questions/1810257/gamma-functions-mean-and-standard-deviation-through-shape-and-rate

# If X is gamma distributed with shape a and rate b, then the mean of X is
# μ=E[X]=a/b,
# and the standard deviation is
# σ = sqrt(Var[X]) =  sqrt(a)/b
# Note that a and b must be positive.
# 
# It follows from the above that, given a desired mean μ and standard deviation σ, 
# the shape and rate that produce a gamma distribution with that desired μ and σ are:
# 
# a=(μ/σ)^2
# b=μ/σ^2

# Input: mean and SD
# Output: shape and rate paramters of the gamma distribution
# this function is useful to obtain the input parameters of rgamma() to draw random numbers from a gamma distribution

gamma_pars <- function(mean, sd){
  
  shape_par <- (mean/sd)^2
  rate_par <- mean/(sd^2)
  return(list("shape" = shape_par, "rate" = rate_par))
}
