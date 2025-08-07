## ADB
##  L-W plotting function
## Sep 24, 2018

# This function takes observed fish length and weight, fits a model of the form
# weight = a*length^b
# and plots it, together with the equation and r-squared

# inputs to this function are:
# length = length of fish in mm
# weight = weight of fish in g
# label = EcoDAT Fish ID. This will be outputed when the Rmarkdown is compiled

# dependencies:
# mround.r

# output is a plot of the fitted model

plot_lw <- function(length, weight, label){
  # create data frame
  dat <- data.frame(length, weight, label)
  names(dat) <- c('length', 'weight', 'label')
  # obtain plot limits
  xylims <- dat %>%
    summarize(xmin =  mround(min(length, na.rm = T), 10), xmax = mround(max(length, na.rm = T), 10),
              ymin =  mround(min(weight, na.rm = T), 1), ymax = mround(max(weight, na.rm = T), 1)
    )
  # fit l-w linear model and get parameters
  lmlenweight <- lm(log(weight)~log(length), data = dat)
  a <- exp(coef(lmlenweight)[1])
  b <- coef(lmlenweight)[2]
  # obtain expected values
  expected <- data.frame(x = xylims$xmin:xylims$xmax, y = exp(coef(lmlenweight)[1])* ((xylims$xmin:xylims$xmax)^ coef(lmlenweight)[2]))
  # equation label
  eq <- as.character(as.expression(
    substitute(italic(y) == a %.% italic(x)^b,
               list(a = format(unname(a), digits = 2),
                    b = format(unname(b), digits = 2)))))
  # r squared label
  eqr <-  as.character(as.expression(
    substitute(italic(R)^2~"="~r2,
               list(r2 = format(summary(lmlenweight)$r.squared, digits = 2)))))
  # n  label
  eqn <-  as.character(as.expression(
    substitute(italic(N)~"="~n,
               list(n = nrow(dat)))))
  # plot
  p <-  ggplot(data = dat, aes(x = length, y = weight)) +
    geom_point(alpha = 0.4, fill = 'grey', aes(text = sprintf("Ecodat Fish ID:%s", label))) +
    geom_line(data = expected, aes(x = x, y = y), col = "red") +
    scale_x_continuous(limits = c(xylims$xmin, xylims$xmax)) +
    scale_y_continuous(limits = c(xylims$ymin, xylims$ymax)) +
    ylab('WEIGHT (g)') + xlab("FORK LENGTH (mm)") + ggtitle("") +
    annotate("text", x = -Inf, y = Inf, label = eq, hjust = -0.3, vjust = 1.5, parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = eqr, hjust = -0.5, vjust = 3, parse = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = eqn, hjust = -0.5, vjust = 6, parse = TRUE)
  # return plot
  return(p)
}
