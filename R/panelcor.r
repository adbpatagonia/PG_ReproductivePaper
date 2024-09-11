  ### FUNCTIONS#######
 # function to create the correlation matrix plot (taken from http://addictedtor.free.fr/graphiques/graphcode.php?graph=137)
panelcor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r<-cor.test(x, y,method="pearson")
    r <- abs(as.double(r$estimate))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste("r(pe) = ",prefix, txt, sep="")
    r2<-cor.test(x, y,method="spearman")
    r2 <- abs(as.double(r2$estimate))
    txt2<- format(c(r2, 0.123456789), digits=digits)[1]
    txt2 <- paste("r(sp) = ",prefix, txt2, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt2)

    test <- cor.test(x,y,method="pearson")
    test2 <- cor.test(x,y,method="spearman")
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
   pvalue<-round(test$p.value,3)
   pvalue<-paste("p[pe] = ",pvalue)
   pvalue2<-round(test2$p.value,3)
   pvalue2<-paste("p[sp] = ",pvalue2)
    text(0.5, 0.9, txt)
    color<- ifelse(test$p.value<0.05,"red","black")
    text(.5, .7, pvalue,col=color )
    text(0.5, 0.3, txt2)
    color2<- ifelse(test2$p.value<0.05,"red","black")
    text(.5, .1, pvalue2,col=color2 )
  }
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
}
  ### END FUNCTIONS#######