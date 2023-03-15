#' @export
cumulcalib <- function(y, p, method=c('twopart','onepart','twopartw','onepartw'), ordered=F, n_sim=0)
{
  out <- list()

  details <- list()

  if(!ordered)
  {
    o <- order(p)
    p <- p[o]
    y <- y[o]
  }

  n <- length(p)

  for(i in 1:length(method))
  {
    mt <- method[i]

    if(mt %in% c('twopart', 'onepart'))
    {
      t <- sum(p*(1-p)) #Total 'time'
      X <- cumsum(p*(1-p))/t
      Y <- cumsum(y-p)/sqrt(t)
    }
    if(mt %in% c('twopartw', 'onepartw'))
    {
      t <- n
      X <- (1:n)/n
      Y <- cumsum((y-p)/(sqrt(p*(1-p))))/sqrt(t)
    }

    if(mt %in% c('twopart','twopartw'))
    {
      stat1 <- Y[n]
      pval1 <- 2*pnorm(-abs(Y[n]),0,1) #Two-sided z test

      Y <- Y- Y[n]*X/X[n]

      loc <- which.max(abs(Y))
      stat2 <- max(abs(Y))
      pval2 <- 1-pKolmogorov(stat2)

      stat <- c(mean=stat1, distance=stat2)
      details[[mt]]$stat <- stat

      fisher <- -2*(log(pval1)+log(pval2))
      pval <- 1-pchisq(fisher,4)
      details[[mt]]$pval <- c(combined=pval, mean=pval1, distnce=pval2)
    }

    if(mt %in% c('onepart','onepartw'))
    {
      stat <- max(abs(Y))
      loc <- which.max(abs(Y))
      pval <- 1-erdos(stat)
      details[[mt]]$stat <- stat
      details[[mt]]$pval <- pval
    }

    if(i==1) #First method is the default and gets to return the X and Y
    {
      out$method <- mt
      out$X <- X
      out$Y <- Y
      out$stat <- stat
      out$pval <- pval
      out$loc <- loc
    }
  }

  out$details <- details
  class(out) <- "cumulcalib"
  return(out)
}






erdos <- function(x, max_m=100)
{
  pi <- 3.1415926535897932384626
  out <-0
  m<-0:max_m
  out <- sum((-1)^m/(2*m+1)*exp(-(2*m+1)^2*pi^2/8/x^2))

  return(4/pi*out)
}










W2CDF <- function(q, n_terms=10)
{
  n <- 0:n_terms
  A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
  A2 <- erfc((4*n+1)/(2*sqrt(2*q)))

  sqrt(2)*sum(A1*A2)
}


W2PDF <- function(x, n_terms=10)
{
  n <- 0:n_terms
  A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
  A2 <- (4*n+1)*exp(-((4*n+1)^2)/(8*x))

  1/(2*sqrt(base:::pi*x^3))*sum(A1*A2)
}



#Taken from the CPAT package
pKolmogorov <- function (q, summands = ceiling(q * sqrt(72) + 3/2))
{
  sqrt(2 * pi) * sapply(q, function(x) {
    if (x > 0) {
      sum(exp(-(2 * (1:summands) - 1)^2 * pi^2/(8 * x^2)))/x
    }
    else {
      0
    }
  })
}




#' @export
plot.cumulcalib <- function(cumulcalib_obj,...)
{
  if(cumulcalib_obj$method %in% c('twopart','twopartw'))
  {
    n <- length(cumulcalib_obj$Y)
    l <- cumulcalib_obj$stat[1]*(1:n)/n
    Y <- cumulcalib_obj$Y + l
    plot(cumulcalib_obj$X,Y,type='l',xlim=c(0,1),ylim=c(min(Y),max(Y)),...)
    lines(c(0,1),c(0,0),col="grey")
    lines(c(0,1),c(0,Y[n]),col="grey")
    lines(c(1,1),c(0,Y[n]),col="blue")
    loc <- cumulcalib_obj$loc
    lines(c(cumulcalib_obj$X[loc],cumulcalib_obj$X[loc]),c(l[loc],Y[loc]),col="red")
  }
  else
  {
    plot(cumulcalib_obj$X,cumulcalib_obj$Y,type='l',xlim=c(0,1),ylim=c(min(cumulcalib_obj$Y),max(cumulcalib_obj$Y)),...)
    lines(c(0,1),c(0,0),col="grey")
    loc <- cumulcalib_obj$loc
    lines(c(cumulcalib_obj$X[loc],cumulcalib_obj$X[loc]),c(0,Y[loc]),col="red")
  }
}
