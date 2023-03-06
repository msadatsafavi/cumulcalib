#' @export
cumulcalib <- function(y, p, inference=T, ordered=F, n_sim=0)
{
  out <- list()

  if(!ordered)
  {
    o <- order(p)
    p <- p[o]
    y <- y[o]
  }
  n <- length(p)

  csya <- cumsum((y-p)/sqrt(p*(1-p)))/sqrt(n)
  out$x <- p
  out$y <- csya

  if(inference)
  {
    loc <- which.max(abs(csya))
    stat <- max(abs(csya))
    stat2 <- mean(csya^2)
    out$inference$stat <- stat
    out$inference$loc <- p[loc]
    out$inference$stat2 <- stat2
    if(n_sim==0)
    {
      out$inference$p_val <- 1-erdos(stat)
      out$inference$p_val2 <- 1-W2CDF(stat2)
    }
    else
    {
      ps <- rep(p,n_sim)
      ys <- rbinom(length(ps),1,ps)
      ys <- (ys-ps)/sqrt(ps*(1-ps))
      dim(ys) <- c(n,n_sim)
      csys <- apply(ys,2,cumsum)/n
      stats <- apply(abs(csys), 2, max)
      stats2 <- apply(csys^2, 2, mean)
      out$inference$p_val <- 1-ecdf(stats)(stat)
      out$inference$p_val2 <- 1-ecdf(stats2)(stat2)
    }
  }
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




c_plot <- function(y, p, type="loess", col="black", replace=T)
{
  o <- order(p)
  p <- p[o]
  y <- y[o]

  if(type=="loess")
  {
    calib_model <- loess(y~p)
    yy <- predict(calib_model)
    if(replace)
    {
      plot(p,yy,type='l',col=col,xlim=c(0,1),ylim=c(0,max(1,yy)), xlab="Predictied probabilities", ylab="True probabilities")
      lines(c(0,1),c(0,1),type='l',col='grey')
    }
    else
    {
      lines(p,yy,type='l',col=col)
    }
  }
}




erdos <- function(x, max_m=100)
{
  pi <- 3.1415926535897932384626
  out <-0
  m<-0:max_m
  out <- sum((-1)^m/(2*m+1)*exp(-(2*m+1)^2*pi^2/8/x^2))

  return(4/pi*out)
}


W2CDF <- function(x, n_terms=10)
{
  n <- 0:n_terms
  A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
  A2 <- erfc((4*n+1)/(2*sqrt(2*x)))

  sqrt(2)*sum(A1*A2)
}


W2PDF <- function(x, n_terms=10)
{
  n <- 0:n_terms
  A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
  A2 <- (4*n+1)*exp(-((4*n+1)^2)/(8*x))

  1/(2*sqrt(base:::pi*x^3))*sum(A1*A2)
}



#' @export
plot.cumulcalib <- function(cumulcalib_obj,...)
{
  plot(cumulcalib_obj$x,cumulcalib_obj$y,type='l',xlim=c(0,1),ylim=c(min(cumulcalib_obj$y),max(cumulcalib_obj$y)),...)
  lines(c(0,1),c(0,0),col="grey")
}
