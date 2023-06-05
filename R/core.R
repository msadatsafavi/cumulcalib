#' @export
cumulcalib <- function(y, p, method=c('twopart','onepart','twopart-bridge', 'onepart-bridge'), ordered=F, n_sim=0)
{
  out <- list()

  methods <- list()

  if(!ordered)
  {
    o <- order(p)
    p <- p[o]
    y <- y[o]
  }

  n <- length(p)

  T_ <- sum(p*(1-p)) #Total 'time'
  t <- cumsum(p*(1-p))/T_
  X <- p
  S <- cumsum(y-p)/sqrt(T_)

  for(i in 1:length(method))
  {
    mt <- method[i]

    if(mt %in% c('twopart','twopart-bridge'))
    {
      stat1 <- S[n]
      pval1 <- 2*pnorm(-abs(S[n]),0,1) #Two-sided z test

      if(mt=='twopart-bridge')
      {
        S2 <- S-S[n]*t/t[n]
        loc <- which.max(abs(S2))
        stat2 <- max(abs(S2))
        pval2 <- 1-pKolmogorov(stat2)
      }
      else
      {
        loc <- which.max(abs(S))
        stat2 <- max(abs(S))
        pval2 <- 1-pAbsuluteDistanceConditional(stat2, w1=S[n])
      }

      fisher <- -2*(log(pval1)+log(pval2))
      pval <- 1-pchisq(fisher,4)

      methods[[mt]]$stat <- fisher
      methods[[mt]]$pval <- pval
      methods[[mt]]$stat_by_component <- c(mean=stat1, distance=stat2)
      methods[[mt]]$pval_by_component <- c(mean=pval1, distnce=pval2)
      methods[[mt]]$loc <- loc
    }

    if(mt %in% c('onepart','onepart-bridge'))
    {
      if(mt=='onepart-bridge')
      {
        S2 <- S- S[n]*t/t[n]
        loc <- which.max(abs(S2))
        stat <- max(abs(S2))
        pval <- 1-pKolmogorov(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      }
      else
      {
        stat <- max(abs(S))
        loc <- which.max(abs(S))
        pval <- 1-erdos(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      }
    }
  }

  out$data <- cbind(X=X,t=t,S=S)
  out$method <- names(methods[1])
  #Copy the first method results to root of the list
  for(nm in names(methods[[1]]))
  {
    out[nm] <- methods[[1]][nm]
  }
  out$by_method <- methods

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




pAbsuluteDistanceConditional <- function(q, w1, method=2, exp_tolerance=-30, summands = ceiling(q * sqrt(72) + 3/2))
{
  out <- c()
  if(1 %in% method)
  {
    A <- -2*q*q
    B <- 2*q*w1
    C <- -exp_tolerance
    D <- sqrt(B*B-4*A*C)

    n <- seq(round((-B+D)/A/2,0), round((-B-D)/A/2,0))
    un <- 2*n*q

    terms <- un*w1-un^2/2
    out <- c(out,sum((-1)^n*exp(terms)))
  }
  if(2 %in% method)
  {
    out <- c(out,sqrt(2*pi)/q*sum(exp(w1^2/2-(2*(1:summands)-1)^2*pi^2/(8*q^2))*(cos((2*(1:summands)-1)*pi*w1/2/q))))
  }
  out
}




#' @export
plot.cumulcalib <- function(cumulcalib_obj, method=NULL, standardized=T, draw_stats=T, ...)
{
  if(is.null(method))
  {
    method <- cumulcalib_obj$method
  }

  S <- cumulcalib_obj$data[,'S']
  n <- length(S)

  #If not standardized, then use predictions instead of times, and also make Y divided by n rather than sqrt(T)
  if(standardized)
  {
    X <- cumulcalib_obj$data[,'t']
  }else{
    X <- cumulcalib_obj$data[,'X']
    S <- S*sqrt(sum(cumulcalib_obj$data[,'t']))/n
  }


  plot(X,S,type='l',xlim=c(0,1),ylim=c(min(S),max(S)),...)
  lines(c(0,1),c(0,0),col="grey")
  if(draw_stats)
  {
    if(method %in% c('twopart','twopart-bridge')) #Mean line only for two-part tests
    {
      lines(c(1,1),c(0,S[n]),col="blue")
    }

    loc <- cumulcalib_obj$loc
    if(method %in% c('onepart-bridge','twopart-bridge')) #If bridge test then adjust the length of the red line and draw the bridge line, BUT only if the graph is standardized
    {
      if(standardized)
      {
        lines(c(0,1),c(0,S[n]),col="gray")
        lines(c(X[loc],X[loc]),c(loc/n*S[n],S[loc]),col="red") #TODO
      }else{
        warning("Bridge statistic cannot be drawn because a non-standardized curve was requested")
      }
    }
    else
    {
      lines(c(X[loc],X[loc]),c(0,S[loc]),col="red")
    }
  }
}
