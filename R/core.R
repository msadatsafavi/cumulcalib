#' Cumulative calibration assessment
#'
#' This is the core function for performing cumulative calibration assessment
#'
#' @return an objective of class cumulcalib that can be printed or plotted
#  @seealso [stringi::stri_length()] which this function wraps.
#' @param y vector of binary responses
#' @param p vector of predicted probabilities.
#' @param method string with either BB (Brownian bridge test, default method), BM (Brownian motion test), BM2p (two-part BM test - experimental), BB1p (one-part BB test wit only the 'bridge' component). Multiple methods can be specified. The first one will be the 'main' method (e.g., when submitting the resulting object to plot()). Default is c("BB","BM")
#' @param ordered if TRUE, y and p are already ordered based on ascending values of p. This is to speed up simulations.
#' @param n_sim if >0, indicates a simulation-based test is requested for inference.
#' @examples
#' pi <- rbeta(1000,1,2)
#' Y <- rbinom(length(pi),1,pi)
#' res <- cumulcalib(Y, pi, method="BB")
#' summary(res)
#' plot(res)
#' @export
cumulcalib <- function(y, p, method=c("BB","BM"), ordered=F, n_sim=0)
{
  out <- list()
  methods <- list()

  if(!ordered) #Order ascendingly based on p, if not already ordered
  {
    o <- order(p)
    p <- p[o]
    y <- y[o]
  }

  n <- length(p)

  #The time component of the random walk
  T_ <- sum(p*(1-p)) #Total 'time'
  t <- cumsum(p*(1-p))/T_ #time values at each p
  X <- p #Predicted values

  #Scaled cumulative sum (divided by n)
  C <- cumsum(y-p)/n  #Scaled partial sum of prediction errors
  C_n <- C[n] #Mean calibration error
  C_star <- max(abs(C))  #Macimum distance

  #This is the S process as described in the paper. The function returns all these metrics regardless of which method is used. But inference is only done per specified method(s)
  scale <- n/sqrt(T_)
  S <- scale*C
  S_star <- scale*C_star
  B_star <- max(abs(S-S[n]*t))
  S_n <- scale*C_n

  #Inference part. We loop over the different methods requested by the user
  for(i in 1:length(method))
  {
    mt <- method[i]

    if(mt %in% c('BM2p','BB')) #Two-part BM and BB both generate component-specific p-values
    {
      stat1 <- S[n]
      pval1 <- 2*pnorm(-abs(S[n]),0,1) #Two-sided z test for mean calibration

      if(mt=='BB')
      {
        loc <- which.max(S-S[n]*t)  #The bridge component of the BB test
        stat2 <- B_star
        pval2 <- 1-pKolmogorov(stat2)
      }
      else
      {
        loc <- which.max(abs(S))
        stat2 <- max(abs(S))
        pval2 <- 1-pMAD_BM_c(stat2, w1=S[n])
      }

      fisher <- -2*(log(pval1)+log(pval2))  #Fisher's method for combining p-values
      pval <- 1-pchisq(fisher,4)

      methods[[mt]]$stat <- fisher
      methods[[mt]]$pval <- pval
      methods[[mt]]$stat_by_component <- c(mean=stat1, distance=stat2)
      methods[[mt]]$pval_by_component <- c(mean=pval1, distance=pval2)
      methods[[mt]]$loc <- loc
    }

    if(mt %in% c('BM','BB1p'))  #These are one-part methods
    {
      if(mt=='BB1p')
      {
        S2 <- S-S[n]*t
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
        pval <- 1-pMAD_BM(stat)
        methods[[mt]]$stat <- stat
        methods[[mt]]$pval <- pval
        methods[[mt]]$loc <- loc
      }
    }
  }

  out$data <- cbind(X=X,t=t,C=C, S=S) #Returns the generated random-walk
  out$method <- names(methods[1]) #The base method is the first requested one

  out$scale <- scale
  out$C_n <- C_n
  out$S_n <- S_n
  out$C_star <- C_star
  out$S_star <- S_star
  out$B_star <- B_star

  #Copy the first method results to root of the list
  for(nm in names(methods[[1]]))
  {
    out[nm] <- methods[[1]][nm]
  }
  out$by_method <- methods

  class(out) <- "cumulcalib"
  return(out)
}





#' CDF of the distribution of the maximum absolute deviation of Brownian motion in [0,1] interval
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param summands maximum number of terms to be evaluated in the infinite series (default=100)
#' @export
pMAD_BM <- function(q, summands=100)
{
  pi <- base::pi
  m<-0:summands
  out <- sum((-1)^m/(2*m+1)*exp(-(2*m+1)^2*pi^2/8/q^2))

  return(4/pi*out)
}


#' Quantile function of the distribution of the maximum absolute deviation of Brownian motion in [0,1] interval
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @export
qMAD_BM <- function(p)
{
  x <- uniroot(function(x) {pMAD_BM(x)-p}, interval=c(0,10))
  unname(x$root)
}




#
# W2CDF <- function(q, n_terms=10)
# {
#   n <- 0:n_terms
#   A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
#   A2 <- erfc((4*n+1)/(2*sqrt(2*q)))
#
#   sqrt(2)*sum(A1*A2)
# }
#
#
# W2PDF <- function(x, n_terms=10)
# {
#   n <- 0:n_terms
#   A1<- gamma(-0.5+1)/(gamma(n+1)*gamma(-0.5-n+1))
#   A2 <- (4*n+1)*exp(-((4*n+1)^2)/(8*x))
#
#   1/(2*sqrt(base:::pi*x^3))*sum(A1*A2)
# }



#' CDF of the Kolmogorov distribution
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param sumands maximum number of terms to be evaluated in the infinite series (default=ceiling(q*sqrt(72)+3/2))
#' @export
pKolmogorov <- function (q, summands=ceiling(q*sqrt(72)+3/2))
{
  if(!is.null(summands)) q <- q*(1+0.12/sqrt(summands)+0.11/summands)

  sqrt(2 * pi) * sapply(q, function(x) {
    if (x > 0) {
      sum(exp(-(2 * (1:summands) - 1)^2 * pi^2/(8 * x^2)))/x
    }
    else {
      0
    }
  })
}


#' Quantile function of the Kolmogorov distribution
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @export
qKolmogorov <- function(p)
{
  x <- uniroot(function(x) {pKolmogorov(x)-p}, interval=c(0,10))
  unname(x$root)
}




#' CDF of the distribution of the maximum absolute deviation of Brownian motion in [0,1] interval, conditional on its terminal value
#' @return a scalar value
#' @param q the quantity at which CDF will be evaluated. Currently accepts only a scalar
#' @param w1 the terminal value
#' @param method different infinite series to use (1,2,3)
#' @param tolerance numerical tolerance as the stopping rule when evaluating the infinite sum (default -30 on the expotential scale)
#' @param summands number of terms to evaluate (default is ceiling(q * sqrt(72) + 3/2))
#' @export
pMAD_BM_c <- function(q, w1, method=1, exp_tolerance=-30, summands = ceiling(q * sqrt(72) + 3/2))
{
  if(q<=max(0,w1)) return(0)
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
  if(3 %in% method)
  { #Wrong
    n <- -1000:1000
    x <- sum(dnorm(w1+4*n*q)-dnorm(w1+4*n*q+2*q))
    out <- c(out,x)
  }

  out
}




#' Quantile function of the distribution of the maximum absolute deviation of Brownian motion in [0,1] interval, conditional on its terminal value
#' @return a scalar value
#' @param p the quantity at which the quantile function will be evaluated. Currently accepts only a scalar
#' @param w1 the terminal value
#' @export
qMAD_BM_c <- function(p, w1)
{
  x <- uniroot(function(x) {pMAD_BM_c(x,w1)-p}, interval=c(0,10))
  unname(x$root)
}






#' Generates cumulative calibration plots
#' @return None
#' @param cumulcalib_obj An object of class cumulcalib generated by cumulcalib()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, returns the default method used in the cumulcalib() call
#' @param draw_stats a list with two elements: sigs=significance levels (a scalar if one-part test, a vector of 2 if two-part test) and cols=colors (a character string if one-part test, a vector of two-character strings if two-part test). Default is list(sigs=c(p1=0.95,p2=0.95), cols=c(p1='blue',p2='red'))
#' @param x2axis=T If true, draws a second x-axis (on top) showing predicted risks
#' @param y2axis=T If true, draws a second y-axis (on right) showing scaled partial sums
#' @export
plot.cumulcalib <- function(cumulcalib_obj, method=NULL, draw_stats=list(sigs=c(p1=0.95,p2=0.95), cols=c(p1='blue',p2='red')), x2axis=T, y2axis=T, ...)
{
  args <- list(...)
  if(is.null(method))
  {
    method <- cumulcalib_obj$method
  }
  else
  {
    if(!("BM" %in% names(cumulcalib_obj$by_method)))
       stop("The requested method has not been requested in the original cumulcalib() call")
  }

  t_ <- cumulcalib_obj$data[,'t']
  X <- cumulcalib_obj$data[,'X']
  W <- cumulcalib_obj$data[,'S']

  n <- length(W)

  sign_p1 <- sign(W[n])
  sign_p2 <- sign(W[cumulcalib_obj$by_method[[method]]$loc])

  if(method=="BCI1p") #The only method that messes with S when drawing it.
  {
    W<-W-t_*W[n]
  }

  sig_p1 <- sig_p2 <- 0 #0 indicates do not draw
  if(!is.null(draw_stats))
  {
    sig_p1 <- qnorm(1-(1-draw_stats$sigs['p1'])/2)

    if(method %in% c("BM"))
    {
      sig_p2 <- qMAD_BM(draw_stats$sigs['p2'])
    }
    if(method %in% c("BM2p"))
    {
      sig_p2 <- qMAD_BM_c(draw_stats$sigs['p2'],w1=unname(cumulcalib_obj$by_method$BM2p$stat_by_component['mean']))
    }
    if(method %in% c("BB1p","BB"))
    {
      sig_p2 <- qKolmogorov(draw_stats$sigs['p2'])
    }
  }

  i <- match("xlim",names(args))
  if(is.na(i))
  {
    args$xlim<-c(-0.03,1.01)
  }
  i <- match("ylim",names(args))
  if(is.na(i))
  {
    args$ylim<- c(min(W,-sig_p1*(sign_p1==-1),-sig_p2*(sign_p2==-1),-1),max(W,sig_p1*(sign_p1==1),sig_p2*(sign_p2==1),1))
  }
  i <- match("type",names(args))
  if(is.na(i))
  {
    args$type<-'l'
  }
  i <- match("xaxt",names(args))
  if(is.na(i))
  {
    args$xaxt<-'n'
  }
  i <- match("yaxt",names(args))
  if(is.na(i))
  {
    args$yaxt<-'n'
  }
  i <- match("xlab",names(args))
  if(is.na(i))
  {
    args$xlab<-'Time'
  }
  i <- match("ylab",names(args))
  if(is.na(i))
  {
    args$ylab<-'Standardized cumulative sum'
  }

  args$x<-t_
  args$y<-W

  args$xaxs<-'i'

  par(mar=c(3,3,ifelse(x2axis,3,1),ifelse(y2axis,3,1)))

  do.call(plot, args)

  #Axes
  axis(1, line=0,  at=(0:10)/10, labels=(0:10)/10, padj=0)
  mtext(args$xlab, 1, line=2)

  axis(side=2, at=pretty(args$ylim), padj=0)
  mtext(args$ylab, 2, line=2)

  if(x2axis)
  {
    xs <- round(unlist(lapply((0:10)/10, function (x) {cumulcalib_obj$data[which.min((cumulcalib_obj$data[,'t']-x)^2),'X']})),3)
    axis(side=3, at=(0:10)/10, labels=xs, padj=0)
    mtext("Predicted values", 3, line=2) #expression(pi) for label will generate the symbol!
  }

  if(y2axis)
  {
    y2p <- pretty(args$ylim/cumulcalib_obj$scale)
    axis(side=4, at=y2p*cumulcalib_obj$scale, labels=y2p, padj=0)
    mtext('Scaled cumulative sum',4, line=2)
  }


  lines(c(0,1),c(0,0),col="grey")

  #Triangle
  {
    polygon(x=c(0,-0.03,-0.03), y=c(0,-1,1), col = 'black')
  }

  #P1 lines
  if(method %in% c("BM2p","BB") && !is.null(draw_stats$cols['p1']))
  {
    lines(c(1,1),c(0,W[n]), col=draw_stats$cols['p1'])
    lines(c(0,1),c(sign_p1*sig_p1,sign_p1*sig_p1),col=draw_stats$cols['p1'],lty=2)
  }

  #P2 lines
  if(!is.null(draw_stats$cols['p2']))
  {
    loc <- cumulcalib_obj$by_method[[method]]$loc
    if(method %in% c('BB')) #If 2p bridge test then adjust the length of the red line and draw the bridge line
    {
      lines(c(0,1),c(0,W[n]),col="gray", lty=2)
      lines(c(t_[loc],t_[loc]),c(t_[loc]/t_[n]*W[n],W[loc]),col=draw_stats$cols['p2'])
      lines(c(0,1),c(sign_p2*sig_p2,sign_p2*sig_p2+W[n]),col=draw_stats$cols['p2'],lty=2)
    }
    else if(method %in% c('BCI1p'))
    {
      lines(c(0,1),c(0,W[n]),col="gray")
      lines(c(t_[loc],t_[loc]),c(0,W[loc]),col=draw_stats$cols['p2'])
      lines(c(0,1),c(sign_p2*sig_p2,sign_p2*sig_p2),col=draw_stats$cols['p2'],lty=2)
    }
    else
    {
      lines(c(t_[loc],t_[loc]),c(0,W[loc]),col=draw_stats$cols['p2'])
      lines(c(0,1),c(sign_p2*sig_p2,sign_p2*sig_p2),col=draw_stats$cols['p2'],lty=2)
    }
  }
}




pKolmogorov2 <- function(x, n=100)
{
  d <- x/sqrt(n)
  X <- seq(from=0, to=1-d, length.out=n)
  1-ks.test(X,punif, exact=F)$p.value
}




#' Summarizes a cumulcalib object
#' @return None
#' @param cumulcalib_obj An object of class cumulcalib generated by cumulcalib()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, returns the default method used in the cumulcalib() call
#' @export
print.cumulcalib <- function(cumulcalib_obj, method=NULL)
{
  if(is.null(method))
  {
    method <- cumulcalib_obj$method
  }
  else
  {
    if(!("BM" %in% names(cumulcalib_obj$by_method)))
      stop("The requested method has not been requested in the original cumulcalib() call")
  }


  n <- dim(cumulcalib_obj$data)[1]
  print(paste("C_n:",cumulcalib_obj$C_n))
  print(paste("S_n:",cumulcalib_obj$S_n))
  print(paste("C_star:",cumulcalib_obj$C_star))
  print(paste("S_star:",cumulcalib_obj$S_star))
  print(paste("B_star:",cumulcalib_obj$B_star))

  if(method=="BM")
  {
    print("Method: One-part Brownian Motion (BM)")
    print(paste("Test statistic value:",cumulcalib_obj$by_method[[method]]$stat))
    print(paste("P-value:",cumulcalib_obj$by_method[[method]]$pval))
    print(paste("Location of maximum drift:",cumulcalib_obj$by_method[[method]]$loc,
                " | time value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'t'],
                " | predictor value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'X']))
  }
  if(method=="BM2p")
  {
    print("Method: Two-part Brownian Motion (BM2p)")
    print("Test statistic values:")
    print(cumulcalib_obj$by_method[[method]]$stat_by_component)
    print("Component-wise p-values:")
    print(cumulcalib_obj$by_method[[method]]$pval_by_component)
    print(paste("Combined value (Fisher'smethod):",cumulcalib_obj$by_method[[method]]$pval))
    print(paste("Terminal value (z score for mean claibration):", cumulcalib_obj$data[n,'S']))
    print(paste("Location of maximum drift:",cumulcalib_obj$by_method[[method]]$loc,
                " | time value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'t'],
                " | predictor value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'X']))
  }
  if(method=="BB1p")
  {
    print("Method: One-part Brownian Bridge (BB1p)")
    print(paste("Test statistic value:",cumulcalib_obj$by_method[[method]]$stat))
    print(paste("P-value:",cumulcalib_obj$by_method[[method]]$pval))
    print(paste("Location of maximum drift:",cumulcalib_obj$by_method[[method]]$loc,
                " | time value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'t'],
                " | predictor value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'X']))
  }
  if(method=="BB")
  {
    print("Method: Two-part Brownian Bridge (BB)")
    print("Test statistic values:")
    print(cumulcalib_obj$by_method[[method]]$stat_by_component)
    print("Component-wise p-values:")
    print(cumulcalib_obj$by_method[[method]]$pval_by_component)
    print(paste("Combined value (Fisher'smethod):",cumulcalib_obj$by_method[[method]]$pval))
    print(paste("Terminal value (z score for mean claibration):", cumulcalib_obj$data[n,'S']))
    print(paste("Location of maximum drift:",cumulcalib_obj$by_method[[method]]$loc,
                " | time value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'t'],
                " | predictor value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'X']))
  }
}


#' Prints a cumulcalib object
#' @return None
#' @param cumulcalib_obj An object of class cumulcalib generated by cumulcalib()
#' @param method Which method to use. Options are BB (Brownian bridge test), BM (Brownian motion test), BB1p (1-part Brownian bridge test), and BM2p (2-part Brownian bridge test). If unspecified, returns the default method used in the cumulcalib() call
#' @export
summary.cumulcalib <- function(cumulcalib_obj, method=NULL)
{
  print.cumulcalib(cumulcalib_obj, method)
}
