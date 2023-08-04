#' @export
cumulcalib <- function(y, p, method=c('BCI2p','CCI2p','BCI1p','CCI1p'), ordered=F, n_sim=0)
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

  C <- cumsum(y-p)/n
  C_n <- C[n]
  C_star <- max(abs(C))
  B_star <- max(abs(C-C[n]*t/t[n]))

  scale <- n/sqrt(T_)
  S <- scale*C
  S_star <- scale*C_star
  S_n <- scale*C_n

  for(i in 1:length(method))
  {
    mt <- method[i]

    if(mt %in% c('CCI2p','BCI2p'))
    {
      stat1 <- S[n]
      pval1 <- 2*pnorm(-abs(S[n]),0,1) #Two-sided z test

      if(mt=='BCI2p')
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
        pval2 <- 1-pMAD_BM_c(stat2, w1=S[n])
      }

      fisher <- -2*(log(pval1)+log(pval2))
      pval <- 1-pchisq(fisher,4)

      methods[[mt]]$stat <- fisher
      methods[[mt]]$pval <- pval
      methods[[mt]]$stat_by_component <- c(mean=stat1, distance=stat2)
      methods[[mt]]$pval_by_component <- c(mean=pval1, distance=pval2)
      methods[[mt]]$loc <- loc
    }

    if(mt %in% c('CCI1p','BCI1p'))
    {
      if(mt=='BCI1p')
      {
        S2 <- S-S[n]*t/t[n]
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

  out$data <- cbind(X=X,t=t,C=C, S=S)
  out$method <- names(methods[1])

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




pMAD_BM <- function(x, max_m=100)
{
  pi <- 3.1415926535897932384626
  out <-0
  m<-0:max_m
  out <- sum((-1)^m/(2*m+1)*exp(-(2*m+1)^2*pi^2/8/x^2))

  return(4/pi*out)
}

qMAD_BM <- function(q)
{
  x <- uniroot(function(x) {pMAD_BM(x)-q}, interval=c(0,10))
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



#Taken from the CPAT package
pKolmogorov <- function (q, summands=ceiling(q*sqrt(72)+3/2), N=NULL)
{
  if(!is.null(N)) q <- q*(1+0.12/sqrt(N)+0.11/N)

  sqrt(2 * pi) * sapply(q, function(x) {
    if (x > 0) {
      sum(exp(-(2 * (1:summands) - 1)^2 * pi^2/(8 * x^2)))/x
    }
    else {
      0
    }
  })
}





qKolmogorov <- function(q, summands = ceiling(q * sqrt(72) + 3/2))
{
  x <- uniroot(function(x) {pKolmogorov(x,summands)-q}, interval=c(0,10))
  unname(x$root)
}



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

qMAD_BM_c <- function(q, w1)
{
  x <- uniroot(function(x) {pMAD_BM_c(x,w1)-q}, interval=c(0,10))
  unname(x$root)
}






#' @export
plot.cumulcalib <- function(cumulcalib_obj, method=NULL, standardized=T, draw_stats=list(sigs=c(p1=0.95,p2=0.95), cols=c(p1='blue',p2='red')), xaxes=c('t','pi') , ...)
{
  args <- list(...)
  if(is.null(method))
  {
    method <- cumulcalib_obj$method
  }

  t_ <- cumulcalib_obj$data[,'t']
  X <- cumulcalib_obj$data[,'X']

  if(standardized)
  {
    W <- cumulcalib_obj$data[,'C']*cumulcalib_obj$scale
    scale <- 1
  }
  else
  {
    W <- cumulcalib_obj$data[,'C']
    scale <- 1/cumulcalib_obj$scale
  }

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
    sig_p1 <- scale*qnorm(1-(1-draw_stats$sigs['p1'])/2)

    if(method %in% c("CCI1p"))
    {
      sig_p2 <- scale*qMAD_BM(draw_stats$sigs['p2'])
    }
    if(method %in% c("CCI2p"))
    {
      sig_p2 <- scale*qMAD_BM_c(draw_stats$sigs['p2'],w1=unname(cumulcalib_obj$by_method$CCI2p$stat_by_component['mean']))
    }
    if(method %in% c("BCI1p","BCI2p"))
    {
      sig_p2 <- scale*qKolmogorov(draw_stats$sigs['p2'])
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
    args$ylim<- c(min(W,-sig_p1*(sign_p1==-1),-sig_p2*(sign_p2==-1),-scale),max(W,sig_p1*(sign_p1==1),sig_p2*(sign_p2==1),scale))
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

  args$x<-t_
  args$y<-W

  args$xlab<-""
  args$ylab<-ifelse(standardized,"Standardized cumulative sum", "Scaled cumulative sum")

  args$xaxs<-'i'

  do.call(plot, args)

  #Axes
  i <- match('t',xaxes)
  if(!is.na(i))
  {
    side<-ifelse(i==1,1,3)
    axis(side,at=(0:10)/10, labels=(0:10)/10)
    mtext('t',side, line=2)
  }
  i <- match('pi',xaxes)
  if(!is.na(i))
  {
    side<-ifelse(i==1,1,3)
    xs <- round(unlist(lapply((0:10)/10, function (x) {cumulcalib_obj$data[which.min((cumulcalib_obj$data[,'t']-x)^2),'X']})),3)
    axis(side=side, at = (0:10)/10, labels=xs , padj = 1)
    mtext(expression(pi),side,line=2)
  }

  axis(side=2, at=pretty(args$ylim))

  lines(c(0,1),c(0,0),col="grey")

  #Triangle
  {
    polygon(x=c(0,-0.03,-0.03), y=c(0,-1,1)*scale, col = 'black')
  }

  #P1 lines
  if(method %in% c("CCI2p","BCI2p") && !is.null(draw_stats$cols['p1']))
  {
    lines(c(1,1),c(0,W[n]), col=draw_stats$cols['p1'])
    lines(c(0,1),c(sign_p1*sig_p1,sign_p1*sig_p1),col=draw_stats$cols['p1'],lty=2)
  }

  #P2 lines
  if(!is.null(draw_stats$cols['p2']))
  {
    loc <- cumulcalib_obj$by_method[[method]]$loc
    if(method %in% c('BCI2p')) #If 2p bridge test then adjust the length of the red line and draw the bridge line, BUT only if the graph is standardized
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




#' @export
summary.cumulcalib <- function(cumulcalib_obj, method=NULL)
{
  if(is.null(method))
  {
    method <- cumulcalib_obj$method
  }

  n <- dim(cumulcalib_obj$data)[1]

  if(method=="CCI1p")
  {
    print("Method: One-part Brownian Motion (CCI1p)")
    print(paste("Test statistic value:",cumulcalib_obj$by_method[[method]]$stat))
    print(paste("P-value:",cumulcalib_obj$by_method[[method]]$pval))
    print(paste("Location of maximum drift:",cumulcalib_obj$by_method[[method]]$loc,
                " | time value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'t'],
                " | predictor value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'X']))
  }
  if(method=="CCI2p")
  {
    print("Method: Two-part Brownian Motion (CCI2p)")
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
  if(method=="BCI1p")
  {
    print("Method: One-part Brownian Bridge (BCI1p)")
    print(paste("Test statistic value:",cumulcalib_obj$by_method[[method]]$stat))
    print(paste("P-value:",cumulcalib_obj$by_method[[method]]$pval))
    print(paste("Location of maximum drift:",cumulcalib_obj$by_method[[method]]$loc,
                " | time value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'t'],
                " | predictor value:", cumulcalib_obj$data[cumulcalib_obj$by_method[[method]]$loc,'X']))
  }
  if(method=="BCI2p")
  {
    print("Method: Two-part Brownian Bridge (BCI2p)")
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
