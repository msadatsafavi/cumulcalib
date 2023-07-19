#We use birth weight data in the MASS package and an exemplary #model for predicting #low birth weight

library(MASS)
data("birthwt")

pi <- 1/(1+exp(-(2.15-0.050*birthwt$age-0.015*birthwt$lwt)))

o <- order(pi)
pi <- pi[o]
Y <- birthwt$low[o]

T <- sum(pi*(1-pi))
t <- 1/T*cumsum(pi*(1-pi))
S <- 1/sqrt(T)*cumsum(Y-pi)

plot(t,S, type='l')

Sn <- S[length(S)]
S__ <- max(S-t*Sn) #This is $**, the test statistics of the bridge test

p1 <- 2*pnorm(-abs(Sn))

#Using ks.test implementation to generate pKolmogorov
#In R, 100 observations is enough to invoke the asymptotic #(rather than exact) test.
n_obs <- 100
d <- S__/sqrt(n_obs)
X <- seq(from=0, to=1-d, length.out=n_obs)
p2 <- ks.test(X,punif, exact=F)$p.value

#Fisher's method for generating a unified p-value
X <- -2*(log(p1)+log(p2))
p <- 1-pchisq(X,4)
