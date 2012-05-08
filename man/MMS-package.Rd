\name{MMS-package}
\alias{MMS-package}
\docType{package}
\title{
Fixed effects Selection in Linear Mixed Models}
\description{
This package performs fixed effects selection in a linear mixed model
}

\details{
\tabular{ll}{
Package: \tab MMS \cr
Type: \tab Package \cr
Version: \tab 1.0 \cr
Date: \tab 2012-05-08 \cr
License: \tab GPL-3 \cr
LazyLoad: \tab yes \cr
}
Two major functions: \code{\link{lassoMM}} and \code{\link{proctestMM}}. The first performs a L1-penalization of the log-likelihood of the marginal model. The  last is a combination of the \code{\link{procbol}} method of the \code{\link{mht-package}} and a multicycle ECM algorithm.}
\author{
Florian Rohart \cr
Maintainer: florian.rohart@gmail.com
}


\references{
Multiple hypotheses testing for variable selection; F. Rohart \cr
Model-consistent sparse estimation through the bootstrap; F. Bach \cr
Fixed effects selection in high dimensional linear mixed models; F. Rohart & al }


\examples{
\dontrun{
N <- 20           # number of groups
p <- 20            # number of covariates (including intercept)
q <- 2            # number of random effect covariates
ni <- rep(6,N)    # observations per group
n <- sum(ni)   # total number of observations

grp <- factor(rep(1:N,ni)) # grouping variable
grp=rbind(grp,grp)

beta <- c(1,2,4,3,rep(0,p-3)) # fixed-effects coefficients
x <- cbind(1,matrix(rnorm(n*p),nrow=n)) # design matrix

u1=rnorm(N,0,sd=sqrt(2))
u2=rnorm(N,0,sd=sqrt(2))
bi1 <- rep(u1,ni) 
bi2 <- rep(u2,ni)
bi <- rbind(bi1,bi2)

z=x[,1:2,drop=FALSE]
   
epsilon=rnorm(120)
y <- numeric(n)
for (k in 1:n) y[k] <- x[k,]\%*\%beta + t(z[k,])\%*\%bi[,k] + epsilon[k]
######
fit1=lassoMM(x,y,z,grp,mu=0.2,var_nonselect_fixed=1,var_nonselect_random=c(1,2),showit=TRUE)
#fit=lassoMM(x,y,z,grp,mu=0.2,var_nonselect_fixed=1,var_nonselect_random=c(1,2),showit=FALSE)
fit2=proctestMM(x,y,z,grp,var_nonselect_fixed=1,var_nonselect_random=c(1,2),alpha=0.1,maxordre=15,showordre=TRUE,showit=TRUE,showresult=TRUE)
#fit2=proctestMM(x,y,z,grp,var_nonselect_fixed=1,var_nonselect_random=c(1,2),alpha=0.1,maxordre=15,showordre=FALSE,showit=FALSE,showresult=FALSE)
}
}