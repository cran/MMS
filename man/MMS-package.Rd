\name{MMS-package}
\alias{MMS-package}
\alias{MMS}
\docType{package}
\title{
\packageTitle{MMS}
}
\description{
\packageDescription{MMS}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MMS}
\packageIndices{MMS}
Two major functions:\cr
\code{\link{lassop}} performs fixed effects selection in linear mixed models through a L1-penalization of the log-likelihood of the marginal model. The method optimizes a criterion via a multicycle ECM algorithm. \cr
\code{\link{mhtp}} performs fixed effects selection in linear mixed models through a modification of the previous algorithm which allows the use of any variable selection method developed in linear models.
It is a combination of the \code{\link{mht}} method from the \code{\link{mht-package}} and the modified algorithm.
}
\author{
\packageAuthor{MMS}

Maintainer: \packageMaintainer{MMS}
}
\references{
Multiple hypotheses testing for variable selection; F. Rohart, 2016, https://doi.org/10.1111/anzs.12157 \cr
Selection of fixed effects in high dimensional linear mixed models using a multicycle ECM algorithm; F. Rohart & al, 2014, https://doi.org/10.1016/j.csda.2014.06.022\cr
}

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
fit0=lmme(x[,1:4],y,z,grp)
fit1=lassop(x,y,z,grp,mu=0.2,fix=1,rand=c(1,2))
#fit1=lassop(x,y,z,grp,mu=0.2,fix=1,rand=c(1,2),showit=TRUE)

fit2=mhtp(x,y,z,grp,fix=1,rand=c(1,2),alpha=0.1,num=15)
#fit2=mhtp(x,y,z,grp,fix=1,rand=c(1,2),alpha=0.1,num=15,show=c(1,1,1))
}
}
