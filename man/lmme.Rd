\name{lmme}
\alias{lmme}
\docType{package}
\title{Linear mixed model estimation}
\description{
Estimation of fixed effects and random effects in a linear mixed model}
\usage{
lmme(data,Y,z,grp,step,showit)}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{z}{Random effects matrix. Of size n*q.}
\item{grp}{Grouping variable of length n.}
\item{step}{The algorithm performs at most \code{maxIT} iterations. Default is 3000.}
\item{showit}{Logical value. If TRUE, shows the convergence process of the algorithm. Default is FALSE.}
}

\value{
\item{data}{List of the user-data: the scaled matrix used in the algorithm, the first column being (1,...,1); Y; z and grp.}
\item{beta}{Estimation of the selected fixed effects.}
\item{sigma_u}{Variance of the random effects. It is a vector of length q.}
\item{sigma_e}{Variance of the noise.}
\item{fitted.values}{Fitted values calculated with the fixed effects and the random effects.}
\item{it}{Number of iterations of the algorithm.}
\item{converge}{Logical. TRUE if the algorithm has converged, FALSE otherwise.}
\item{call}{The call that produced this object.}
}

\details{\code{lmme} performs an ML-estimation of fixed effects and random effects in a linear mixed model when no selection is involved.}

\examples{
\dontrun{
N <- 20           # number of groups
p <- 80            # number of covariates (including intercept)
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
for (k in 1:n) y[k] <- x[k,]%*%beta + t(z[k,])%*%bi[,k] + epsilon[k]

########
fit=lmme(x,y,z,grp)
}
}
