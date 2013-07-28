\name{lassoMM}
\alias{lassoMM}
\docType{package}
\title{
 L1-penalization in linear mixed models}
\description{
Performs a  L1-penalization in linear mixed models}
\usage{
lassoMM(data,Y,z,grp,D,mu,step,fix,rand,penalty.factor,alpha,showit)
}


\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{z}{Random effects matrix. Of size n*q.}
\item{grp}{Grouping variable of length n.}
\item{D}{Logical value. If TRUE, the random effects are considered to be independent, i.e. \code{Psi} is a diagonal matrix. D=TRUE should be used with nested grouping factors.}
\item{mu}{Positive regularization number to be used for the Lasso.}
\item{step}{The algorithm performs at most \code{step} iterations. Default is 3000.}
\item{fix}{Number of variables which are not submitted to selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{rand}{A vector of length q: each entry k is the position of the random effects number k in the data matrix, 0 otherwise. If z contains variables that have both a fixed and a random effect, it is advised to not submit them to selection.}
\item{penalty.factor}{Argument of 'glmnet'. Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables that are not in 1:fix.}
\item{alpha}{Argument of 'glmnet'. The elasticnet mixing parameter, with \eqn{ 0\le \alpha \le 1}. \code{alpha}=1 is the lasso penalty, and \code{alpha}=0 the ridge penalty.}
\item{showit}{Logical value. If TRUE, shows the iterations of the algorithm. Default is FALSE.}
}


\details{
This function performs fixed effects selection in linear mixed models through a L1-penalization of the log-likelihood of the marginal model. The method optimizes a criterion via a multicycle ECM algorithm at the regularization parameter mu. \cr
Two algorithms are available: one when the random effects are assumed to be independent (D=TRUE) and one when they are not (D=FALSE).
Selection on the random is only performed when D=TRUE.}

\value{A 'lassoMM' object is returned.
\item{data}{List of the user-data: the scaled matrix used in the algorithm, the first column being (1,...,1); Y and Z, which is the design matrix of the random effects.}
\item{beta}{Estimation of the fixed effects.}
\item{fitted.values}{Fitted values calculated with the fixed effects and the random effects.}
\item{Psi}{Variance of the random effects. Matrix of dimension q*q.}
\item{sigma_e}{Variance of the residuals.}
\item{it}{Number of iterations of the algorithm.}
\item{converge}{Logical. TRUE if the algorithm has converged, FALSE otherwise.}
\item{u}{Vector of the concatenation of the estimated random effects (u_1',...,u_q')'.}
\item{call}{The call that produced this object.}
\item{mu}{The penalty used in the algorithm.}
}

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
for (k in 1:n) y[k] <- x[k,]\%*\%beta + t(z[k,])\%*\%bi[,k] + epsilon[k]

########
fit=lassoMM(x,y,z,grp,D=1,mu=0.2,fix=1,rand=c(1,2))
}
}
