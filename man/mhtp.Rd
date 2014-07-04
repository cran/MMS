\name{mhtp}
\alias{mhtp}
\docType{package}
\title{
Multiple testing procedure for variable selection in linear mixed models}
\description{
Performs a multiple hypotheses testing in linear mixed models}
\usage{mhtp(data,Y,z,grp,D,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{z}{Random effects matrix. Of size n*q.}
\item{grp}{Grouping variable of length n.}
\item{D}{Logical value. If TRUE, the random effects are considered to be independent, i.e. \code{Psi} is a diagonal matrix. D=TRUE should be used with nested grouping factors.}
\item{fix}{Number of variables which are not submitted to selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{rand}{A vector of length q: each entry k is the position of the random effects number k in the data matrix, 0 otherwise. If \code{z} contains variables that have both a fixed and a random effect, it is advised to not submit them to selection.}
\item{alpha}{A user supplied type I error sequence. Default is (0.1,0.05).}
\item{step}{The algorithm performs at most \code{step} iterations. Default is 3000.}
\item{num}{Number of variables one wishes to order. Default is min(n-1,p-1,30).}
\item{ordre}{Several possible algorithms to order the variables, ordre=c("bolasso","pval","pval_hd","FR"). "bolasso" uses the dyadic algorithm with the Bolasso techniaue, "pval" uses the p-values obtained with a regression on the full set of variables (only when p<n), "pval_hd" uses marginal regression, "FR" uses Forward Regression. Default is "bolasso".}
\item{m}{Number of bootstrapped iteration of the Lasso. Only use if the algorithm is set to "bolasso". Default is m=10.}
\item{show}{Vector of logical values, \code{show=(showordre,showresult,showit)}. Default is (0,0,0). If \code{showordre}=TRUE, show the ordered variables at each step of the algorithm.
if \code{showresult}=TRUE, show the value of the statistics and the estimated quantile at each step of the procedure; if \code{ordre}=bolasso.
if \code{showit}=TRUE, show the iterations of the algorithm.} 
\item{IT}{Number of simulations in the calculation of the quantile. Default is 10000.}
\item{maxq}{Number of maximum multiple hypotheses testing to do. Default is min(log(min(n,p)-1,2),5).}
\item{speed}{Logical value. If TRUE, the algorithm is speeded up once the criteria convergence in \code{beta} and \code{u} are fulfilled. }
}


\details{\code{mhtp} performs fixed effects selection in linear mixed models.
It is a combination of the \link{mht} function from the \link{mht-package} with the algorithm used in the \link{lassop} function;
so you might want to have a look at the help of \link{mht}.

}

\value{A 'mhtp' object is returned for which \code{refit} is available.
\item{data}{List of the user-data: the scaled matrix used in the algorithm, the first column being (1,...,1); Y; z and grp.}
\item{beta}{The estimated vector of fixed effects coefficients. Each row concern a specific user level alpha.}
\item{fitted.values}{Fitted values calculated with the fixed effects and the random effects.}
\item{u}{Matrix with #alpha columns. Each column is the concatenation of the estimated random effects (u_1',...,u_q')' for the user level alpha.}
\item{Psi}{Variance of the random effects. Matrix of dimension q*q.}
\item{sigma_e}{Variance of the residuals.}
\item{it}{Number of iterations of the algorithm.}
\item{quantile}{Array of all the estimated quantiles calculated during the procedure.}
\item{ordrebeta}{All different order that has been used during the procedure.}
\item{converge}{Did the algorithm converge?}
\item{call}{The call that produced this object.}
\item{arg}{List of all the arguments of the function (used to refit the function).}
}

\seealso{\code{\link{mht}}, \code{\link{refit.mhtp}}}

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

########
fit=mhtp(x,y,z,grp,D=0,fix=1,rand=c(1,2),alpha=0.1,num=15)
#fit=mhtp(x,y,z,grp,D=0,fix=1,rand=c(1,2),alpha=0.1,num=15,show=c(1,1,1))
}
}