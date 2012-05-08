\name{refit.proctestMM}
\alias{refit.proctestMM}
\docType{package}
\title{Linear mixed model estimation}
\description{
Estimation of fixed effects and random effects in a linear mixed model}
\usage{
\method{refit}{proctestMM}(object,Ynew,z,grp,var_nonselect_fixed,var_nonselect_random,alpha,step,maxordre,choix_ordre=c("bolasso","pval","pval_hd"),m,showordre,IT,maxq,showtest,showresult,showit,speed,...)
}

\arguments{
\item{object}{Object of class "proctestMM".}
\item{Ynew}{Response variable of length n.}
\item{z}{Random effects matrix. Of size n*q.}
\item{grp}{Grouping variable of length n.}
\item{var_nonselect_fixed}{Number of variables which are not submitted to selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{var_nonselect_random}{A vector of length q: each entry k is the position of the random effects number k in the data matrix, 0 otherwise. If z contains variables that have both a fixed and a random effect, it is advised to not submit them to selection.}
\item{alpha}{A user supplied type I error sequence. Default is (0.1,0.05).}
\item{step}{The algorithm performs at most \code{maxIT} iterations. Default is 3000.}
\item{maxordre}{Number of variables one wishes to order. Default is min(n-1,p-1,30).}
\item{choix_ordre}{Three possible algorithms to order the variables. "bolasso" uses the regularization path of the Bolasso technique, "pval" uses the p-values from the test of nullity of the coefficients obtained from a linear regression on the full set of variables (only when p<n), the last "pval_hd" makes regression on one variable at a time. Default is "bolasso".}
\item{m}{Number of bootstrapped iteration of the Lasso. Only use if the algorithm is set to "bolasso". Default is m=10.}
\item{showordre}{Logical value. If TRUE, show the ordered variables at each step of the algorithm. Default is FALSE.}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 10000.}
\item{maxq}{Number of maximum multiple hypotheses testing to do. Default is min(log(min(n,p)-1,2),5).}
\item{showtest}{Logical value. Only use if the algorithm is set to "bolasso". if TRUE, show the number of regularization parameters tested to show the advancement of the dyadic algorithm. Default is FALSE.}
\item{showresult}{Logical value. if TRUE, show the value of the statistics and the estimated quantile at each step of the procedure. Default is FALSE.}
\item{showit}{Logical value. If TRUE, shows the iterations of the algorithm. Default is FALSE.}
\item{speed}{Logical value. If TRUE, the algorithm is speeded up once the criterion convergence in \code{beta} and \code{u} is fulfilled. }
\item{\dots}{not used}
}

\value{A 'proctestMM object' is returned.
\item{data}{List of the user-data: the scaled matrix used in the algorithm, the first column being (1,...,1); Y; z and grp.}
\item{beta}{Estimation of the selected fixed effects.}
\item{fitted.values}{Fitted values calculated with the fixed effects and the random effects.}
\item{u}{Matrix with #alpha columns. Each column is the concatenation of the estimated random effects (u_1',...,u_q')' for the user level alpha.}
\item{sigma_u}{Variance of the random effects. It is a vector of length q.}
\item{sigma_e}{Variance of the noise.}
\item{it}{Number of iterations of the algorithm.}
\item{quantile}{Array of all the estimated quantiles calculated during the procedure.}
\item{ordrebeta}{All different order that has been used during the procedure.}
\item{call}{The call that produced this object.}
\item{arg}{List of all the arguments of the function.}
}

\seealso{\code{\link{proctestMM}}}

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
for (k in 1:n) y[k] <- x[k,]%*%beta + t(z[k,])%*%bi[,k] + epsilon[k]

########
fit=proctestMM(x,y,z,grp,var_nonselect_fixed=1,var_nonselect_random=c(1,2),showordre=TRUE,showit=TRUE,showresult=TRUE,alpha=0.1,maxordre=15)
fit2=refit(fit)
}
}