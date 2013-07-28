lassoMM=function(data,Y,z,grp,D,mu,step,fix,rand,penalty.factor,alpha,showit)
{
	q=ncol(z)
if(missing(D)&(q==1)){D=1}
if(missing(D)&(q>1)){D=0}

if(D==0){
out=lassoMM_nodiag(data,Y,z,grp,mu,step,fix,rand,penalty.factor,alpha,showit)
	}else{
out=lassoMM_diag(data,Y,z,grp,mu,step,fix,rand,penalty.factor,alpha,showit)
	}
out$call=match.call()

out
}
