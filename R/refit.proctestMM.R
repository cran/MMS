refit.proctestMM=function(object,Ynew,z,grp,D,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
{
if(missing(D))
{
	if("proctestMM_nodiag"%in%class(object))
	{			out=refit.proctestMM_nodiag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}
	if("proctestMM_diag"%in%class(object))
	{					out=refit.proctestMM_diag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}
}else{
	if(D==0){
out=refit.proctestMM_nodiag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}else{
out=refit.proctestMM_diag(object,Ynew,z,grp,fix,rand,alpha,step,num,ordre,m,show,IT,maxq,speed,...)
	}
	
	
}
	
out$call=match.call()
		
out	
}