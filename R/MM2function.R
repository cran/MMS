MM2=function(data2,yhat,m,sigma,maxqdep,maxordre,var_nonselect,random,showtest,showit,showresult,correspondance,ORDREBETA2,indice,aV2,choix_ordre,alpha,alph,IT)
{#####################
# MM2 commence ###
#####################
compteurordre=nrow(ORDREBETA2)
ntot=nrow(data2)
p=ncol(data2)
#on ordonne les variables pour avoir les var_nonselect_fixed en premier et ensuite les var_nonselect_random!=0, pour pouvoir faire les test correctement sans selection var_nonselect=sum des deux

ordre=dyadiqueordre(data2,yhat,m,maxordre,var_nonselect=var_nonselect,showtest=showtest,showordre=FALSE,random=random)# donne l'ordre (dans ordre) et le nombre de fois ou l'algo a redemarré (dans prob)	

	b=ordre$ordre

bb=correspondance[2,b]
#mean that data2[ordre]=data[correspondance[2,ordre]]=data[,bb]

if(showit){print(bb)}
	ORDREBETA=bb[1:maxordre]

#on complete l'ordre par les variables restantes
	a=match(1:p,b)
	b=c(b,(1:p)[which(is.na(a))])
	
XI_ord=data2[,b] #on a ainsi les XI ordonnées dans XI_ord

dec=decompbaseortho(XI_ord)
#on rajoute dans nonind2 les dernieres variables, celle qui n'ont pas d'utilitÈs puisque dans Rn
nonind=dec$nonind
U=dec$U
if(p>(ntot+length(nonind)))
{nonind2=c(nonind,(ntot+1+length(nonind)):p)	
	Uchap=U[,-nonind2]
}else{
	if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}
		}
dim_X=ncol(Uchap)			#nombre de variables utiles


beta2=lm(yhat~Uchap-1)$coefficients #decomposition de Y2 dans la base orthonormal (X(1),..,X(p))
beta2[-which(beta2!=0)]=0

indice2=var_nonselect-1-sum(nonind<=var_nonselect)

ktest=var_nonselect-sum(nonind<=var_nonselect)
nbr_test=var_nonselect

calcul=numeric(0)
maxq=min(log(min(ntot,dim_X)-1,2),maxqdep)
aV=array(0,c(length(alpha),maxq,dim_X)) #on y met tous les quantiles
	T=1
	while((T>0)&&(dim_X>ktest+1))
		{
		maxq=min(log(min(ntot,dim_X)-ktest-1,2),maxqdep)
	#on est a ktest fixé, on regarde dans tout ce qu'on a fait avant ou indice>ktest si on a deja calculé le quantile

	#on a indice de longueur compteur-1, aV_compteur avec les quantiles, var_select de dim compteur-1, et ORDREBETA2 
if(ktest>indice2)
{abc=compteurordre#dim de indice/ORDREBETA2/var_select
	i=0
	I=numeric(0)
	TT=0
	#A=numeric(0)
	while((TT==0)&&(i<abc))#dim(ORDREBETA2)=abc*maxordre
	{i=i+1
	a=numeric(0)
	K=numeric(0)
	for(j in 1:nbr_test)
	{a=c(a,sum(ORDREBETA[j]==ORDREBETA2[i,1:nbr_test]))}
	
	if((sum(a)==nbr_test)&&(indice[2,i]>=ktest)&&(indice[1,i])<=ktest){TT=1
	#	K=c(K,ktest)
		I=c(I,i)}
	}

	if(length(I)!=0)
	{	aV[,,ktest]=aV2[,,ktest,I]#get(paste("aV",I,sep="_"))[,,ktest]
		calcul=c(calcul,0)#on met 0 si le quantile est deja calculé
		}else{
		if(showresult){print(paste("ktest=",ktest))}
		quant=quantileprocbol(XI_ord,ktest,alpha,IT=IT,maxq=maxq,sigma=sigma)
		aV[,1:maxq,ktest]=quant$quantile	
		calcul=c(calcul,1)#on met 1 si on a calculé un quantile manquant
		if(showresult){print(aV[,,ktest])}
		}
indice2=indice2+1
}
		#test S	
		bb=numeric(0)
	
		for(m in 0:(maxq-1))
			{
			a=sum(beta2[(ktest+1):(ktest+2^m)]^2)/sigma
			b=a>aV[alph,m+1,ktest]#F #1 si on doit rejeter le test, 0 sinon
			bb=c(bb,b) #on met tous les tests de Hk
			}
		if(length(which(bb!=0))>0) #securitÈ a la base, inutile maintenant? a verifier
			{bb[-which(bb!=0)]=0}else{bb=matrix(0,1,m)}
		
		if(sum(bb)>0){T=1
		ktest=ktest+1
			
		nbr_test=nbr_test+1
		if(length(nonind)>0)
		{for(i in 1:length(nonind))
			if(sum(nonind==(nbr_test+1))==1){nbr_test=nbr_test+1}}
											
		}else{
			T=0
			#print(ktest)
		} #on rejete Hk s'il y a au moins un test qui rejete
	}#fin while
if(ktest==dim_X){k0=dim_X}else{k0=ktest}

NBR=nbr_test #rÈsultat contenant le nombre de variables sÈlectionnÈes
NBR_effect=k0


if(sum(calcul)!=0)
{#on rajoute que si on a calculer au moins un quantile
	aV3=array(0,c(dim(aV2)[1],dim(aV2)[2],dim(aV2)[3],nrow(ORDREBETA2)+1))
	aV3[,,,1:nrow(ORDREBETA2)]=aV2
	compteurordre=compteurordre+1
	ORDREBETA2=rbind(ORDREBETA2,ORDREBETA)
	aV3[,,,nrow(ORDREBETA2)]=aV
	aV2=aV3

	#aV3[,,,nrow(ORDREBETA2)]=aV#assign(paste("aV",compteurordre,sep="_"),aV)
	indice=cbind(indice,c(var_nonselect-sum(nonind<=var_nonselect),NBR_effect))#c(indice,NBR)
	
	#je veux qu'on SI on a calculer un quantile d'un ktest en plus, on remplace l'existant par le nouveau avec indice superieur
	}	
	ind=ORDREBETA[1:NBR]

out=list(ind=ind,indice=indice,ORDREBETA2=ORDREBETA2,aV2=aV2,compteurordre=compteurordre,ORDREBETA=ORDREBETA,aV=aV)
###########################
#  fin MM2 ###
##########################
}