SURE<-function(Dlist,rank=2,K=2,modname=c("mDNA","RNA","miRNA","Protein"))
{
    cat("\nExecuting SURE Algorithm\n")
	Data=Dlist
	tauMax=0.95
	M=length(Dlist)
	Datac=list()
	Datas=list()
	Thetas=list()
	mr=vector(mode="numeric",length=M)
	order=vector(mode="numeric",length=M)
	n=dim(Data[[1]])[1]
	cluster=matrix(0,M,n)
	nsc=rank
	for(i in 1:M)
	{
		Datas[[i]]=scale(Data[[i]],center=T,scale=F)
		sv=svd(Datas[[i]])
		Thetas[[i]]=list(U=as.matrix(sv$u[,1:nsc]),D=sv$d[1:nsc],V=as.matrix(sv$v[,1:nsc]))
		Ids=as.matrix(Thetas[[i]]$U)
		kmbest=kmeans(Ids,K,iter.max=100,nstart=30)
		pve=(kmbest$betweenss/kmbest$totss)
		sm=sum(Thetas[[i]]$D)
		mr[i]=pve
		cat("\nModule",modname[i],"Relevance= ",mr[i])
		cluster[i,]=kmbest$cluster		
	}
	NMImat=matrix(0,M,M)
    for(i in 1:M)
    {
        for(j in 1:M)
            NMImat[i,j]=nmi(cluster[i,],cluster[j,])
        NMImat[i,i]=0
    }
    NMImat=NMImat/max(NMImat)
	order=sort(mr,index.return=TRUE,decreasing=TRUE)$ix
	Eindex=-Inf
	for(th in seq(0,tauMax,0.05))
	{
        seq=order[1]
	    mod=numeric(M)
	    mod[order[1]]=1
	    ThetaAs=list()
	    ThetaAs=list(U=as.matrix(Thetas[[order[1]]]$U),D=Thetas[[order[1]]]$D,V=as.matrix(Thetas[[order[1]]]$V),seq=seq,th=th)
	    taken=order[1]
	    for(j in 1:M)
	    {
	        intIndex=rep(0,M)
	        for(i in 1:M)
		    {
			    if(mod[i]==0)
			    {
				    nmiseq=NMImat[i,seq]
				    intIndex[i]=mean(nmiseq)
			    }
		    }
	        ind=which.max(intIndex)
            nmival=intIndex[ind]
		    if(nmival>th)
		    {
			    G=as.matrix(t(ThetaAs$U)%*%Thetas[[ind]]$U)
			    P=as.matrix(ThetaAs$U%*%G)
			    H=Thetas[[ind]]$U-P
			    Hqr=qr(H)
			    v=qr.Q(Hqr)
				none=length(which(apply(v, 2, var) == 0))
			    if(length(Thetas[[ind]]$D)==1)
				    DVi=as.matrix(Thetas[[ind]]$D*t(Thetas[[ind]]$V))
			    else
				    DVi=diag(Thetas[[ind]]$D)%*%t(Thetas[[ind]]$V)
			    if(length(ThetaAs$D)==1)
				    ns11=as.matrix(ThetaAs$D*t(ThetaAs$V))
			    else
				    ns11=diag(ThetaAs$D)%*%t(ThetaAs$V)
			    ns12=G%*%DVi
			    if(none==dim(v)[2])
			    {
				    Umrg=ThetaAs$U
				    new=cbind(ns11,ns12)
			    }
			    else
			    {
				    t=dim(v)[2]-none
				    vd=as.matrix(v[,1:t]) 
				    Umrg=cbind(ThetaAs$U,vd)
				    ns21=matrix(0,t,dim(ThetaAs$V)[1])
				    ns22=(t(vd)%*%Thetas[[ind]]$U)%*%DVi
				    new=rbind(cbind(ns11,ns12),cbind(ns21,ns22))
			    }		
			    new.sv=svd(new)
			    Urot=Umrg%*%new.sv$u				
			    seq=c(seq,ind)		
			    mn=min(nsc,length(sv$d))
			    ThetaAs=list(U=as.matrix(Urot[,1:mn]),D=new.sv$d[1:mn],V=as.matrix(new.sv$v[,1:mn]),seq=seq,th=th)	
                mod[ind]=1
		        taken=ind			    
		    }
		    else
		        break
	    }	
	    km=kmeans(ThetaAs$U%*%diag(ThetaAs$D),K,nstart=30)
	    pve=km$betweenss/km$totss
	    if(pve>=Eindex)
		{
			Eindex=pve
			FinalSVD=ThetaAs		
		}
	}
	if(length(FinalSVD$D)==1)
		pcU=as.matrix(FinalSVD$U%*%(FinalSVD$D))
	else
		pcU=FinalSVD$U%*%diag(FinalSVD$D)
	write.table(FinalSVD$U,row.names=F,col.names=F,quote=F,file="JointU.txt")
	write.table(diag(FinalSVD$D),row.names=F,col.names=F,quote=F,file="JointS.txt")
	write.table(FinalSVD$V,row.names=F,col.names=F,quote=F,file="JointV.txt")
	write.table(pcU,row.names=F,col.names=F,quote=F,file="JointPCs.txt")
	cat("\nFinal Joint Eigenspace components written to file: JointU.txt, JointS.txt, JointV.txt")
    cat("\nFinal principal components written to file: JointPCs.txt")
    cat("\nSelected Modalities = ",modname[FinalSVD$seq])
    cat("\nValue of threshold tau= ",FinalSVD$th,"\n\n"); 
	return(FinalSVD)
}

nmi<-function(clust1,clust2)
{
	n=length(clust1)
	K=max(clust1)
	ContgCT=matrix(0,K,K)
	for(i in 1:K)
	{
		for(j in 1:K)
		{
			ci=grep(i,clust1)
			tj=grep(j,clust2)
			Iij=intersect(ci,tj)
			ContgCT[i,j]=length(Iij)
		}
	}
	InfCT=0
	HC=0
	HT=0
	for(i in 1:K)
	{
		for(j in 1:K)
		{
		  wl=(n*ContgCT[i,j])/(sum(ContgCT[i,])*sum(ContgCT[,j]))
		  if(wl!=0)
				InfCT=InfCT+(ContgCT[i,j]/n)*log(wl)
		}
	}
	for(i in 1:K)
	{
		ci=sum(ContgCT[i,])/n
		if(ci!=0)
			HC=HC+ci*log(ci)
	}
	HC=-HC
	for(j in 1:K)
	{
		tj=sum(ContgCT[,j])/n
		if(tj!=0)
			HT=HT+tj*log(tj)
	}
	HT=-HT
	NMIsum=2*InfCT/(HC+HT)
	return(NMIsum)
}
