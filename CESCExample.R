DataSet="CESC"
n=124
K=3
rank=K
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/mDNA"), sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/RNA"), sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/miRNA"), sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/Protein"), sep=" ",header=TRUE,row.names=1))
modname=c("mDNA","RNA","miRNA","Protein")
#Log Transform of Sequence based Gene and miRNA modality
LogData=Data
LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)

#****************************************************************************** End of Data Import ***********************************************************************


cat("\nDataset=",DataSet)
cat("\n#Samples=",n)
cat("\nModalities=",modname)
source("SURE.R")
JointSVD=SURE(LogData,rank=rank,K=K,modname=modname)

Dmat=JointSVD$U
clust=kmeans(Dmat,K,iter.max=100,nstart=30)$cluster
cat("\n Clustering: \n",clust)

