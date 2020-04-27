DataSet="GBM"
n=168
K=4
rank=K
Data<-list()
Data[[1]] <- as.matrix(read.table("DataSets/GBM/RNA", sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table("DataSets/GBM/miRNA", sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table("DataSets/GBM/CNV", sep=" ",header=TRUE,row.names=1))
modname=c("RNA","miRNA","CNV")

#****************************************************************************** End of Data Import ***********************************************************************


cat("\nDataset=",DataSet)
cat("\n#Samples=",n)
cat("\nModalities=",modname)
source("SURE.R")
JointSVD=SURE(Data,rank=rank,K=K,modname=modname)

Dmat=JointSVD$U
clust=kmeans(Dmat,K,iter.max=100,nstart=30)$cluster
cat("\n Clustering: \n",clust)
