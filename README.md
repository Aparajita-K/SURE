# SURE
Selective Update of Relevant Eigenspaces   

### Article:
A. Khan and P. Maji, "Selective Update of Relevant Eigenspaces for Integrative Clustering of Multimodal Data," *IEEE Transactions on Cybernetics*, pp. 1--13, 2020.    
URL: https://ieeexplore.ieee.org/document/9097445

The SURE algorithm permforms integrative clustering on high-dimensional multimodal data sets. A multimodal data set consists of ``M`` modalities X<sub>1</sub>, ..., X<sub>m</sub>, ..., X<sub>M</sub>. Each modality X<sub>m</sub> represents the observations for same set of ``n`` samples from the ``m``-th data source.

![Multimodal Data](/Images/Multimodal-Data.jpg)
The algorithm extracts a ``r`` dimensional joint eigenspace of the integrated data matrix **X** given by   
**X** = [X<sub>1</sub> X<sub>2</sub> ... X<sub>M</sub>],   
from rank ``r`` eigenspaces of the individual modalities X<sub>m</sub>'s.   
The algorithm first the extracts ``r`` eigenspaces of the individual modalities using their singular value decomposition.

![Individual Eigenspaces](/Images/Individual-Eigenspaces.jpg)  
The algorithm constructs the joint eigenspace of **X** by merging the individual left and right singular subspaces, taking into account the intersection and residual of the overlapping subspaces.

![Joint Eigenspace](/Images/Joint-Eigenspace.jpg)   
After extracting the joint eigenspace, *K*-means clustering can be performed on the first ``k`` columns of the ``(n x r)`` joint left-subspace U(**X**) to obtain the clusters of the multimodal data set. 

#### Execution Instructions 
R demo code for GBM and CESC data sets are given in `GBMExample.R` and `CESCExample.R` files, respectively. To run the SURE algorithm on these data sets within the R environment execute:
>source("GBMExample.R")   

>source("CESCExample.R")   

Components of the joint eigenspace are written to the following files:
- ``JointU.txt`` - contains ``(n x r)`` joint left subspace, where ``n`` is the number of samples in the data set and ``r`` is  the rank of the joint subspace.
- ``JointS.txt`` - contains ``(r x r)`` diagonal matrix of singular values.
- ``JointV.txt`` - contains ``(d x r)`` joint right subspace, where ``d`` is the total number of features in the integrated data.    

The principal components of the integrated data matrix can be obtained by multiplying the ``(n x r)`` and the ``(r x r)`` matrices in ``JointU`` and ``JointS``. The joint Principal Components are written to file ``JointPCs.txt``.   
The *K*-means clustering can then be perfomed on the principal components.   
File ``SURE.R`` contains the R implementation of the proposed method as a function ``SURE``. Details of the fuctions is as follows:

Function Name: ``SURE``   
#### #Usage
``SURE(Data, rank, K, modname)``    
Arguments-

``Data``:  A list object containing ``M`` data matrices representing ``M`` different omic data types measured in a set of ``n`` samples. For each matrix, the rows represent samples, and the columns represent genomic features. The matrices in the list can have variable number of columns (features), but they all must have the same number of ``n`` rows(samples).

``rank``: The rank of the individual and joint eigenspaces.

``K``: The number of clusters in the data set.

``modname``: A string array of names of the modalities. Required for modality selection. Example: ``mod=c("RNA","miRNA","CNV")``

#### Example Call:   
```r
#SURE on glioblastoma multiforme (GBM) data set
Data<-list()
Data[[1]] <- as.matrix(read.table("DataSets/GBM/RNA", sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table("DataSets/GBM/miRNA", sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table("DataSets/GBM/CNV", sep=" ",header=TRUE,row.names=1))
K=4
modname=c("RNA","miRNA","CNV")
source("SURE.R")
JointSVD=SURE(Data,rank=K,K=K,modname=modname)
```

For the CESC Data set, log transform sequence based RNA and miRNA modalities before execution of SURE Algorithm.   
#### Example Call:   
```r
#SURE on cervical carcinoma (CESC) data set
DataSet="CESC"
n=124
K=3
rank=K
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/mDNA"), sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/RNA"), sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/miRNA"), sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/Protein"), sep=" ",header=TRUE,row.names=1))
modname=c("mDNA","RNA","miRNA","RPPA")
#Log Transform of Sequence based Gene and miRNA modality
LogData=Data
LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)
source("SURE.R")
JointSVD=SURE(LogData,rank=rank,K=K,modname=modname)
```

Description of the multimodal data sets used in this work, along with their statistical power and survival analyses results is given in ``Appendix.pdf`` file. 

##### Contact Information

Aparajita Khan   
Senior Research Fellow   
Machine Intelligence Unit    
Indian Statistical Institute    
203, B. T. Road    
Kolkata- 700108    
E-mail: aparajitak_r@isical.ac.in,  aparajitakhan1107@gmail.com
