# SURE
Selective Update of Relevant Eigenspaces

### Article:
A. Khan and P. Maji, "Selective Update of Relevant Eigenspaces for Integrative Clustering of Multimodal Data," *IEEE Transactions on Cybernetics*, pp. 1--13, 2020. 

#### Execution Instructions for R
R demo code for GBM and CESC data sets are given in `GBMExample.R` and `CESCExample.R` files, respectively. To run the SURE algorithm on these data sets within the R environment execute:
>source("GBMExample.R")

>source("CESCExample.R")

Components of the Joint eigenspace are written to the following files:
- ``JointU.txt`` - contains ``(n x r)`` joint left subspace, where ``n`` is the number of samples in the data set and ``r`` is  the rank of the joint subspace.
- ``JointS.txt`` - contains ``(r x r)`` diagonal matrix of singular values.
- ``JointV.txt`` - contains ``(d x r)`` joint right subspace, where ``d`` is the total number of features in the integrated data.
 The principal components of the integrated data matrix can be obtained by multiplying the ``(n x r)`` and the ``(r x r)`` matrices in ``JointU`` and ``JointS``. The joint Principal Components are written to file ``JointPCs.txt``.
The *K*-means clustering can then be perfomed on the principal components.
