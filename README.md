# SLScoda Description Document

This is a description document to introduce the usage for SLScoda based on a paper titled "Sparse Laplacian Shrinkage for High-dimensional Regression with Correlated Compositional Data" (Zhang et al. 2024+). 
We first generalize the sparse Laplacian shrinkage
(SLS; Huang et al. 2011) method to compositional data and propose a new method, SLScoda, to address variable selection in high-dimensional linear log-contrast model with correlated compositional data. 
SLScoda can simultaneously consider the compositional effect,
the high dimensionality, and the component correlation of compositional data. SLScoda outperforms its competitor in terms of prediction error, estimation accuracy, and variable selection. 
The related R code and the example are summarized in a single folder SLScodaRcode. Please set this folder as the working directory.

 
## 1 SLScoda function 

There are the main function of SLScoda method in **SLScoda.R**.

```{r , eval=FALSE}
SLScoda(X,Z,y,lamb.ls=NULL,nlamb=20, beta0=NULL, penalty="MCP", gam.para=3, A.type="N.1", A=NULL)
```

**Input parameter**

* X                
  - compositional data; the element of X is not 0.
  - n x p data matrix; n: sample size, p: compositional variables size. 
* Z
  - log(X)
* y    
  - Response data
* lamb.ls
  - The list of Lambda: one is the vector of lambda_1, and the other is the vector of lambda_2.
* nlamb  
  - If lamb.ls is null, nlamb is the number of lamdata generated.
* beta0
  - Initial estimation for regression coefficients.
* penalty
  - Penalty function: "MCP" or "SCAD".
* gam.para
  - Regularization parameter $\gamma$: if the penalty function is "MCP", $\gamma=3$; if the penalty function is "SCAD", $\gamma=3.27$.
* A.type
  - The type of constructing a adjacency matrix: 
  - "N.1": Pearson’s correlation.
  - "N.2": Sparcc (Friedman and Alm, 2012).
* A
  - Adjacency matrix. 

**Output result**

* beta       
  - Estimation for regression coefficients
* lambda1.vec, lambda2.vec      
  - the vectors of lambda_1 and lambda_2
* lambda.opt       
  - Optimal lambda1 and lambda2.

# 2 Example 

There is a simple example for using the SLScoda function in **SLScoda.R**. First load the required R package and R files.

```{r , eval=FALSE}
library(Rcpp)
sourceCpp("SLScodaRcode/SLScoda.cfunc.cpp")
source("SLScodaRcode/SLScoda.R")
source("SLScodaRcode/sparcc.modify.R")
```

The user of SLScoda should prepare two data sets: compositional data and response data. For example,  generate the simulation data.
```{r , eval=FALSE}
library(MASS)
set.seed(1212)
n=500; p=50
sig.eps <- 0.5; rho <- 0.5
beta <- c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2,rep(0,p-8)) # the true coefficients 

Sigma <- matrix(0,p,p)
for(i in 2:p){
  for(j in 1:(i-1)){
    Sigma[i,j] <- rho^{i-j}
  }
}
Sigma <- Sigma + t(Sigma) + diag(p)
theta <- c(rep(log(0.5*p),each=5),rep(0,each=p-5))

W <- mvrnorm(n,mu=theta, Sigma=Sigma)  
X <-  exp(W) / rowSums(exp(W))  # compositional data
eps <- rnorm(n, mean=0, sd=sig.eps)
Z <- log(X)
y <- Z %*% beta + eps  #  response data
```

Call ``SLScoda`` function in **SLScoda.R**.
```{r , eval=FALSE}
res.est <- SLScoda(X,Z,y,lamb.ls=NULL, nlamb=20, penalty="MCP", gam.para=3, A.type="N.1", A=NULL);
res.est$beta[1:10]
```

# 3 Reference

* Zhang, S. et al. (2024+), "Sparse Laplacian Shrinkage for High-dimensional
Regression with Correlated Compositional Data",  *Submit to journal*.
* Huang, J., Ma, S., Li, H. and Zhang, C.-H. (2011). The sparse Laplacian shrinkage estimator for high-dimensianl regression, The Annals of Statistics 39(4): 2021–2046.
* Friedman, J., and Alm, E. J. (2012), "Inferring correlation networks from genomic survey data", *PLOS Computational Biology*, 8(9), e1002687.
