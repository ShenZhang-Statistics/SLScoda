library(Rcpp)

##  sparse Laplacian shrinkage for compositional data
SLScoda <- function(X,Z,y,lamb.ls=NULL,nlamb=20,
                   beta0=NULL,penalty="MCP",
                   gam.para=3,A.type="N.1",A=NULL){
  n <- dim(Z)[1]; p <- dim(Z)[2]
  
  ## lambta
  if(is.null(lamb.ls)){
    lambda.factor = 1e-3
    lam0 <- t(Z) %*% (y )/(n)
    lam0 <- max(abs(lam0))  
    lamb1.vec <- exp(seq(from = log(lam0 ), 
                         to = log(lambda.factor * lam0), length = nlamb))
    
    lamb2.vec<- c(exp(seq(from = log(lam0/2 ),
                to = log(lambda.factor * lam0/2), length = round(nlamb/2)) ),0)
  }else{
    lamb1.vec <- lamb.ls[[1]];
    lamb2.vec <- lamb.ls[[2]];
  }
  nlamb <- length(lamb1.vec)
  ## adjacency matrix: A
  if(is.null(A)){
    if(A.type=="N.1"){
      A <- N.1(Z,n,p,c0=3.29)
    }else if(A.type=="N.2"){
      A <- N.2(X,n,p,c0=3.29)
    }else{
      return(message("This A.type is not exist, please check it! "))
    }
  }else{
    A <- A
  }
  beta0 <- rep(0,p)
  output <- SLScoda_lamb(y,Z, A,beta0,lamb1.vec,lamb2.vec,penalty,
                        gam.para,1000,1000, 1e-08, 1e-08, 1.01,1,1e-10 )
 
  opt.res <- list()
  opt.res$lambda1.vec <- lamb1.vec;  
  opt.res$lambda2.vec <- lamb2.vec;
  opt.res$lambda.opt <- c(output$lam1_opt,output$lam2_opt);
  opt.res$beta <- output$beta_opt;
  return(opt.res)
  
}

##  A: adjacency matrix by Pearson’s correlation
N.1 <- function(Z,n,p,c0=3.29)
{
  A <-  abs(cor(Z))
  cutoff <- (exp(2*c0/(n-3)^0.5)-1)/(exp(2*c0/(n-3)^0.5)+1)
  A[abs(A)<cutoff] <- 0
  return(A)
}
##  A: adjacency matrix by  SparCC
N.2 <- function(X,n,p,c0=3.29)
{
  result <- compute.corr.mod(fracs=X, iter=10, th=0.1)
  A <- abs(result$Cor.mat)
  cutoff <- (exp(2*c0/(n-3)^0.5)-1)/(exp(2*c0/(n-3)^0.5)+1)
  A[A<cutoff] <- 0
  return(A)
}
