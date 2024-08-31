library(Rcpp)
##  A: adjacency matrix by Pearson’s correlation
A1 <- function(Z,n,p,c0=3.29){
  A <-  abs(cor(Z))
  cutoff<- (exp(2*c0/n^0.5)-1)/(exp(2*c0/n^0.5)+1)
  A[abs(A)<cutoff] <- 0
  return(A)
}

##  A: adjacency matrix by  SparCC
A2 <- function(X,n,p,c0=3.29){
  result <- compute.corr.mod(fracs=X, iter=10, th=0.1)
  A <- abs(result$Cor.mat)
  cutoff <- (exp(2*c0/n^0.5)-1)/(exp(2*c0/n^0.5)+1)
  A[A<cutoff] <- 0
  return(A)
}

##  sparse Laplacian shrinkage for compositional data
SLScoda <- function(X,Z,y,lambda_ls=NULL,nlambda=20,
                   beta0=NULL, penalty="MCP",
                   gamma_para=3, A_type="A1",A=NULL){
  n <- dim(Z)[1]; p <- dim(Z)[2]
  
  ## lambta
  if(is.null(lambda_ls)){
    lambda_factor = 1e-3
    lambda0 <- t(Z) %*% (y )/(n)
    lambda0 <- max(abs(lambda0))  
    lambda1_vec <- exp(seq(from = log(lambda0 ), 
                         to = log(lambda_factor * lambda0), length = nlambda))
    
    lambda2_vec<- exp(seq(from = log(lambda0/2 ),
                to = log(lambda_factor * lambda0/2), length = round(nlambda/2)) )
  }else{
    lambda1_vec <- lambda_ls[[1]];
    lambda2_vec <- lambda_ls[[2]];
  }
  nlambda <- length(lambda1_vec)
  ## adjacency matrix: A
  if(is.null(A)){
    if(A_type=="A1"){
      A <- A1(Z,n,p,c0=3.29)
    }else if(A_type=="A2"){
      A <- A2(X,n,p,c0=3.29)
    }else{
      return(message("This A_type is not exist, please check it! "))
    }
  }
  
  beta0 <- rep(0,p)
  output <- SLScoda_lambda2(y,Z, A, beta0,lambda1_vec,lambda2_vec,penalty,
                        gamma_para,1000,1000,
                        1e-08, 1e-08, 1.01,1,1e-10 )
 
  opt_res <- list()
  opt_res$lambda1_vec <- lambda1_vec;  
  opt_res$lambda2_vec <- lambda2_vec;
  opt_res$lambda_opt <- c(output$lambda1_opt,output$lambda2_opt);
  opt_res$beta <- output$beta_opt;
  return(opt_res)
}

