#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::interfaces(r, cpp)]]
// Augmented Lagrangian method (ALM) nested with coordinate descent (CD) for zero-sum constraints
// with MCP and SLS regression problem. 
// [[Rcpp::export]]
Rcpp::List MCPSLS_lamb2(arma::vec y,
                       arma::mat Z,
					   arma::mat A,
                       arma::vec beta,
                       arma::vec lamb1,	
                       int inner_maxiter,
                       int outer_maxiter,
					   double lamb2,
					   double gam_para,
                       double inner_eps,
                       double outer_eps,
                       double mu_ratio,
                       double u_ini,
                       double tol
) {
  int inner_inter, outer_inter, i, j; /* index for loop */
  double inner_err, outer_err; /* convergence for loop */
  int n = y.size();
  int p = Z.n_cols;
  int nlam = lamb1.size();
  int npass = 0;  
  
  double u;
  double S;
  double tt;
  double C_linear;
  double alpha_mu = 0;
  double eta;
  
  arma::vec gamma(p);
  arma::vec xi(p);
  arma::vec r(n);
  arma::vec d(p);
  arma::vec mse(nlam);
  arma::vec df(nlam);

  arma::vec beta_old_inner(p), beta_old_outer(p), Ldiff(p);
  arma::mat path(p, nlam);
  arma::uvec ids;

  for(i = 0; i < p; i++) {
    gamma(i) = sum(Z.col(i) % Z.col(i));
	  d(i)= sum(abs(A.row(i)))  - abs(A(i,i));	
  }
  gamma = gamma / n;
  
  if(mu_ratio < 1) {
    u_ini = 0;
    outer_maxiter = 1;
  }

  //Ni.zeros();
  //df.zeros();

  /*----------- outer loop for path starts-------------*/

  for(j = 0; j < nlam; j++) {
    /*----------- middle loop ALM starts-------------*/
    //laml = pf * lamb1(j);
    u = u_ini;
    //r = y - join_rows(Z, Zc) * beta;
    r = y - Z * beta;
    C_linear = sum(beta);
    xi = gamma + u + lamb2 * d;	
    outer_inter = 1;
    outer_err = 1;

    while((outer_inter <= outer_maxiter) && (outer_err > outer_eps)) {

      /*----------- inner loop Block-wise-lasso starts-------------*/
      beta_old_outer = beta;
      inner_inter = 1;
      inner_err = 1;

      while( (inner_inter <= inner_maxiter) && (inner_err > inner_eps)) {
        npass++;
        beta_old_inner = beta;
	
        //update coeffients for composition variables
        for(i = 0; i < p; i++) {
          r = r + Z.col(i) * beta(i);
          C_linear = C_linear - beta(i);
		      eta = lamb2 * ( sum(A.col(i) % beta)- A(i,i) * beta(i));	
		  S = sum(Z.col(i) % r) / n + eta - u * C_linear - alpha_mu; 
		  
          if( std::abs(S) >  lamb1(j) * gam_para * xi(i) ) {
            beta(i) =  S / xi(i);
            r = r - Z.col(i) * beta(i);
            C_linear = C_linear + beta(i);
          }else{
			      tt = std::abs(S) - lamb1(j);
           
      			if(tt >0.0){
      				if(S > 0.0) {
      				  beta(i) =  gam_para * (tt)/(gam_para * xi(i) -1);
      				}else if(S < 0.0){
      				  beta(i) =  -gam_para * (tt)/(gam_para * xi(i) -1);
      				} else {
      				  beta(i) = 0.0;
      				}
      				r = r - Z.col(i) * beta(i);
      				C_linear = C_linear + beta(i); 
      			}else{				
              beta(i) = 0.0;
            }  			
    		  }	  
        }
        
        inner_inter++;
        //check convergence for inner loop
        Ldiff = beta - beta_old_inner;
        inner_err = sum((beta_old_inner - beta) % (beta_old_inner - beta));				
        
	      if( inner_err < (inner_eps * p) ) {
          inner_err = inner_err / std::max(sum(beta_old_inner % beta_old_inner), 1e-30);
        } 
        
      }
      /*----------- inner loop Block-wise-lasso ends-------------*/
      outer_inter++;

      //update lagrange multiplier and penalty term
      alpha_mu = alpha_mu + u * C_linear;
      u = u * mu_ratio;
      xi = gamma + u + lamb2 * d; 
      
      //check convergence for middle loop
      outer_err = std::abs(C_linear);
      if(outer_err <= outer_eps) {
        Ldiff = beta - beta_old_outer;
        outer_err = sum(Ldiff % Ldiff) / std::max(sum(beta_old_outer % beta_old_outer), 1e-30);
      }      
    }
   
    /*----------- middle loop ALM ends-------------*/    
    ids = find(abs(beta.subvec(0, p - 1)) < tol);
    beta.elem(ids).fill(0.0);
    path.col(j) = beta;
	
	mse(j) = mean((y- Z * beta)%(y- Z * beta));
	df(j) = sum(abs(beta) > 0);
  }
  
  /*----------- outer loop for path ends-------------*/
  //formulate output
  arma::vec output_lambda;
  arma::vec output_df;
  arma::mat output_path(p, j);
  int output_ll;
  output_lambda = lamb1.subvec(0, j-1);
  output_path = path.cols(0, j-1);
  output_ll = npass;
  return Rcpp::List::create(Rcpp::Named("beta") = output_path,
                            Rcpp::Named("lamb1") = output_lambda,
							Rcpp::Named("lamb2") = lamb2,
                            Rcpp::Named("npass") = output_ll,
							Rcpp::Named("mse") = mse,
							Rcpp::Named("df") = df,
                            Rcpp::Named("penalty") = "MCP");
}


// Augmented Lagrangian method (ALM) nested with coordinate descent (CD) for zero-sum constraints
// with SCAD and SLS regression problem. 
// [[Rcpp::export]]
Rcpp::List SCADSLS_lamb2(arma::vec y,
                        arma::mat Z,
                        arma::mat A,
                        arma::vec beta,
                        arma::vec lamb1,	
                        int inner_maxiter,
                        int outer_maxiter,
                        double lamb2,
                        double gam_para,
                        double inner_eps,
                        double outer_eps,
                        double mu_ratio,
                        double u_ini,
                        double tol
) {
  int inner_inter, outer_inter, i, j; /* index for loop */
  double inner_err, outer_err; /* convergence for loop */
  int n = y.size();
  int p = Z.n_cols;
  int nlam = lamb1.size();
  int npass = 0;  
  
  double u;
  double S;
  double tt;
  double C_linear;
  double alpha_mu = 0;
  double eta;
  
  arma::vec gamma(p);
  arma::vec xi(p);
  arma::vec r(n);
  arma::vec d(p);
  arma::vec mse(nlam);
  arma::vec df(nlam);
    
  arma::vec beta_old_inner(p), beta_old_outer(p), Ldiff(p);  
  arma::mat path(p, nlam);
  arma::uvec ids;
  
  for(i = 0; i < p; i++) {
    gamma(i) = sum(Z.col(i) % Z.col(i));
    d(i)= sum(abs(A.row(i)))  - abs(A(i,i));	
  }
  gamma = gamma / n;
  
  if(mu_ratio < 1) {
    u_ini = 0;
    outer_maxiter = 1;
  }
  
  /*----------- outer loop for path starts-------------*/
  
  for(j = 0; j < nlam; j++) {
    /*----------- middle loop ALM starts-------------*/
    u = u_ini;
    r = y - Z * beta;
    C_linear = sum(beta);
    xi = gamma + u + lamb2 * d;	
    outer_inter = 1;
    outer_err = 1;
    
    while((outer_inter <= outer_maxiter) && (outer_err > outer_eps)) {
      
      /*----------- inner loop Block-wise-lasso starts-------------*/
      beta_old_outer = beta;
      inner_inter = 1;
      inner_err = 1;
      
      while( (inner_inter <= inner_maxiter) && (inner_err > inner_eps)) {
        npass++;
        beta_old_inner = beta; 
        
        //update coeffients for composition variables
        for(i = 0; i < p; i++) {
          r = r + Z.col(i) * beta(i);
          C_linear = C_linear - beta(i);
          eta = lamb2 * ( sum(A.col(i) % beta)- A(i,i) * beta(i));	
          S = sum(Z.col(i) % r) / n + eta - u * C_linear - alpha_mu; 
          
          if( std::abs(S) >  lamb1(j) * gam_para * xi(i) ) {
            beta(i) =  S / xi(i);
          }else if((std::abs(S) <=  lamb1(j) * gam_para * xi(i)) & 
                   (std::abs(S) > 2* lamb1(j)  * xi(i))){
            if(S > 0.0) {
              beta(i) =  ((gam_para-1)* S - gam_para * lamb1(j))/((gam_para-1)*xi(i) -1);
            }else{
              beta(i) =  ((gam_para-1)* S + gam_para * lamb1(j))/((gam_para-1)*xi(i) -1);
            }
          }else {
            tt = std::abs(S) - lamb1(j);
            if(tt > 0.0){
              if(S > 0.0) {
                beta(i) =  tt/xi(i);
              }else {
                beta(i) =  -tt/ xi(i);
              }
            }else{
              beta(i) =0.0;
            }
          }
          r = r - Z.col(i) * beta(i);
          C_linear = C_linear + beta(i);	  
        }
        
        inner_inter++;
        //check convergence for inner loop
        Ldiff = beta - beta_old_inner;
        inner_err = sum((beta_old_inner - beta) % (beta_old_inner - beta));				
        
        if( inner_err < (inner_eps * p) ) {
          inner_err = inner_err / std::max(sum(beta_old_inner % beta_old_inner), 1e-30);
        } 
        
      }
      /*----------- inner loop Block-wise-lasso ends-------------*/
      outer_inter++;
      
      //update lagrange multiplier and penalty term
      alpha_mu = alpha_mu + u * C_linear;
      u = u * mu_ratio;
      xi = gamma + u + lamb2 * d; 
      
      //check convergence for middle loop
      outer_err = std::abs(C_linear);
      if(outer_err <= outer_eps) {
        Ldiff = beta - beta_old_outer;
        outer_err = sum(Ldiff % Ldiff) / std::max(sum(beta_old_outer % beta_old_outer), 1e-30);
      }      
    }
   
    /*----------- middle loop ALM ends-------------*/    
    ids = find(abs(beta.subvec(0, p - 1)) < tol);
    beta.elem(ids).fill(0.0);
    path.col(j) = beta;	
	
	mse(j) = mean((y- Z * beta)%(y- Z * beta));
	df(j) = sum(abs(beta) > 0);
  }
  
  /*----------- outer loop for path ends-------------*/
  //formulate output
  arma::vec output_lambda;
  arma::vec output_df;
  arma::mat output_path(p, j);
  int output_ll;
  output_lambda = lamb1.subvec(0, j-1);
  output_path = path.cols(0, j-1);
  output_ll = npass;
  return Rcpp::List::create(Rcpp::Named("beta") = output_path,
                            Rcpp::Named("lamb1") = output_lambda,
							Rcpp::Named("lamb2") = lamb2,
                            Rcpp::Named("npass") = output_ll,
							Rcpp::Named("mse") = mse,
							Rcpp::Named("df") = df,
                            Rcpp::Named("penalty") = "SCAD");
}


// [[Rcpp::export]]
Rcpp::List SLScoda_lamb(arma::vec y,
                     arma::mat Z,
                     arma::mat A,
                     arma::vec beta,
                     arma::vec lamb1,
                     arma::vec lamb2,
                     std::string penalty="MCP",
                     double gam_para=3,
                     int inner_maxiter=1000,
                     int outer_maxiter=1000,
                     double inner_eps=1e-08,
                     double outer_eps=1e-08,
                     double mu_ratio=1.01,
                     double u_ini=1,
                     double tol=1e-10  ) {
  Rcpp::List output;
  int nlam2 = lamb2.size();
  double gic_init = 1/0.0;
  arma::vec beta_opt;
  double lam1_opt = lamb1(0);
  double lam2_opt = lamb2(0);
  double gic_min;
  arma::vec gic_new;
  int idmin;
  
  int n = y.size();
  int p = Z.n_cols; 
  for(int i=0; i < nlam2; i++){
	if(penalty == "MCP"){
	  output = MCPSLS_lamb2(y,Z, A, beta,lamb1, inner_maxiter,outer_maxiter,
							  lamb2(i),gam_para, inner_eps, outer_eps,
							  mu_ratio, u_ini, tol);
	}else if(penalty == "SCAD"){
	  output = SCADSLS_lamb2(y,Z, A, beta,lamb1, inner_maxiter,outer_maxiter,
							   lamb2(i),gam_para, inner_eps, outer_eps,
							   mu_ratio, u_ini, tol);
	}else{
		Rcpp::stop("The penalty parameter is incorrect. Please select 'MCP' or 'SCAD'!");
	}    
    
    arma::mat beta_mat = output["beta"];
    arma::vec mse =output["mse"] ;
    arma::vec df = output["df"];
    int max_pn =  (p+n + abs(p-n))/2;
    df = (df>=2)  % (df-1) + (df<2) * 0;
    gic_new = log(mse) + (log(max_pn) * log(log(n)) / n) * df;
    gic_min = min(gic_new);
    
    if(gic_min < gic_init ){
      idmin = arma::index_min(gic_new);
      beta_opt = beta_mat.col(idmin);
      lam1_opt = lamb1(idmin);
      lam2_opt = lamb2(i);
    }
  }
  return Rcpp::List::create(Rcpp::Named("lam1_opt") = lam1_opt,
                            Rcpp::Named("lam2_opt") = lam2_opt,
                            Rcpp::Named("beta_opt") = beta_opt,
                            Rcpp::Named("penalty") = penalty);

}
