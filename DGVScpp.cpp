#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double updatephi1(mat beta,
                  mat pstar,
                  double lamb1,
                  double a,
                  double b)
{
  int T = beta.n_cols;
  
  rowvec phi1seq = {0.5, 0.8, 0.85, 0.90, 0.95, 0.99};
  int l = phi1seq.n_rows;
  
  double sumpstar0;
  double sumbeta0;
  double sumbetabeta;
  double sumbeta;
  rowvec meanlogpost(l);
  double phi1;
  
  sumpstar0 = accu(pstar.col(0));
  sumbeta0 = accu(pstar.col(0) % beta.col(0) % beta.col(0));
  sumbetabeta = accu(pstar.cols(1,T-1) % beta.cols(0,T-2) % beta.cols(1,T-1));
  sumbeta = accu(pstar.cols(1,T-1) % beta.cols(0,T-2) % beta.cols(0,T-2));
  meanlogpost = 0.5*log(1-phi1seq%phi1seq)*sumpstar0 - 0.5*(1-phi1seq%phi1seq)/lamb1*sumbeta0 +
              phi1seq/lamb1*sumbetabeta - phi1seq%phi1seq/(2*lamb1)*sumbeta +
             (a-1)*log(1 +phi1seq) + (b-1)*log(1-phi1seq);
  phi1 = phi1seq(index_max(meanlogpost));
  
  return phi1;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double updatelamb(mat beta,
                  mat pstar,
                  double phi1,
                  double c,
                  double a,
                  double b)
{
  int T = beta.n_cols;
  
  double A = 0;
  double B = 0;
  
  double lamb1;
  

  A = accu((1-pstar) % beta % beta) +
      (1-phi1*phi1)*c*accu(pstar.col(0) % beta.col(0) % beta.col(0)) +
      c*accu(pstar.cols(1,T-1) % (beta.cols(1,T-1) - phi1*beta.cols(0,T-2)) % (beta.cols(1,T-1) - phi1*beta.cols(0,T-2))) +
      2*c*b;
  B = c*(accu(1-pstar) + accu(pstar) + 2*(a+1));
  
  lamb1 = A/B;
  
  return lamb1;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List EMVScpp(rowvec Y,
             mat Ft,
             double Theta,
             mat beta0,
             double lamb10, double c,
             bool estlamb) 
{
  int ite = 10000;
  int p = Ft.n_rows;
  int T = Ft.n_cols;
  
  //mixing weights
  cube theta = zeros<cube>(p, T, ite);
  cube pstar = zeros<cube>(p, T, ite);
  double psist;
  double psi;
  double mu;
  double m0 = 0;
  
  //precisions
  rowvec dn = zeros<rowvec>(T);
  dn(0) = 1;
  rowvec d = zeros<rowvec>(T);
  d(0) = 1;
  double delta = 0.95;
  rowvec dr = zeros<rowvec>(T);
  dr(0) = 1;
  mat nustar(T, ite);

  
  //autoregressive parameter
  double phi0 = 0;
  rowvec phi1 = zeros<rowvec>(ite+1);
  phi1(0) = 0.8;
  double aphi = 20;
  double bphi = 1.5;
  
  //slab and spike variances
  rowvec lamb1 = zeros<rowvec>(ite);
  lamb1(0) = lamb10;
  rowvec lamb0 = zeros<rowvec>(ite);
  lamb0(0) = c*lamb1(0);
  double alamb = 3;
  double blamb = 500;
  
  //regression coefficients
  cube beta = zeros<cube>(p,T,ite);
  beta.slice(0) = beta0;
  colvec eD(p);
  mat D = zeros<mat>(p, p);
  colvec mu2(p); 
  mat Dinv = zeros<mat>(p, p);
  mat Siginv = zeros<mat>(p, p);
  colvec eSig0(p);
  
  //stop criterion
  double e1 = 0.001;
  
  //double mse;
  List output;
  
  bool open = true;
  int i = 0;
  
  while(open){
    
    //-----E-step
    for(int t=1; t < T; t++){
      for(int j=0; j < p; j++){
        //mixing weights 
        psist = normpdf(beta(j,(t-1),i), phi0, sqrt(lamb1(i)/(1-phi1(i)*phi1(i))));
        theta(j,t,i) = Theta*psist/(Theta * psist + (1-Theta)*normpdf(beta(j,(t-1),i), m0, sqrt(lamb0(i))));
        
        
        mu = phi0 + phi1(i)*(beta(j,(t-1),i) - phi0); 
        psi = normpdf(beta(j,t,i), mu, sqrt(lamb1(i)));
        pstar(j,t,i) = theta(j,t,i)*psi/(theta(j,t,i)*psi +
          (1-theta(j,t,i))*normpdf(beta(j,t,i), m0, sqrt(lamb0(i))));
       
      }
      
      //precisions
      dn(t) = delta * dn(t-1) + 1;
      dr(t) = Y(t) - as_scalar(Ft.col(t).t()*beta.slice(i).col(t));
      d(t) = delta * d(t-1) + dr(t)*dr(t);
    }
    //mixing weights 
    pstar.slice(i).col(0) = theta.slice(i).col(1);
    pstar.slice(i).row(0) =  ones<rowvec>(T);
    theta.slice(i).row(0) =  ones<rowvec>(T); 
    
    //precisions
    nustar(T-1,i) = dn(T-1)/d(T-1);
    for(int t = (T-2); t > 0 ; t--){
      nustar(t,i) = (1-delta)*dn(t)/d(t) + delta*nustar((t+1),i);
    }
    
    //-----M-Step
    
    //autoregressive parameter
    phi1(i+1) = updatephi1(beta.slice(i), pstar.slice(i), lamb1(i), aphi, bphi);
    
    //slab variance
    if(estlamb){
      lamb1(i+1) = updatelamb(beta.slice(i),pstar.slice(i), phi1(i+1), c, alamb, blamb);
      lamb0(i+1) = c*lamb1(i+1);
    }else{
      lamb1(i+1) = lamb10;
      lamb0(i+1) = c*lamb1(i+1);
    }
    
    //regression coefficients
    for(int t=1; t < (T-1); t++){
      eD = pstar.slice(i).col(t)/lamb1(i+1) + (1-pstar.slice(i).col(t))/(lamb0(i+1)) +
            (phi1(i+1)*phi1(i+1)*pstar.slice(i).col(t+1)/lamb1(i+1));
      D = diagmat(eD);
      mu2 = nustar(t,i)*Y(t)*Ft.col(t) +
            phi1(i+1)/lamb1(i+1)*beta.slice(i).col(t-1)%pstar.slice(i).col(t) +
            (phi1(i+1)/lamb1(i+1)*beta.slice(i).col(t+1)%pstar.slice(i).col(t+1));
      
      Dinv = diagmat(1/eD);
      Siginv = Dinv - nustar(t,i)*Dinv*
        (Ft.col(t)*Ft.col(t).t()/(1+nustar(t,i)*as_scalar(Ft.col(t).t()*Dinv*Ft.col(t))))*Dinv;
      beta.slice(i+1).col(t) = Siginv*mu2; 
    }
    //t=T
    eD = pstar.slice(i).col(T-1)/lamb1(i+1) + (1-pstar.slice(i).col(T-1))/(lamb0(i+1));
    D = diagmat(eD);
    mu2 = nustar(T-1,i)*Y(T-1)*Ft.col(T-1) +
          phi1(i+1)/lamb1(i+1)*beta.slice(i).col(T-2)%pstar.slice(i).col(T-1);
    
    Dinv = diagmat(1/eD);
    Siginv = Dinv - nustar(T-1,i)*Dinv*
      (Ft.col(T-1)*Ft.col(T-1).t()/(1+nustar(T-1,i)*as_scalar(Ft.col(T-1).t()*Dinv*Ft.col(T-1))))*Dinv;
    beta.slice(i+1).col(T-1) = Siginv*mu2;
      
    //t=0
    eSig0 = (1-phi1(i+1)*phi1(i+1))*pstar.slice(i).col(0)/lamb1(i+1) +
      (1-pstar.slice(i).col(0))/(lamb0(i+1)) + phi1(i+1)*phi1(i+1)*pstar.slice(i).col(1)/lamb1(i+1);
    beta.slice(i+1).col(0) = phi1(i+1)/lamb1(i+1)*diagmat(1/eSig0)*
                             beta.slice(i+1).col(1)%pstar.slice(i).col(1);
      
    open = ((abs(beta.slice(i+1)-beta.slice(i))/(beta.slice(i) + 0.01)).max() > e1 || 
           (abs(lamb1(i+1)-lamb1(i))/(lamb1(i) + 0.01)) > e1) && i < (ite-2);
    i++;
  }
  
  
  
  output["beta"] = beta.slice(i);
  output["pstar"] = pstar.slice(i-1);
  output["nustar"] = nustar.col(i-1);
  output["lamb1"] = lamb1(i);
  output["c"] = c;
  output["phi1"] = phi1(i);
  output["ite"] = i-1;
  //output["mse"] = mse;
  
  return output;
}

