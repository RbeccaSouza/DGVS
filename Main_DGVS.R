#-------------------------------------------#
#   Dynamic Graphical Variable Selection    #
#-------------------------------------------#
#--- DGVS

#-- Setions:

# Simulation
# Estimation
# Auxiliary Functions
# Results
#   DAG 1 - DAG Simple
#   DAG 2 - DAG Complicated
#   DAG 3 - DAG Complex

#---- Packages -----
library(mvtnorm)
library(Matrix)
library(lqmm)
library(coda)
library(igraph)
library(R.matlab)
library(invgamma)

library(multicon)
library(testit)
library(DGM)
#library(multdyn)
library(cowplot)
library(png)
library(ggplot2)
library(pacman)
library(reshape2)
library(grid)
library(ggplotify)
library(dplyr) # easier data wrangling 
library(viridis) # colour blind friendly palette, works in B&W also
library(ggExtra)
library(animation)
library(hrbrthemes)

library(parallel)
library(foreach)
library(doParallel)

library(eegkit)

library("Rcpp")


#-------------------------------------------------------#
#                  Dynamic EMVS                         #
#-------------------------------------------------------#
#----EMVS cpp----
#Rcpp
sourceCpp("DGVScpp.cpp") 


#----------------------------------------------------#
#                  simulation                        #
#----------------------------------------------------#
#----- Simulation Graph -----

simu <- function(p, T., adj=array(dim=c(p,p,T.)), phi0, phi1, lamb0, lamb1, v){
  Y <- array(0,dim=c(p,T.))
  Ft <- array(dim=c(p,T.,p))
  if(length(v)>1){
    V <- v
  }else{
    V <- rep(v, T.)
  }
  
  beta <- array(dim=c(p,T.,p))
  gama <- array(dim=c(p,p,T.))
  W <- array(dim=c(p,p,T.))
  G <- array(dim=c(p,p,T.))
  H <- matrix(NA,p,T.)
  
  C <- array(dim=c(p,p,T.))
  m <- matrix(NA,p,T.)
  
  for(t in 1:T.){
    for(r in 1:p){
      gama[r,,t] <- c(1,t(adj[-r,r,t]))
    }
  }
  
  #--Sampling beta
  for(r in 1:p){
    #t=0
    m[,1] <- phi0*gama[r,,1]
    C[,,1] <- diag(gama[r,,1]*lamb1/(1-phi1^2) + (1-gama[r,,1])*lamb0)
    beta[,1,r] <- rmvnorm(1, m[,1], C[,,1])
    W[,,1] <- diag(gama[r,,1]*lamb1 + (1 - gama[r,,1])*lamb0)
    G[,,1] <- diag(gama[r,,1]*phi1)
    H[,1] <- phi0*gama[r,,1]
    
    for(t in 2:T.){
      W[,,t] <- diag(gama[r,,t]*lamb1 + (1 - gama[r,,t])*lamb0)
      G[,,t] <- diag(gama[r,,t]*phi1)
      H[,t] <- phi0*gama[r,,t]
      
      beta[,t,r] <- H[,t] + G[,,t] %*% (beta[,(t-1),r]-H[,t]) +
                    t(rmvnorm(1, rep(0,p), W[,,t]))
      Ft[,t,r] <- c(1,Y[-r,t])
      Y[r,t] <- t(Ft[,t,r]) %*% beta[,t,r] + rnorm(1, 0, sqrt(V[T.]))
    }
  }
  
  return(list(Y=t(Y), beta=beta))
}

simu.sample <- function(N, p, T., adj=array(dim=c(p,p,T.)), phi0, phi1, lamb0, lamb1, v){
  Y <- array(dim=c(T., p, N))
  beta <- array(dim=c(p, T., p, N))
  
  for(s in 1:N){
    suj <- simu(p, T., adj, phi0, phi1, lamb0, lamb1, v)
    Y[,,s] <- suj$Y
    beta[,,,s] <- suj$beta
  }
  
  return(list(Y=Y, beta=beta))
}

#------------------------------------------------------#
#                    Estimation                         #
#------------------------------------------------------#
#----Estimation----

grid_c<- function(Y, F., p, T., Theta, beta0, lamb1.0, c, estlamb=T){
  nc <- length(c)
  
  resul <- list()
  filtro <- list()
  
  msfe <- NULL
  resul[[1]] <- EMVScpp(Y, F., Theta, beta0, lamb1.0, c[1], estlamb)
  
  filtro[[1]] <- dss.lpl(Y, F.,(resul[[1]]$pstar>=0.5)*1, resul[[1]]$phi1,
                         resul[[1]]$lamb1, c[1]*resul[[1]]$lamb1, resul[[1]]$nustar)
  error <- filtro[[1]]$et
  msfe[1] <- sum(error^2)
  
  if(nc > 1){
    for(l in 2:nc){
      resul[[l]] <- EMVScpp(Y, F., Theta, beta0, lamb1.0, c[l], estlamb)
      filtro[[l]] <- dss.lpl(Y, F.,(resul[[l]]$pstar>=0.5)*1, resul[[l]]$phi1,
                             resul[[l]]$lamb1, c[l]*resul[[l]]$lamb1, resul[[l]]$nustar)
      error <- filtro[[l]]$et
      msfe[l] <- sum(error^2)
    }
  }
  menor <- which(msfe==min(msfe))
  output <- resul[[menor]]
  
  return(output)
}


EMVS.subj <- function(Y, p, T., Theta, beta0=array(NA,dim=c(p,T.,p)), lamb1.0, c, estlamb=T){
  
  resultados.beta <- array(dim=c(p,T.,p))
  resultados.pstar <- array(dim=c(p,T.,p))
  
  ite <- NULL
  lamb1 <- NULL
  phi1<- NULL
  nustar <- array(dim=c(T., p))
  c_ <- NULL
  
  for(r in 1:p){ 
    if(is.na(beta0[1,1,r])){ 
      beta0[1:p,,r] <- t(rmvnorm(1, rep(2,p), diag(rep(lamb1.0/(1-0.9^2),p))))
    }
  }
  
  #-----------EMVS em paralelo-------
  #oEMVS <- EMVS(Y.r, F.r, p, T., Theta, beta0=beta0[,,r], lamb1.0, lamb0)
  #Detecta o n?mero de n?cleos dispon?veis e cria cluster

  #cl <- parallel :: makeCluster (detectCores ())
  #doParallel :: registerDoParallel (cl)
  #time_foreach <- system.time ({
  #  oEMVS <-foreach(r=1:p,.packages = "Rcpp", .export=c("grid_c","dss.lpl"),
  #                  .noexport=c("EMVScpp","updatephi1","updatelamb") )%dopar%{
  #    sourceCpp("EMVScpp.cpp")
  #    grid_c(Y[,r], F.=rbind(rep(1,T.+1),t(Y[,-r])), p, T., Theta, beta0=beta0[,,r], lamb1.0, c, estlamb)
  #  }
  #})
  #  parallel :: stopCluster (cl)
  
  time_foreach <- NULL
  for(r in 1:p){
      time_foreach[r] <- system.time(
      oEMVS <- grid_c(Y[,r], F.=rbind(rep(1,T.),t(Y[,-r])), p, T., Theta,
                           beta0=beta0[,,r], lamb1.0, c, estlamb)
      )[3]   
   
      resultados.beta[,,r]  <- oEMVS$beta
      resultados.pstar[,,r] <- oEMVS$pstar
      ite[r] <- oEMVS$ite
      lamb1[r] <- oEMVS$lamb1
      phi1[r] <- oEMVS$phi1
      nustar[,r] <- oEMVS$nustar
      c_[r] <- oEMVS$c
      
      remove(oEMVS)
      gc()
  }
  
  output <- list(beta=resultados.beta, pstar=resultados.pstar,
                 time= time_foreach, nustar=nustar, lamb1=lamb1, c=c_,
                 phi1=phi1, ite=ite)
  return(output)
}

EMVS.sample <- function(Y, p, T., N, Theta=1, beta0=array(NA,dim=c(p,T.,p,N)), lamb1.0=100, c, estlamb=T){
  
  est.beta <- array(dim=c(p, T., p, N))
  est.pstar <- array(dim=c(p, T., p, N))
  lamb1 <- array(dim=c(p, N))
  ite <- array(dim=c(p, N))
  time <- array(dim=c(N,p))
  phi1 <- array(dim=c(p, N))
  nustar <- array(dim=c(T., p, N))
  c_ <- array(dim=c(p,N))
  
  for(s in 1:N){
    Y.s <- Y[,,s]
    resu.s <- EMVS.subj(Y.s, p, T., Theta, beta0[,,,s],lamb1.0,c, estlamb)
    est.beta[,,,s] <- resu.s$beta
    est.pstar[,,,s] <- resu.s$pstar
    lamb1[,s] <- resu.s$lamb1
    ite[,s] <- resu.s$ite
    time[s,] <- resu.s$time
    phi1[,s] <- resu.s$phi1
    nustar[,,s] <- resu.s$nustar
    c_[,s] <- resu.s$c
    
    remove(Y.s)
    remove(resu.s)
  }
  
  output <- list(beta=est.beta, 
                 pstar=est.pstar,
                 nustar=nustar,
                 lamb1=lamb1,
                 c=c_,
                 phi1=phi1,
                 Time=time,
                 ite=ite)
  
  return(output)
}

#------------------------------------------------------#
#              Auxiliary Functions                     #
#------------------------------------------------------#
#---- Auxiliary Functions ----

get.adj <- function(est.pstar, limiar, limiar.t){
  d <- dim(est.pstar)
  T. <- d[2]
  p <- d[1]
  
  if(length(d) == 4){
    N <- d[4]
    Adjw.t <- array(0,dim=c(p, p, T., N))
    Adjw <- array(0,dim=c(p, p, N))
    Adj.t <- array(0,dim=c(p, p, T., N))
    Adj <- array(0,dim=c(p, p, N))
  
    for(s in 1:N){
      #ADJACENCy MATRIX at time
      for(t in 1:(T.)){
        for(r in 1:p){
          parent.t <- (1:p)[-r]
          Adjw.t[parent.t, r, t, s] <- est.pstar[-1,t,r,s]
        }
      }
      Adj.t <- (abs(Adjw.t)>=limiar)*1

      #ADJACENCY MATRIX 
      Adjw[,, s] <- apply(Adjw.t[, , -1, s], c(1,2), mean)
      
      sum.adj.t <- apply(Adj.t[,,-1, s],c(1,2), sum)
      Adj[,,s] <- (sum.adj.t>=((T.-1)*limiar.t))*1
    }
  
  }else{
    Adjw.t <- array(0,dim=c(p, p, T.))
    Adjw <- array(0,dim=c(p, p))
    Adj.t <- array(0,dim=c(p, p, T.))
    Adj <- array(0,dim=c(p, p))
    
    #ADJACENCy MATRIX at time
    for(t in 1:(T.)){
      for(r in 1:p){
        parent.t <- (1:p)[-r]
        Adjw.t[parent.t, r, t] <- est.pstar[-1,t,r]
      }
    }  
    Adj.t <- (abs(Adjw.t)>=limiar)*1
    
    #ADJACENCY MATRIX 
    Adjw <- apply(Adjw.t[,,-1], c(1,2), mean) #without time 0

    sum.adj.t <- apply(Adj.t[,,-1],c(1,2), sum)
    Adj <- (sum.adj.t>=((T.-1)*limiar.t))*1
  }
  return(list(Adjw.t=Adjw.t,
              Adjw=Adjw,
              Adj.t=Adj.t,
              Adj=Adj))
}

binom.nettest <- function(adj, alter="two.sided", fdr=0.05){
  mydim=dim(adj)
  M = sum(adj, na.rm=T) # total edges over all N subjects, all R(R-1) edges
  
  if (length(mydim) == 3) { # without runs
    N=mydim[3] # No. of subjects
    N_Comp=mydim[1]
    adj_ = apply(adj, c(1,2), sum, na.rm=T)
    
  } else if (length(mydim) == 4) { # with runs
    N=mydim[4] # No. of subjects
    N_runs=mydim[3]
    N_Comp=mydim[1]
    N=N*N_runs # adjust N by for no. of runs
    adj_ = apply(adj, c(1,2), sum, na.rm=T) # sum acrosss subjects and runs
  }
  
  # binom test for every edge occurance
  p0 = M/(N*N_Comp*(N_Comp-1)) # H0 edge probability
  
  p = array(NA, dim=c(N_Comp,N_Comp))
  for (i in 1:N_Comp) {
    for (j in 1:N_Comp) {
      if (i!=j) {
        tmp=binom.test(adj_[i,j],N,p=p0, alternative=alter)
        p[i,j]=tmp$p.value
      }
    }
  }
  
  # FDR
  p_fdr=matrix(p.adjust(p, method = "fdr"),N_Comp,N_Comp)
  adj_fdr=adj_
  adj_fdr[p_fdr>=fdr]=0
  
  store=list()
  store$p0=p0
  store$p=p
  store$p_fdr=p_fdr
  store$adj=adj_/N
  store$adj_fdr=adj_fdr/N # significant proportions
  
  return(store)
}


get.change <- function(Adj){
  d <- dim(Adj)
  T. <- d[3]
  p <- d[1] 
  change <- NULL
  fdr_ <- array(dim=c(p,p,T.))
  
  if(length(d)==3){
    i <- 1
    fdr_[,,1] <- Adj[,,1]
    for(t in 2:T.){
      if(!identical(Adj[,,t], Adj[,,(t-1)])){
        change[i] <- t
        fdr_[,,i] <- Adj[,,t]
        i <- i+1
      }
    }  
  }else{
    i <- 1
    t1 <- binom.nettest(Adj[,,1,], alter = "greater", fdr = 0.05)
    adj.t1 <- (t1$adj>0)*1
    fdr_[,,1] <- t1$adj
    for(t in 2:T.){
      t2 <- binom.nettest(Adj[,,t,], alter = "greater", fdr = 0.05)
      adj.t2 <- (t2$adj>0)*1
      if(!identical(adj.t1, adj.t2)){
        change[i] <- t
        fdr_[,,i] <- t2$adj
        i <- i+1
      }
      adj.t1 <- adj.t2
    }
  }
  Nn<- length(change) + 1
  fdr <- array(dim=c(p,p,Nn))
  fdr <- fdr_[,,1:Nn]
  
  return(list(time=c(1,change), fdr=fdr))
}

dss.lpl <- function(Yt, Ft, gama, phi1, lamb1, lamb0, nu){
  m0 = 0
  CS0 = lamb1/(1-phi1^2)
  nt = NULL
  nt[1] <- 3
  delta = 0.90
  
  p = nrow(Ft)     
  CS0 = CS0*diag(p)
  m0 = rep(m0,p)
  Nt = length(Yt) # the length of the time series + t0
  
  Wt = array(0, dim=c(p,p,Nt))
  Wt[,,1] = diag(rep(lamb1,p))
  
  Gt = array(0, dim=c(p,p,Nt))
  Gt[,,1] = diag(rep(phi1,p))
  
  # Set up allocation matrices, including the priors
  mt = array(0,dim=c(p,Nt))
  Ct = array(0,dim=c(p,p,Nt))
  Ct[,,1] = CS0
  
  at = array(0,dim=c(p,Nt))
  at[,1]=m0
  Rt = array(0,dim=c(p,p,Nt))
  
  ft = numeric(Nt)
  qt = numeric(Nt)
  et = numeric(Nt)
  
  
  lpl = numeric(Nt)
  
  for (t in 2:Nt){
    Wt[,,t] = diag(gama[,t]*lamb1/(1-phi1^2) + (1-gama[,t])*lamb0)
    Gt[,,t] = diag(gama[,t]*phi1)
    
    #Seguindo West and Harrison paginas 111 e 112 e Rockova(2020)
    
    # Posterior at {t-1}:
    Rt[,,t] = Gt[,,t]%*%Ct[,,(t-1)]%*%t(Gt[,,t]) + Wt[,,t]
    at[,t] = Gt[,,t]%*%mt[,t-1] 
    
    # One-step forecast: 
    ft[t] = t(Ft[,t]) %*% at[,t]
    qt[t] = t(Ft[,t])%*% Rt[,,t] %*%Ft[,t] + nu[t]
    et[t] = Yt[t] - ft[t]
    
    # Posterior at t: 
    At = (Rt[,,t] %*% Ft[,t])/qt[t]
    mt[,t] = at[,t] + (At*et[t])
    Ct[,,t] = Rt[,,t]-(At %*% t(At))*qt[t]
    
    nt[t] = delta*nt[t-1] + 1
    
    # Log Predictive Likelihood 
    lpl[t] = lgamma((nt[t]+1)/2)-lgamma(nt[t]/2)-0.5*log(pi*nt[t]*qt[t])-((nt[t]+1)/2)*log(1+(1/nt[t])*et[t]^2/qt[t])
  }
  output <- list(et = et, lpl=lpl)
  return(output)
}

# pruning function (one subject)
poda1 <- function(adj, Y, phi1, lamb1, lamb0, nu, e = 20, limiar=0.5){ #for one subject

  p <- dim(adj)[1]
  T.<- dim(adj)[3]
  am <- array(0, dim=c(p, p, T.))
  
  # determine bidirectional edges
  bi <- array(0, dim=c(p, p, T.))
  for(t in 1:T.){
    am[,,t]=adj[,,t]
    
    for (i in 1:p){
      for (j in 1:p) {
        if (adj[i,j,t] == 1 & adj[j,i,t] == 1){
          bi[i,j,t] <- 1
          
          lpls <- array(NA, dim=c(p, p, 2))
          #model for i with father j
          gamai <- rbind(rep(1,T.),adj[-i,i,])
          Fi <- rbind(rep(1,T.),t(Y[,-i]))
          lpli <- dss.lpl(Y[,i], Fi, gamai, phi1[i], lamb1[i], lamb0[i], nu[,i])
          lpli <- lpli$lpl
          
          #model for j with father i
          gamaj <- rbind(rep(1,T.),adj[-j,j, ])
          Fj <- rbind(rep(1,T.),t(Y[,-j]))
          lplj <- dss.lpl(Y[,j], Fj , gamaj, phi1[j], lamb1[j], lamb0[j], nu[,j])
          lplj <- lplj$lpl
          
          lpls[i,j,1] = lpls[j,i,1] = lpli[t] + lplj[t]
          
          #model for i without father j 
          adjaux <- adj[,i,]
          adjaux[j,t] <- 0 # removing the father j
          gamai_j <- rbind(rep(1,T.),adjaux[-i,])
          lpli_j <- dss.lpl(Y[,i], Fi , gamai_j, phi1[i], lamb1[i], lamb0[i], nu[,i])
          lpli_j <- lpli_j$lpl
          
          lpls[i,j,2] = lplj[t] + lpli_j[t]
          
          #model for j without father i
          adjaux <- adj[,j,]
          adjaux[i,t] <- 0 # removing the father i
          gamaj_i <- rbind(rep(1,T.),adjaux[-j,])
          lplj_i <- dss.lpl(Y[,j], Fj , gamaj_i, phi1[j], lamb1[j], lamb0[j], nu[,j])
          lplj_i <- lplj_i$lpl
          
          lpls[j,i,2] = lpli[t] + lplj_i[t]
          
          if (lpls[i,j,1] - e <= max(lpls[i,j,2], lpls[j,i,2]) ) {
            # if unidirectional lpl is larger than bidirectional with a
            # bayes factor penatly, take the simpler unidirectional model.
            if (lpls[i,j,2] > lpls[j,i,2]) {
              am[i,j,t] = 1; am[j,i,t] = 0
            } else if (lpls[i,j,2] < lpls[j,i,2])  {
              am[i,j,t] = 0; am[j,i,t] = 1
            }
          }
          
        }
      }
    }
    
  }
  sum.adj.t <- apply(am[,,-1],c(1,2), sum)
  Adj <- (sum.adj.t>=(T.*limiar))*1

  resu <- list(Adj.t=am, Adj=Adj)
 
  return(resu)
}

# pruning function (multiple subjects)
poda <- function(Adj, Y, phi1, lamb1, lamb0, nu, e=20, limiar=0.5){
  d <- dim(Adj)
  p <- d[1]
  T. <- d[3]
  
  if(length(d)==4){
    N <- d[4]
    
    cl <- parallel :: makeCluster (detectCores ())
    doParallel :: registerDoParallel (cl)
    poda.sample <-foreach(s=1:N, .export =c("poda1","dss.lpl"))%dopar%{
      poda1(Adj[,,,s],Y[,,s], phi1[,s], lamb1[,s], lamb0[,s], nu[,,s], e, limiar)
    }
    parallel :: stopCluster (cl)
  
    g.poda<-array(NA, dim=c(p,p,N))
    g.poda.t<-array(NA, dim=c(p,p,T.,N))
    for(s in 1:N){
      g.poda[,,s] <- poda.sample[[s]]$Adj
      g.poda.t[,,,s] <- poda.sample[[s]]$Adj.t
    }
    
  }else{
    poda.1 <- poda1(Adj,Y, phi1, lamb1, lamb0, nu, e, limiar)
    g.poda <- poda.1$Adj
    g.poda.t <- poda.1$Adj.t
  }
  
  output<-list(Adj=g.poda, Adj.t=g.poda.t)
  return(output)
}

binom.nettest.t <- function(adj, alter="two.sided", fdr=0.05){
  d <- dim(adj)
  p <- d[1]
  T. <- d[3]
  fdr_t <- array(dim=c(p,p,T.))
  adj_t <- array(dim=c(p,p,T.))
  for(t in 1:T.){
    stat <- binom.nettest(adj[,,t,], alter="two.sided", fdr=0.05)
    fdr_t[,,t] <- stat$adj_fdr
    adj_t[,,t] <- stat$adj
  }
  store <- list(adj.t=adj_t, adj.fdr.t=fdr_t)
  return(store)
}

perf <- function(x, true) {
  
  d = dim(x)
  Nn=d[1]
  if (length(d) == 3) {
    N=d[3]
  } else if (length(d) == 2) {
    x=array(x, dim=c(Nn,Nn,1))
    N=1
  }
  
  subj=array(NA,dim=c(N,8))
  cases=array(NA,dim=c(N,4))
  
  for (i in 1:N) {
    TP = sum(x[,,i] & true)
    FP = sum((x[,,i] - true) == 1)
    FN = sum((true - x[,,i]) == 1)
    TN = sum(!x[,,i] & !true) - ncol(x[,,i])
    
    cases[i,]=c(TP,FP,FN,TN)
    
    # see https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    tpr = TP/(TP+FN) # 1
    spc = TN/(TN+FP) # 2
    ppv = TP/(TP+FP) # 3
    npv = TN/(TN+FN) # 4
    fpr = FP/(FP+TN) # 5
    fnr = FN/(TP+FN) # 6
    fdr = FP/(TP+FP) # 7
    acc = (TP+TN)/(TP+FP+FN+TN) # 8
    
    subj[i,]=c(tpr,spc,ppv,npv,fpr,fnr,fdr,acc)
  }
  colnames(cases) = c("TP", "FP", "FN", "TN")
  colnames(subj) = c("tpr", "spc", "ppv", "npv", "fpr", "fnr", "fdr", "acc")
  
  p = list()
  p$subj=subj
  p$cases=cases
  
  p$tpr = sum(cases[,1])/(sum(cases[,1]) + sum(cases[,3]))
  p$spc = sum(cases[,4])/(sum(cases[,4]) + sum(cases[,2]))
  p$acc = (sum(cases[,1]) + sum(cases[,4]))/(sum(cases[,1]) + sum(cases[,2]) + sum(cases[,3]) + sum(cases[,4]))
  p$ppv = sum(cases[,1])/(sum(cases[,1]) + sum(cases[,2]))
  return(p)
}

perf.t <- function(adj, adj.true){
  tpr <- NULL
  spc <- NULL
  acc <- NULL
  ppv <- NULL
  d <- dim(adj)
  if(length(d)==3){
    T. <- d[3]
    for(t in 1:T.){
      tpr[t] <- perf(adj[,,t], adj.true[,,t])$tpr
      spc[t] <- perf(adj[,,t], adj.true[,,t])$spc
      acc[t] <- perf(adj[,,t], adj.true[,,t])$acc
      ppv[t] <- perf(adj[,,t], adj.true[,,t])$ppv
    }
  }
  if(length(d)==4){
    T. <- d[4]
    for(t in T.){
      tpr[t] <- perf(adj[,,,t], adj.true[,,t])$tpr
      spc[t] <- perf(adj[,,,t], adj.true[,,t])$spc
      acc[t] <- perf(adj[,,,t], adj.true[,,t])$acc
      ppv[t] <- perf(adj[,,,t], adj.true[,,t])$ppv
    }
  }
  return(list(tpr=tpr, spc=spc, acc=acc, ppv=ppv))
}


#------------------------------------------------------#
#                   Results                            #
#------------------------------------------------------#
#---- Results ----

#---- Simulations - G1 - DAG simple----
N <- 100
T. <- 220
p <- 8
adj1 <- array(0,dim=c(p,p, T.))
adj1[1,2,] = adj1[2,3,] = 1
adj1[1,2,] = adj1[1,8,] = adj1[3,4,] = adj1[5,6,] = 1
adj1[2,3,] = 1
adj1[4,5,] = adj1[6,7,] = adj1[7,8,] = 1
gplotMat(adj1[,,1],"true network", hasColMap = F)
phi0 <- 0
phi1 <- 0.9
lamb0 <- 0.001
lamb1 <- 1
v <- 0.25

g1.dag <- simu.sample(N, p, T., adj1, phi0, phi1, lamb0, lamb1, v)

#run model
tempo <- Sys.time()
g1.09<-EMVS.sample(Y=g1.dag$Y, p, T., N, Theta=0.9,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
g1.05<-EMVS.sample(Y=g1.dag$Y, p, T., N, Theta=0.5, beta0=g1.09$beta, lamb1.0=2, c=c(0.001,0.01), estlamb=T)
g1.01<-EMVS.sample(Y=g1.dag$Y, p, T., N, Theta=0.1, beta0=g1.05$beta, lamb1.0=2, c=c(0.001,0.01), estlamb=T)
Sys.time() - tempo

#---- Connections Analysis ----
g1.adj <- get.adj(g1.01$pstar, limiar=0.50, limiar.t=0.5)
poda1.adj <- poda(g1.adj$Adj.t, g1.dag$Y,g1.01$phi1, g1.01$lamb1, g1.01$c*g1.01$lamb1, g1.01$nustar)


#-- Plots

amosN <- sample(1:N, 90)
amosN1 <- sample(1:N, 50)
amosN2 <- sample(1:N, 30)
amosN3 <- sample(1:N, 10)

nameN <- c("N = 90", 'N = 50','N = 20','N = 10')

p1g1<-as.ggplot(~plot(graph_from_adjacency_matrix(adj1[,,1]), edge.arrow.size=0.5,edge.curved=0,vertex.size=10,
            edge.lty=1, edge.color="black",layout=layout.circle,vertex.frame.color="black",
            vertex.color="white", main="G1"))
p2g1<- gplotMat(adj1,"True Network", hasColMap = F)

g1.adj.N <- get.adj(g1.01$pstar[,,,amosN], limiar=0.50, limiar.t=0.5)
g1.adj.N1 <- get.adj(g1.01$pstar[,,,amosN1], limiar=0.50, limiar.t=0.5)
g1.adj.N2 <- get.adj(g1.01$pstar[,,,amosN2], limiar=0.50, limiar.t=0.5)
g1.adj.N3 <- get.adj(g1.01$pstar[,,,amosN3], limiar=0.50, limiar.t=0.5)
g1.adj <- get.adj(g1.01$pstar, limiar=0.50, limiar.t=0.5)

poda1.adj.N <- poda1.adj$Adj[,,amosN]
poda1.adj.N1 <- poda1.adj$Adj[,,amosN1]
poda1.adj.N2 <- poda1.adj$Adj[,,amosN2]
poda1.adj.N3 <- poda1.adj$Adj[,,amosN3]

stats.g1.N <- binom.nettest(poda1.adj.N, alter = "greater", fdr = 0.05)
stats.g1.N1 <- binom.nettest(poda1.adj.N1, alter = "greater", fdr = 0.05)
stats.g1.N2 <- binom.nettest(poda1.adj.N2, alter = "greater", fdr = 0.05)
stats.g1.N3 <- binom.nettest(poda1.adj.N3, alter = "greater", fdr = 0.05)

stats.g1.t <- binom.nettest.t(poda1.adj$Adj.t, alter = "greater", fdr = 0.05)
perf.g1.t <- perf.t((stats.g1.t$adj.t>=0.3)*1, adj1)
stats.g1.t.N1 <- binom.nettest.t(poda1.adj$Adj.t[,,,amosN1], alter = "greater", fdr = 0.05)
perf.g1.t.N1 <- perf.t((stats.g1.t.N1$adj.t>=0.3)*1, adj1)
stats.g1.t.N2 <- binom.nettest.t(poda1.adj$Adj.t[,,,amosN2], alter = "greater", fdr = 0.05)
perf.g1.t.N2 <- perf.t((stats.g1.t.N2$adj.t>=0.3)*1, adj1)
stats.g1.t.N3 <- binom.nettest.t(poda1.adj$Adj.t[,,,amosN3], alter = "greater", fdr = 0.05)
perf.g1.t.N3 <- perf.t((stats.g1.t.N3$adj.t>=0.3)*1, adj1)

data.perf<- data.frame(Proportion=c(perf.g1.t$tpr,perf.g1.t$spc,perf.g1.t$acc,
                                    perf.g1.t.N1$tpr,perf.g1.t.N1$spc,perf.g1.t.N1$acc,
                                    perf.g1.t.N2$tpr,perf.g1.t.N2$spc,perf.g1.t.N2$acc,
                                    perf.g1.t.N3$tpr,perf.g1.t.N3$spc,perf.g1.t.N3$acc),
                       Stat=c(rep("Sensitivity",220), rep("Specificity",220),
                              rep("Accuracy",220)),
                       Time=1:220,
                       Sample=c(rep("N = 90", 660),rep("N = 50", 660),
                                rep("N = 30", 660),rep("N = 10", 660)))
p3g1<-ggplot(data.perf, aes(Time, Proportion, color = Stat))+
  geom_point() + geom_line() +
  facet_wrap(~Sample)+ labs(color='') 
 
pt1 = gplotMat(stats.g1.N$adj, 'N = 90', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt2 = gplotMat(stats.g1.N1$adj, 'N = 50', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt3 = gplotMat(stats.g1.N2$adj, 'N = 30', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt4 = gplotMat(stats.g1.N3$adj, 'N = 10', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)

Top1 <- plot_grid(p2g1,  nrow=1, ncol=1, rel_widths = c(1,2)) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Mid <- plot_grid(pt4, pt3, pt2, pt1,  nrow=2, ncol=2) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Top <- plot_grid(Top1,Mid, nrow=1, ncol=2, rel_widths = c(1,2),labels = c("A","B"))
Bott <- plot_grid(p3g1,  nrow=1, ncol=1,labels = c("C")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
G1.plot <-plot_grid(Top, Bott, ncol = 1, nrow = 2, rel_heights = c(0.5,0.6),
          vjust = 1, hjust = -0.2)


#---- Simulations - G2 - DAG complicated ----
N <- 90
T. <- 220
p <- 8
adj2 <- array(0,dim=c(p,p, T.))
adj2[1,2,1:120] = adj2[1,8,1:120] = 1
adj2[2,3,1:120] = 1
adj2[4,5,1:120] = adj2[6,7,1:120] = adj2[7,8,1:120]= 1
gplotMat(adj2[,,1],"true network", hasColMap = F)

adj2[1,2,121:220] = adj2[1,8,121:220] = adj2[2,3,121:220] = 1
adj2[4,5,121:220] = adj2[6,7,121:220] = adj2[7,8,121:220] = 1
adj2[3,4,121:220] = adj2[5,6,121:220] = 1  
gplotMat(adj2[,,121],"true network", hasColMap = F)

phi0 <- 0
phi1 <- 0.9
lamb0 <- 0.001
lamb1 <- 1
v <- 0.25

g2.dag <- simu.sample(N, p, T., adj2, phi0, phi1, lamb0, lamb1, v)

tempo <- Sys.time()
g2.09<-EMVS.sample(Y=g2.dag$Y, p, T., N, Theta=0.9,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
g2.05<-EMVS.sample(Y=g2.dag$Y, p, T., N, Theta=0.5, beta0=g2.09$beta,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
g2.01<-EMVS.sample(Y=g2.dag$Y, p, T., N, Theta=0.1, beta0=g2.05$beta,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
Sys.time() - tempo


#---- Connections Analysis ----
g2.adj <- get.adj(g2.01$pstar, limiar=0.50, limiar.t=0.5)
poda2.adj <- poda(g2.adj$Adj.t, g2.dag$Y,g2.01$phi1, g2.01$lamb1, g2.01$c*g2.01$lamb1, g2.01$nustar)


#-- Plots

l<-1.6
h<-0.2
p1g2 <- gplotMat(adj2[,,1],"True Network \nt = 1:120", hasColMap = F, titleTextSize = 10, axisTextSize=10, textSize = 10)+
  theme(plot.margin = unit(c(h, l, h, 0.2), "cm"))
p2g2 <- gplotMat(adj2[,,121],"True Network \nt = 121:220", hasColMap = F, titleTextSize = 10, axisTextSize=10, textSize = 10)+
  theme(plot.margin = unit(c(h, l, h, 0.2), "cm"))

amosN <- sample(1:N, 90)
amosN1 <- sample(1:N, 50)
amosN2 <- sample(1:N, 30)
amosN3 <- sample(1:N, 10)

nameN <- c("N = 90", 'N = 50','N = 30','N = 10')

g2.adj <- get.adj(g2.01$pstar , limiar=0.50, limiar.t=0.5)

auxsum.t1 <- array(dim=c(p,p,N))
auxadj.t1 <- array(dim=c(p,p,N))
auxsum.t2 <- array(dim=c(p,p,N))
auxadj.t2 <- array(dim=c(p,p,N))
for(s in 1:N){
  auxsum.t1[,,s] <- apply(poda2.adj$Adj.t[,,2:120,s], c(1,2),sum)
  auxadj.t1[,,s] <- (auxsum.t1[,,s]>=(119*0.5))*1
  auxsum.t2[,,s] <- apply(poda2.adj$Adj.t[,,121:220,s], c(1,2),sum)
  auxadj.t2[,,s] <- (auxsum.t2[,,s]>=(120*0.5))*1
}

poda2.adj.t1 <- auxadj.t1[,,amosN]
poda2.adj.t2 <- auxadj.t2[,,amosN]
poda2.adj.t1.N1 <- auxadj.t1[,,amosN1]
poda2.adj.t2.N1 <- auxadj.t2[,,amosN1]
poda2.adj.t1.N2 <- auxadj.t1[,,amosN2]
poda2.adj.t2.N2 <- auxadj.t2[,,amosN2]
poda2.adj.t1.N3 <- auxadj.t1[,,amosN3]
poda2.adj.t2.N3 <- auxadj.t2[,,amosN3]

stats.g2.t1.N <- binom.nettest(poda2.adj.t1, alter = "greater", fdr = 0.05)
stats.g2.t1.N1 <- binom.nettest(poda2.adj.t1.N1, alter = "greater", fdr = 0.05)
stats.g2.t1.N2 <- binom.nettest(poda2.adj.t1.N2, alter = "greater", fdr = 0.05)
stats.g2.t1.N3 <- binom.nettest(poda2.adj.t1.N3, alter = "greater", fdr = 0.05)
stats.g2.t2.N <- binom.nettest(poda2.adj.t2, alter = "greater", fdr = 0.05)
stats.g2.t2.N1 <- binom.nettest(poda2.adj.t2.N1, alter = "greater", fdr = 0.05)
stats.g2.t2.N2 <- binom.nettest(poda2.adj.t2.N2, alter = "greater", fdr = 0.05)
stats.g2.t2.N3 <- binom.nettest(poda2.adj.t2.N3, alter = "greater", fdr = 0.05)

stats.g2.t <- binom.nettest.t(poda2.adj$Adj.t, alter = "greater", fdr = 0.05)
perf.g2.t <- perf.t((stats.g2.t$adj.t>=0.3)*1, adj2)
stats.g2.t.N1 <- binom.nettest.t(poda2.adj$Adj.t[,,,amosN1], alter = "greater", fdr = 0.05)
perf.g2.t.N1 <- perf.t((stats.g2.t.N1$adj.t>=0.3)*1, adj2)
stats.g2.t.N2 <- binom.nettest.t(poda2.adj$Adj.t[,,,amosN2], alter = "greater", fdr = 0.05)
perf.g2.t.N2 <- perf.t((stats.g2.t.N2$adj.t>=0.3)*1, adj2)
stats.g2.t.N3 <- binom.nettest.t(poda2.adj$Adj.t[,,,amosN3], alter = "greater", fdr = 0.05)
perf.g2.t.N3 <- perf.t((stats.g2.t.N3$adj.t>=0.3)*1, adj2)

data.perf<- data.frame(Proportion=c(perf.g2.t$tpr,perf.g2.t$spc,perf.g2.t$acc,
                                    perf.g2.t.N1$tpr,perf.g2.t.N1$spc,perf.g2.t.N1$acc,
                                    perf.g2.t.N2$tpr,perf.g2.t.N2$spc,perf.g2.t.N2$acc,
                                    perf.g2.t.N3$tpr,perf.g2.t.N3$spc,perf.g2.t.N3$acc),
                       Stat=c(rep("Sensitivity",220), rep("Specificity",220),
                              rep("Accuracy",220)),
                       Time=1:220,
                       Sample=c(rep("N = 90", 660),rep("N = 50", 660),
                                rep("N = 30", 660),rep("N = 10", 660)))
p3g2<-ggplot(data.perf, aes(Time, Proportion, color = Stat))+
  geom_point() + geom_line() +
  facet_wrap(~Sample)+ labs(color='') 

pt1.1 = gplotMat(stats.g2.t1.N$adj, 'N = 90 \nt = 1:120', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt2.1 = gplotMat(stats.g2.t1.N1$adj, 'N = 50 \nt = 1:120', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt3.1 = gplotMat(stats.g2.t1.N2$adj, 'N = 30 \nt = 1:120', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt4.1 = gplotMat(stats.g2.t1.N3$adj, 'N = 10 \nt = 1:120', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt1.2 = gplotMat(stats.g2.t2.N$adj, 'N = 90 \nt = 121:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt2.2 = gplotMat(stats.g2.t2.N1$adj, 'N = 50 \nt = 121:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt3.2 = gplotMat(stats.g2.t2.N2$adj, 'N = 30 \nt = 121:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt4.2 = gplotMat(stats.g2.t2.N3$adj, 'N = 10 \nt = 121:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)

Top1 <- plot_grid(p1g2, pt4.1, pt3.1, pt2.1, pt1.1,  nrow=1, ncol=5,labels = c("A","B")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Mid <- plot_grid(p2g2,pt4.2, pt3.2, pt2.2, pt1.2,  nrow=1, ncol=5,labels = c("C","D")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Top <- plot_grid(Top1,Mid, nrow=2, ncol=1)
Bott <- plot_grid(p3g2,  nrow=1, ncol=1,labels = c("E")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
G2.plot <-plot_grid(Top, Bott, ncol = 1, nrow = 2,rel_widths = c(1,1),
                    vjust = 1, hjust = -0.2)



#---- Simulations - G3 - DAG Complex----
N <- 90
T. <- 220
p <- 8
adj3 <- array(0,dim=c(p,p, T.))
adj3[1,2,1:75] = adj3[1,8,1:75] = 1
adj3[2,3,1:75] = 1
adj3[4,5,1:75] = adj3[6,7,1:75] = adj3[7,8,1:75]= 1
gplotMat(adj3[,,1],"true network", hasColMap = F)

adj3[1,2,76:150] = adj3[1,8,76:150] = adj3[2,3,76:150] = 1
adj3[4,5,76:150] = adj3[6,7,76:150] = adj3[7,8,76:150] = 1
adj3[3,4,76:150] = adj3[5,6,76:150] = 1  
gplotMat(adj3[,,76],"true network", hasColMap = F)

adj3[1,2,151:220] = adj3[2,3,151:220] = 1
adj3[6,7,151:220] = adj3[7,8,151:220] = 1
adj3[2,6,151:220] = adj3[3,7,151:220] = 1  
gplotMat(adj3[,,151],"true network", hasColMap = F)


phi0 <- 0
phi1 <- 0.9
lamb0 <- 0.001
lamb1 <- 1
v <- 0.25

g3.dag <- simu.sample(N, p, T., adj3, phi0, phi1, lamb0, lamb1, v)

tempo <- Sys.time()
g3.09<-EMVS.sample(Y=g3.dag$Y, p, T., N, Theta=0.9,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
g3.05<-EMVS.sample(Y=g3.dag$Y, p, T., N, Theta=0.5, beta0=g3.09$beta,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
g3.01<-EMVS.sample(Y=g3.dag$Y, p, T., N, Theta=0.1, beta0=g3.05$beta,lamb1.0=2, c=c(0.001,0.01), estlamb=T)
Sys.time() - tempo

#---- Connections Analysis ----
g3.adj <- get.adj(g3.01$pstar, limiar=0.50, limiar.t=0.5)
poda3.adj <- poda(g3.adj$Adj.t, g3.dag$Y,g3.01$phi1, g3.01$lamb1, g3.01$c*g3.01$lamb1, g3.01$nustar)


#-- Plots
l<-1.6
h<-0.2
p1g3 <- gplotMat(adj3[,,1],"true network\nt=1:75", hasColMap = F, titleTextSize = 10, axisTextSize=10, textSize = 10)+
  theme(plot.margin = unit(c(h, l, h, 0.2), "cm"))
p2g3 <- gplotMat(adj3[,,76],"true network\nt=76:150", hasColMap = F, titleTextSize = 10, axisTextSize=10, textSize = 10)+
  theme(plot.margin = unit(c(h, l, h, 0.2), "cm"))
p3g3 <- gplotMat(adj3[,,151],"true network\nt=151:220",hasColMap = F, titleTextSize = 10, axisTextSize=10, textSize = 10)+
  theme(plot.margin = unit(c(h, l, h, 0.2), "cm"))

amosN <- sample(1:N, 90)
amosN1 <- sample(1:N, 50)
amosN2 <- sample(1:N, 30)
amosN3 <- sample(1:N, 10)


auxsum.t1 <- array(dim=c(p,p,N))
auxadj.t1 <- array(dim=c(p,p,N))
auxsum.t2 <- array(dim=c(p,p,N))
auxadj.t2 <- array(dim=c(p,p,N))
auxsum.t3 <- array(dim=c(p,p,N))
auxadj.t3 <- array(dim=c(p,p,N))
for(s in 1:N){
  auxsum.t1[,,s] <- apply(poda3.adj$Adj.t[,,2:75,s], c(1,2),sum)
  auxadj.t1[,,s] <- (auxsum.t1[,,s]>=(74*0.5))*1
  auxsum.t2[,,s] <- apply(poda3.adj$Adj.t[,,76:150,s], c(1,2),sum)
  auxadj.t2[,,s] <- (auxsum.t2[,,s]>=(75*0.5))*1
  auxsum.t3[,,s] <- apply(poda3.adj$Adj.t[,,151:220,s], c(1,2),sum)
  auxadj.t3[,,s] <- (auxsum.t3[,,s]>=(70*0.5))*1
}

poda3.adj.t1 <- auxadj.t1[,,amosN]
poda3.adj.t2 <- auxadj.t2[,,amosN]
poda3.adj.t3 <- auxadj.t3[,,amosN]
poda3.adj.t1.N1 <- auxadj.t1[,,amosN1]
poda3.adj.t2.N1 <- auxadj.t2[,,amosN1]
poda3.adj.t3.N1 <- auxadj.t3[,,amosN1]
poda3.adj.t1.N2 <- auxadj.t1[,,amosN2]
poda3.adj.t2.N2 <- auxadj.t2[,,amosN2]
poda3.adj.t3.N2 <- auxadj.t3[,,amosN2]
poda3.adj.t1.N3 <- auxadj.t1[,,amosN3]
poda3.adj.t2.N3 <- auxadj.t2[,,amosN3]
poda3.adj.t3.N3 <- auxadj.t3[,,amosN3]

stats.g3.t1.N <- binom.nettest(poda3.adj.t1, alter = "greater", fdr = 0.05)
stats.g3.t1.N1 <- binom.nettest(poda3.adj.t1.N1, alter = "greater", fdr = 0.05)
stats.g3.t1.N2 <- binom.nettest(poda3.adj.t1.N2, alter = "greater", fdr = 0.05)
stats.g3.t1.N3 <- binom.nettest(poda3.adj.t1.N3, alter = "greater", fdr = 0.05)
stats.g3.t2.N <- binom.nettest(poda3.adj.t2, alter = "greater", fdr = 0.05)
stats.g3.t2.N1 <- binom.nettest(poda3.adj.t2.N1, alter = "greater", fdr = 0.05)
stats.g3.t2.N2 <- binom.nettest(poda3.adj.t2.N2, alter = "greater", fdr = 0.05)
stats.g3.t2.N3 <- binom.nettest(poda3.adj.t2.N3, alter = "greater", fdr = 0.05)
stats.g3.t3.N <- binom.nettest(poda3.adj.t3, alter = "greater", fdr = 0.05)
stats.g3.t3.N1 <- binom.nettest(poda3.adj.t3.N1, alter = "greater", fdr = 0.05)
stats.g3.t3.N2 <- binom.nettest(poda3.adj.t3.N2, alter = "greater", fdr = 0.05)
stats.g3.t3.N3 <- binom.nettest(poda3.adj.t3.N3, alter = "greater", fdr = 0.05)

stats.g3.t <- binom.nettest.t(poda3.adj$Adj.t, alter = "greater", fdr = 0.05)
perf.g3.t <- perf.t((stats.g3.t$adj.t>=0.3)*1, adj3)
stats.g3.t.N1 <- binom.nettest.t(poda3.adj$Adj.t[,,,amosN1], alter = "greater", fdr = 0.05)
perf.g3.t.N1 <- perf.t((stats.g3.t.N1$adj.t>=0.3)*1, adj3)
stats.g3.t.N2 <- binom.nettest.t(poda3.adj$Adj.t[,,,amosN2], alter = "greater", fdr = 0.05)
perf.g3.t.N2 <- perf.t((stats.g3.t.N2$adj.t>=0.3)*1, adj3)
stats.g3.t.N3 <- binom.nettest.t(poda3.adj$Adj.t[,,,amosN3], alter = "greater", fdr = 0.05)
perf.g3.t.N3 <- perf.t((stats.g3.t.N3$adj.t>=0.3)*1, adj3)

data.perf<- data.frame(Proportion=c(perf.g3.t$tpr,perf.g3.t$spc,perf.g3.t$acc,
                                    perf.g3.t.N1$tpr,perf.g3.t.N1$spc,perf.g3.t.N1$acc,
                                    perf.g3.t.N2$tpr,perf.g3.t.N2$spc,perf.g3.t.N2$acc,
                                    perf.g3.t.N3$tpr,perf.g3.t.N3$spc,perf.g3.t.N3$acc),
                       Stat=c(rep("Sensitivity",220), rep("Specificity",220),
                              rep("Accuracy",220)),
                       Time=1:220,
                       Sample=c(rep("N = 90", 660),rep("N = 50", 660),
                                rep("N = 30", 660),rep("N = 10", 660)))
p4g3<-ggplot(data.perf, aes(Time, Proportion, color = Stat))+
  geom_point() + geom_line() +
  facet_wrap(~Sample)+ labs(color='') 

pt1.1 = gplotMat(stats.g3.t1.N$adj, 'N = 90 \nt = 1:75', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt2.1 = gplotMat(stats.g3.t1.N1$adj, 'N = 50 \nt = 1:75', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt3.1 = gplotMat(stats.g3.t1.N2$adj, 'N = 30 \nt = 1:75', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt4.1 = gplotMat(stats.g3.t1.N3$adj, 'N = 10 \nt = 1:75', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt1.2 = gplotMat(stats.g3.t2.N$adj, 'N = 90 \nt = 76:150', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt2.2 = gplotMat(stats.g3.t2.N1$adj, 'N = 50 \nt = 76:150', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt3.2 = gplotMat(stats.g3.t2.N2$adj, 'N = 30 \nt = 76:150', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt4.2 = gplotMat(stats.g3.t2.N3$adj, 'N = 10 \nt = 76:150', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt1.3 = gplotMat(stats.g3.t3.N$adj, 'N = 90 \nt = 151:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt2.3 = gplotMat(stats.g3.t3.N1$adj, 'N = 50 \nt = 151:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt3.3 = gplotMat(stats.g3.t3.N2$adj, 'N = 30 \nt = 151:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)
pt4.3 = gplotMat(stats.g3.t3.N3$adj, 'N = 10 \nt = 151:220', '%', titleTextSize = 10, axisTextSize=10, textSize = 10, barWidth = 0.2)

Top1 <- plot_grid(p1g3, pt4.1, pt3.1, pt2.1, pt1.1,  nrow=1, ncol=5,labels = c("A","B")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Mid1 <- plot_grid(p2g3,pt4.2, pt3.2, pt2.2, pt1.2,  nrow=1, ncol=5,labels = c("C","D")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Mid2 <- plot_grid(p3g3,pt4.3, pt3.3, pt2.3, pt1.3,  nrow=1, ncol=5,labels = c("E","F")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
Top <- plot_grid(Top1,Mid1, Mid2, nrow=3, ncol=1)
Bott <- plot_grid(p4g3,  nrow=1, ncol=1,labels = c("G")) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
G3.plot <-plot_grid(Top, Bott, ncol = 1, nrow = 2,rel_heights = c(1,0.5),
                    vjust = 1, hjust = -0.2)


