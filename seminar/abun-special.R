#################################################################################
####   Function used to calculate     log( gamma(N+1)/gamma(N-n+1) )  
digam<-function(n_big, n_small)
{ out=0
  if(n_small>0&n_big>=n_small)  out=sum(log(n_big+1-c(1:n_small)))    
  out
}

#################################################################################
###   minimize (minus)likelihood of the first part:    
###           - log( gamma(N+1)/gamma(N-n+1) )  - (N-n) log(alpha)   
like_part1<-function(n_small, alpha)
{ fgamma <- function(n_big)
  {   -digam(n_big, n_small)-(n_big-n_small)*log(alpha)   }   
  optimize(fgamma, interval=c(n_small, 100*n_small), tol=0.01)
}

#################################################################################
###    ELR function with  t  a vector: 
###        min_{lambda}   - sum_i  log(1+lambda t_i) 
elr_fun<-function(t)
{ emplik<-function(lam)
  { eps=1e-5;  z=1+lam*t; z1=z[z<eps]/eps
    -sum(log(z[z>eps]))-sum(2*z1-0.5*z1^2+log(eps)-1.5)
   }
  optimize(f=emplik, interval=c(-1000, 1000))
}
#################################################################################
###   Function used to calculate    pi(t) = exp(t)/(1+exp(t))
pi_fun<-function(t){ sign(t)/(1+exp(-abs(t))) + (1-sign(t))/2 }

#################################################################################
############# functions used in model ##########################
qf_tmp <- function(t, dm)
{  if(dm==2){ out =  c(1, qfun1(t))
   }else{     out =  c(1, qfun1(t), qfun2(t))}
   out 
}
################################################################################# 
qf <- function(t, dm)
{  t=as.matrix(t)
   n_small = dim(t)[1]
   out = matrix(rep(0, n_small*dm), ncol=dm)
   for(i in 1:n_small)
   { out[i, ] = qf_tmp(t[i, ], dm)   }
   out  
} 
#################################################################################
###   Function used to calculate G1
G1_exponent<-function(u, beta)
{  dm=length(beta)
   qm=qf(u, dm)
   qm%*%beta  
}
#################################################################################
###   Function used to calculate G1
G1_fun<-function(u, beta)
{  dm=length(beta)
   qm=qf(u, dm)
   pi_fun(qm%*%beta) 
}

#################################################################################
like_part23<-function(x, d, k, alpha)
{  n_small = length(d)

  fun23<-function(beta)
  {   G1_ex = G1_exponent(x, beta); 
      phi_ex=0.5*(sign(G1_ex)+1)*G1_ex+log(1+exp(-abs(G1_ex)))
      phi_ex = -k*phi_ex 
      phi=exp(phi_ex)
      tt=sum(d*G1_ex+phi_ex) 
      fv = elr_fun(phi-alpha)$objective  + tt
     # print(c(beta,fv))
      -fv
   }
 
  gr<-function(beta)
  { G1 = G1_fun(x, beta);  phi  = (1-G1)^k
    elr_res=elr_fun(phi-alpha);     lam= elr_res$minimum 

    V1=rep(0, dim_q) 
    qall = qf(x, dim_q)
    for(i in 1:n_small)
    { qfxi = qall[i, ]   
      tmp = (lam*alpha-1)/(1+lam*(phi[i]-alpha)) 
     
      V1 = V1 +(tmp*k*G1[i] + d[i])*qfxi 
     } 

    V1=V1/n_small
    -V1
   }

#  beta0=rep(0, dim_q) 
   beta0=beta_init   
   lower=rep(-1000, dim_q)
   out=nlminb(beta0, fun23, gradient=gr, lower=lower, upper=-lower) 
   list(par=out$par, value=out$objective)
}
#################################################################################
###   calculate empirical likelihood function   #################################
###   Input
###      x:   a vector (dimension: n) or matrix (number of rows: n). 
###           Here n denotes the number of observations.
###      d:   vector (dimenstion: n), capture history. 
###      k:   the number of captures.
###      n_big_0:   a given value for the population size,  at which 
###           the empirical log-likelihood function is evaluated.
###
###   Output
###      out$objective:     empirical log-likelihood function at n_big_0.
#################################################################################
likelihood.null<-function(x,d, k, n_big_0)
{  n_small = length(d)
   falpha.null<-function(alpha) 
   {  tmp1=-(n_big_0-n_small)*log(alpha)-digam(n_big_0, n_small) 
      tmp2=like_part23(x, d, k, alpha)$value
      tmp1+tmp2 
    } 
   eps=1e-6
   out=optimize(falpha.null, interval=c(eps, 1-eps), tol=0.0001)
   out$objective
} 
#################################################################################
###   calculate empirical likelihood estimates   ################################
###   Input
###      x:   a vector (dimension: n) or matrix (number of rows: n). 
###           Here n denotes the number of observations.
###      d:   vector (dimenstion: n), capture history. 
###      k:   the number of captures.
###      n_big_0:   a given value for the population size,  used to 
###           calculate  empirical likelihood ratio
###
###   Output
###      lrt:      minus twice empirical log-likelihood ratio
###      n_hat:    Maximum empirical likelihood estimate of N
###      alpha_hat: Maximum empirical likelihood estimate of alpha
###      beta_hat:  Maximum empirical likelihood estimate of beta
###      n_small:  number of observations. 
###
#################################################################################
likelihood<-function(x,d, k, n_big_0)
{  n_small=length(d)

   falpha<-function(alpha) 
   {  tmp1=like_part1(n_small, alpha)$objective
      tmp2=like_part23(x, d, k, alpha)$value 
      tmp1+tmp2  
    } 

   eps=1e-6
   out=optimize(falpha, interval=c(eps, 1-eps), tol=0.0001)
   like= out$objective
   alpha_est=out$minimum 

##############  EL estimator of N  ##############
   out1=like_part1(n_small, alpha_est)
   n_hat = out1$minimum 
#################################################   
############## EL estimator of beta  ############               
   out23=like_part23(x,d, k, alpha_est) 
   beta_hat = out23$par                 
#################################################  

  like0 =likelihood.null(x,d, k, n_big_0)
  
  lrt=2*(like0-like)
  list(lrt=lrt, n_hat=n_hat, alpha_est=alpha_est, beta_hat=beta_hat, 
       n_small=n_small, like=-like)
}
#################################################################################
###   calculate conditional likelihood estimates   ##############################
###   Input
###      x:   a vector (dimension: n) or matrix (number of rows: n). 
###           Here n denotes the number of observations.
###      d:   vector (dimenstion: n), capture history. 
###      k:   the number of captures.
###
###   Output
###      n_tilde:  Maximum conditional likelihood estimate of N
###      beta_tilde:  Maximum conditional likelihood estimate of beta
###      sigma_hat:  Estimate of the asymptotic variance of n_tilde. 
#################################################################################
likelihood.condition<-function(x, d, k)
{  n_small=length(d)
 
   fun.c<-function(beta)
   {  G1_ex = G1_exponent(x, beta); 
      phi_ex=0.5*(sign(G1_ex)+1)*G1_ex+log(1+exp(-abs(G1_ex)))
      phi_ex = -k*phi_ex 
      phi=exp(phi_ex) 
      newphi=phi
      for(i in 1:n_small) 
      { if(phi[i]==1){ 
            newphi[i] = G1_ex[i]+log(k)  
        }else{ newphi[i] = log(1-phi[i]) }
        }
      tt= d*G1_ex+phi_ex 
      
      -sum(tt)+sum(newphi)
   } 
 
  gr.c<-function(beta)
  { dm=length(beta); G1 = G1_fun(x, beta);  phi  = (1-G1)^k
    
    V1=rep(0, dm)
    qall = qf(x, dm)
    for(i in 1:n_small)
    { tmp = d[i] - k*G1[i]/(1-phi[i])
      V1 = V1 +tmp*qall[i, ] 
     } 
    -V1/n_small
   }
   
   beta0=rep(0, dim_q)      
   lower=rep(-1000, dim_q)
   out=nlminb(beta0, fun.c, gradient=gr.c, lower=lower, upper=-lower) 
  
   beta_tilde=out$par
   beta_init<<- beta_tilde
   V23 = rep(0, dim_q);  V22= V23%*%t(V23)
   G1 = G1_fun(x, beta_tilde);  phi  = (1-G1)^k
   for(i in 1:n_small){ if(phi[i]>1-1e-10){ phi[i]=1-1e-10}  }
   n_tilde = sum(1/(1-phi)); phistar=sum(1/(1-phi)^2)/n_tilde

   qall = qf(x, dim_q) 
   for(i in 1:n_small)
   {   q=qall[i, ]
       aa = phi[i]/(1-phi[i])^2*k*G1[i]
       V23=V23+ (aa*q) 
       bb = d[i] - k*G1[i]/(1-phi[i]) ; 
       tmp=bb^2*(q%*%t(q))   
       V22=V22+ tmp
      }

   V23=-V23/n_tilde;    V22=V22/n_tilde 
   ev=eigen(V22, symmetric=1);  ev.d=ev$values; ev.v=ev$vectors
   for( j in 1:length(ev.d)){  
      if(abs(ev.d[j])>1e-10){ev.d[j]=1/ev.d[j]
        }else{ ev.d[j] = 1e10}    
    } 
   V22.inv=ev.v%*%diag(ev.d)%*%t(ev.v); 
   
   tmp = t(V23)%*%V22.inv%*%V23 
   sigma_hat =  phistar - 1 + tmp 
   list(n_tilde=n_tilde, beta_tilde= out$par,  sigma_hat =sigma_hat )
}



#################################################################################
###   calculate empirical likelihood confidence interval for abundance   ######## 
###   Input
###      x:   a vector (dimension: n) or matrix (number of rows: n). 
###           Here n denotes the number of observations.
###      d:   vector (dimenstion: n), capture history. 
###      k:   the number of captures.
###      crit:  critical value to determine the interval
###      like0: the global maximum empirical log-likelihood 
###      n_hat: maximum empirical likelihood estimate of population size
###
###   Output
###      interval:   vector of length 2,  consisting of the lower and upper bound
###             of the resultant EL confidence interval 
#################################################################################
CI_EL<-function(x, d, k, crit, like0, n_hat)
{ interval=c(0,0); n_small = length(d)
  c1=n_small;   c2=c1
  
  like=-likelihood.null(x,d, k, n_small)  
  if(like<crit)
  {  c2=n_hat
     while(c2-c1>1)
     { c.tmp= (c1+c2)/2 
       like=-likelihood.null(x,d, k,c.tmp)
       if(like>crit){
         c2=c.tmp
       }else{c1=c.tmp}
      #  cat('lower:', c(c1,c2,c.tmp, like, crit), '\n')
     }
  }
  interval[1] = (c1+c2)/2
  
  c1=n_hat 
  like = like0
  c3=1; temp=c1
  while(like>crit)
  { c1=temp; temp=temp+10*c3
    like=-likelihood.null(x,d, k,temp)
    c3=c3*2
   #  cat('For upper:', c(c1, like, crit, c3), '\n')
  }
  c1=n_hat
  c2=temp
  while(c2-c1>1)
  { c.tmp= (c1+c2)/2 
    like=-likelihood.null(x,d, k,c.tmp)
    if(like>crit){ c1=c.tmp
     }else{c2=c.tmp}
    # cat('upper:', c(c1,c2,c.tmp, like, crit), '\n')
  }
  interval[2] = (c1+c2)/2
  interval
}


#################################################################################
###  Simulation in one scenario
#################################################################################
simu<-function(n_big_0, k, beta0, nrep)
{     
  result=NULL 
  for(i in 1:nrep)
  {  set.seed(i*2377-500)
    out=dat.gen(n_big_0, k, beta0)
    x=out$x;  d=out$d 
    out.CL=likelihood.condition(x, d, k) 
    out.EL=likelihood( x,  d, k,  n_big_0)  

n_hat = out.EL$n_hat
like0 = out.EL$like
crit  = like0-0.5*3.841459
interval_EL=CI_EL(x, d, k, crit, like0, n_hat)


    stat=c(i, out.EL$n_small, out.EL$n_hat, out.CL$n_tilde, out.EL$alpha_est,
         out.EL$lrt, out.CL$sigma_hat, interval_EL)
    stat=round(stat,4)
    result=rbind(result, stat) 
    print(stat)
  }  
  
  colnames(result)=c('No.', 'n_small', 'n_hat', 'n_tilde', 'alpha_hat', 'lrt',
     'sigma_hat', 'EL-low', 'EL-upper')
  result
}

