################################################################################# 
###  function  q(x) used in the logistic model  ################################# 
qfun1<-function(x){ x   }
qfun2<-function(x){ x^2 } 
#################################################################################


#################################################################################
###   Generate data    ##########################################################
###   Input
###      n_big_0:   population size 
###      k:   the number of captures
###      beta0:  true value of beta

###
###   Output
###      x:   oberved covariates
###      d:   capture history
#################################################################################
dat.gen<-function(n_big_0, k, beta0)
{ xa=rnorm(n_big_0) 
  xa=matrix(xa, ncol=1)
  G1=G1_fun(xa, beta0)
  da  = rbinom(n_big_0, k, G1)    
  x=xa[da>0, ]
  d=da[da>0]
  list(x= x, d= d)
}

#################################################################################
###       main program               ############################################
#################################################################################
start.time=proc.time()    #### record the beginning time
source("abun-special.R")


beta00=c(-1,  2, -0.2)    
dim_q = 3; NN=c(200, 400, 1000);  nrep=2000;   
flnm=c('D-200-2.txt', 'D-200-3.txt', 'D-400-2.txt', 
       'D-400-3.txt', 'D-1000-2.txt', 'D-1000-3.txt')

kk=c(2, 8)

for(i in 1:2)
{ n_big_0=NN[i]
  for(j in 1:2)
  { k=kk[j]
    print(c(NN[i], k))
    beta0=beta00[1:dim_q]
    result=simu(n_big_0, k, beta0, nrep) 

    m = (i-1)*2+j
    write.table(result,  flnm[m], col.names=T, row.names=F) 
   }
}

#################### Report CPU time consumed ############################## 
end.time=proc.time()           
cputime=end.time-start.time
print(cputime)  
