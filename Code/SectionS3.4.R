## This code investigates the sensitivity to the violation of positivity assumption. 
## This code can be used to reproduce Table S5 in the Supporting Information.
library(MASS)
library(survival)
library(xtable)
source("Code/aux.R")

set.seed(1234)
randoms = sample(1:100, 100)
phi = 0.015 - 0.0001*randoms # generate the center effects 
tau = 1.5
time.list = c(5, 10)
sigma = c(0,2,4,6,8,10) # under sigma = 0, no violation of the positivity assumption 
# as a value of sigma increases, the amount of nonoverlap in the baseline covariates across centers increases 
betas = c(1, 1, 1, 1, 1)

cox_results = list()
k = 6 # k=0,1,2,3,4,5
for(ii in 1:500){
  set.seed(ii)
  J = 100 # the number of centers 
  nJ = 200 # the number of subjects per each center
  N = J*nJ # total number of subjects in the simulated data
  
  # generate baseline covariates distribution 
  X1 = X2 = X3 = X4 = X5 = matrix(0, nrow = J, ncol = nJ)
  T0 = base.T0 = Delta = C0 = Censor = Y0 = matrix(0, nrow = J, ncol = nJ);
  
  for(j in 1:J){
    X1[j,] = pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k]))
    X2[j,] = pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k]))
    X3[j,] = pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k]))
    X4[j,] = pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k]))
    X5[j,] = pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k]))
    
    W1 = runif(nJ, 0, 1)
    T0[j,] = (-log(W1) / (exp(betas[1]*X1[j,] + betas[2]*X2[j,] + 
                                betas[3]*X3[j,] + betas[4]*X4[j,] + betas[5]*X5[j,])*(phi[j])^(tau)))^(1 / tau)
    
    W2 = runif(nJ, 0, 1)
    C0[j,] = runif(nJ, 1, 100)
    Delta[j,] = T0[j,] <= C0[j,] # censoring indicator
    Y0[j,] = pmin(T0[j,], C0[j,]) # observed failure times
  }
  
  sim1 = data.frame(timeto = as.numeric(t(Y0)), X1 = as.numeric(t(X1)), 
                    X2 = as.numeric(t(X2)), X3 = as.numeric(t(X3)),
                    X4 = as.numeric(t(X4)), X5 = as.numeric(t(X5)),
                    center = rep(1:J, each = nJ),
                    id = 1:N, Delta = as.numeric(t(Delta)))
  sim1 = sim1[order(sim1$timeto),]
  
  ## fit a stratified Cox proportional model
  fit.cox = coxph( Surv(timeto, Delta) ~ X1 + X2 + X3 + X4 + X5 + strata(center), data = sim1) 
  sim1$scores.cox = predict(fit.cox, reference = "sample")
  S.cox =  S.naive = S.cox.var = S.naive.var = matrix(0, J, length(time.list))
  S.cox[,1] = 1; S.naive[,1] = 1
  
  ## risk classes
  cox.class = c(quantile(sim1$scores.cox , seq(0.2, 1, 0.2))[1:4], max(sim1$scores.cox) + 0.1) 
  all.cox.class = apply(as.matrix(sim1$scores.cox), 1, function(x) sum(x > cox.class)) + 1
  sim1$all.cox.class = all.cox.class
  
  # iterate across J centers:
  for(j in 1:J){
    tmp.center = sim1[sim1$center == j,]
    nJ = nrow(tmp.center)
    order.center = tmp.center[order(tmp.center$timeto),]
    ## cox weight
    table.mat = as.numeric(table(sim1$all.cox.class)) / nrow(sim1)
    center.observed.cox = c()
    for(q in 1:5){
      center.observed.cox[q] = sum(order.center$all.cox.class == q)
    }
    weights.cox = table.mat/center.observed.cox
    weights.cox = ifelse(is.infinite(weights.cox), 0, weights.cox)
    weights.cox = ifelse(is.na(weights.cox), 0, weights.cox)
    weights.cox = weights.cox*nJ/sum(table.mat[center.observed.cox > 0])
    
    ## (i) using the weighted Nelson-Aalen estimator (prognostic scores from a cox model)
    count = atrisk = Lambda =tmp.S = c()
    order.center$weights = weights.cox[order.center$all.cox.class]
    for(t in 1:nrow(order.center)){
      if(order.center$Delta[t] == 1){
        count[t] = weights.cox[(order.center$all.cox.class[t])]
      }else{
        count[t] = 0
      }
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      atrisk[t] = sum(weights.cox[(order.center$all.cox.class[y.index])])
      
      # cumulative cause-specific hazards
      if(t==1){
        Lambda[t] = 0 + (count[t]/atrisk[t])
      }else{
        Lambda[t] = Lambda[(t-1)] + (count[t]/atrisk[t])
      }
      
      # survival function
      tmp.S[t] = exp(-Lambda[t])
    }
    
    for(k in 1:length(time.list)){
      S.cox[j,k] = tmp.S[which.min(abs(time.list[k] - order.center$timeto))]
      S.cox.var[j,k] = (S.cox[j,k])^2*var.Lambda(obs.Y = order.center$timeto, delta = order.center$Delta,
                                                 time.point = time.list[k], weights = order.center$weights)
    }
    
    ## (ii) using the Naive Nelson-Aalen estimator without weighting
    count = atrisk = Lambda =tmp.S = c()
    order.center$weights = rep(1, nrow(order.center))
    for(t in 1:nrow(order.center)){ 
      if(order.center$Delta[t] == 1){
        count[t] = 1
      }else{
        count[t] = 0
      }
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      atrisk[t] = length(y.index)
      
      # cumulative cause-specific hazards
      if(t==1){
        Lambda[t] = 0 + (count[t]/atrisk[t])
      }else{
        Lambda[t] = Lambda[(t-1)] + (count[t]/atrisk[t])
      }
      
      # survival function
      tmp.S[t] = exp(-Lambda[t])
    }
    for(k in 1:length(time.list)){
      S.naive[j,k] = tmp.S[which.min(abs(time.list[k] - order.center$timeto))]
      S.naive.var[j,k] = (S.naive[j,k])^2*var.Lambda(obs.Y = order.center$timeto, delta = order.center$Delta,
                                                     time.point = time.list[k], weights = order.center$weights)
    }
  }
  
  cox_results[[ii]] =  list(S.cox = S.cox, S.naive = S.naive,
                            S.cox.var = S.cox.var, S.naive.var = S.naive.var)
}

## calculate true S
M = 10^4; true.S = matrix(NA, M, 100); k=6
for(ii in 1:M){
  set.seed(ii)
  
  J = 100 # the number of centers 
  nJ = 200 # the number of subjects per each center
  N = J*nJ # total number of subjects in the simulated data
 
  X1 = X2 = X3 = X4 = X5 = c()
  T0 = base.T0 = Delta = C0 = Censor = Y0 = matrix(0, nrow = J, ncol = nJ);
  
  # generate the baseline covariates for "a reference population"
  for(j in 1:J){
    X1 = c(X1, pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k])))
    X2 = c(X2, pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k])))
    X3 = c(X3, pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k])))
    X4 = c(X4, pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k])))
    X5 = c(X5, pmax(0.1+0.001*j*sigma[k], pmin(rnorm(nJ, 0.45+0.001*(j-1), 0.5), 1+0.001*j*sigma[k])))
  }
  for(j in 1:J){
    W1 = runif(N, 0, 1)
    T0 = (-log(W1) / (exp(betas[1]*X1 + betas[2]*X2 + 
                            betas[3]*X3 + betas[4]*X4 + betas[5]*X5)*(phi[j])^(tau)))^(1 / tau)
    true.S[ii,j] = mean(T0 >= 5)
  }
}

tau.true = colMeans(true.S) - mean(colMeans(true.S))

dat = cox_results
S.cox = S.naive = S.cox.var = S.naive.var = matrix(NA, length(dat), 100)
for(i in 1:length(dat)){
  S.cox[i,] = dat[[i]]$S.cox[,1]
  S.naive[i,] = dat[[i]]$S.naive[,1]
  S.cox.var[i,] = dat[[i]]$S.cox.var[,1]
  S.naive.var[i,] = dat[[i]]$S.naive.var[,1]
}

## the excess probability (tau) results
S.cox.tau = S.naive.tau = matrix(0, length(dat), 100)
for(i in 1:length(dat)){
  S.cox.tau[i,] = S.cox[i,] - mean(S.cox[i,])
  S.naive.tau[i,] = S.naive[i,] - mean(S.naive[i,])
}

S.cox.tau.var = S.naive.tau.var = 
  S.cox.tau.cr = S.naive.tau.cr = matrix(NA, length(dat), 100)
for(i in 1:length(dat)){
  for(k in 1:100){
    S.cox.tau.var[i,k] = ((100-1)/100)^2*(S.cox.var[i,k]) +
      sum((1/100)^2*(S.cox.var[i,-k])) 
    S.naive.tau.var[i,k] = ((100-1)/100)^2*(S.naive.var[i,k]) +
      sum((1/100)^2*(S.naive.var[i,-k])) 
  }
  S.cox.tau.cr[i,] = (S.cox.tau[i,] - 1.96*sqrt(S.cox.tau.var[i,]) <= tau.true & 
                        S.cox.tau[i,] + 1.96*sqrt(S.cox.tau.var[i,]) >= tau.true)
  S.naive.tau.cr[i,] = (S.naive.tau[i,] - 1.96*sqrt(S.naive.tau.var[i,]) <= tau.true & 
                          S.naive.tau[i,] + 1.96*sqrt(S.naive.tau.var[i,]) >= tau.true)
}

ind = order(tau.true)[seq(1, 100, 10)]
mat = matrix(NA, 5, 10)
colnames(mat) = c(1:10)
rownames(mat) = c("True S", "True tau", "Bias in tau", "SE", "Coverage rate")
mat[1,] = colMeans(true.S)[ind]
mat[2,] = tau.true[ind]
mat[3,] = apply(S.cox.tau, 2,mean)[ind]-tau.true[ind]
mat[4,] = apply(S.cox.tau, 2, sd)[ind]
mat[5,] = apply(S.cox.tau.cr, 2, mean)[ind]
print(xtable(mat, digits = 3))


ind = order(tau.true)[seq(1, 100, 10)]
mat = matrix(NA, 4, 10)
colnames(mat) = c(1:10)
rownames(mat) = c("True tau", "Bias in tau", "SE", "Coverage rate")
mat[1,] = tau.true[ind]
mat[2,] = apply(S.naive.tau, 2,mean)[ind]-tau.true[ind]
mat[3,] = apply(S.naive.tau, 2, sd)[ind]
mat[4,] = apply(S.naive.tau.cr, 2, mean)[ind]
print(xtable(mat, digits = 3))

print(round(summary(apply(S.cox.tau, 2,mean)-tau.true),3))
print(round(summary(apply(S.cox.tau.cr, 2, mean)),3))