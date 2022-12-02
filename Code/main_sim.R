## This code generates the simulated data and estimates the center effects using prognostic scores. 
## This code can be used to reproduce Tables 1 and 2 in the main manuscript and Table S1 in the Supporting Information.

library(MASS)
library(survival)
source("Code/aux.R")

## (1) generates the survival outcomes from a Cox model 
set.seed(1234)
randoms = sample(1:100, 100)
phi = 0.015 - 0.0001*randoms # generate the center effects 
tau = 1.5
time.list = c(5, 10)

cox_results = list()
for(ii in 1:500){
  set.seed(ii)

  J = 100 # the number of centers 
  nJ = 200 # the number of subjects per each center
  N = J*nJ # total number of subjects in the simulated data
  
  # generate baseline covariates distribution 
  X1 = X2 = X3 = X4 = X5 = matrix(0, nrow = J, ncol = nJ)
  T0 = base.T0 = Delta = C0 = Censor = Y0 = matrix(0, nrow = J, ncol = nJ);
  betas = c(1, 1, 1, 1, 1)

  for(j in 1:J){
    prob = c(rep(0.1-0.001*(j-1), 3), rep(0.1, 4), rep(0.1+0.001*(j-1), 3))
    X1[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X2[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X3[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X4[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X5[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    
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


###############################################
## (2) generates the survival outcomes from an additive hazard model
set.seed(1234)
randoms = sample(1:100, 100)
time.list = c(5,10)

additive_results = list()
for(ii in 1:500){
  set.seed(ii)
  ##### 
  J = 100 # the number of centers 
  nJ = 200 # the number of subjects per each center
  N = J*nJ # total number of subjects in the simulated data
  
  # generate cneter-specific baseline covariates
  X1 = X2 = X3 = X4 = X5 = matrix(0, nrow = J, ncol = nJ)
  for(j in 1:J){
    prob = c(rep(0.1-0.001*(j-1), 3), rep(0.1, 4), rep(0.1+0.001*(j-1), 3)) 
    X1[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X2[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X3[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X4[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
    X5[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
  }
  
  T0 = base.T0 = Delta = C0 = Censor = Y0 = matrix(0, nrow = J, ncol = nJ);
  betas = c(-0.1, -0.1, -0.1, -0.1, -0.1)/5
  inverse.Lambda = function(input, X, coeff.Lambda){
    (1/(coeff.Lambda[3]))*(-(coeff.Lambda[2] + as.numeric(X %*% betas))/2 +  
                             sqrt(coeff.Lambda[3]*input + (coeff.Lambda[2] + as.numeric(X %*% betas))^2/4 -coeff.Lambda[1]*coeff.Lambda[3]))
  }
  coeff = cbind(rep(0, 100), 0.05 + 0.0003*randoms, 0.003)
  
  # generate survival outcomes from an additive hazards model
  for(j in 1:J){
    W1 = runif(nJ, 0, 1)
    T0[j,] = inverse.Lambda(-log(W1), cbind(X1[j,], X2[j,], X3[j,], X4[j,], X5[j,]), coeff[j,])
    C0[j,] = runif(nJ, 1, 30) # random censoring 
    Delta[j,] = ifelse(is.na(T0[j,]), 0, T0[j,] <= C0[j,]) # censoring indicator
    Y0[j,] = ifelse(is.na(T0[j,]), C0[j,], pmin(T0[j,], C0[j,])) # observed failure times
  }
  
  sim1 = data.frame(timeto = as.numeric(t(Y0)), X1 = as.numeric(t(X1)), 
                    X2 = as.numeric(t(X2)), X3 = as.numeric(t(X3)),
                    X4 = as.numeric(t(X4)), X5 = as.numeric(t(X5)),
                    center = rep(1:J, each = nJ),
                    id = 1:N, Delta = as.numeric(t(Delta)))
  sim1 = sim1[order(sim1$timeto),]
  
  ## fit a stratified Cox proportional hazards model
  fit.cox = coxph( Surv(timeto, Delta) ~ X1 + X2 + X3 + X4 + X5 + strata(center), data = sim1) 
  sim1$scores.cox = predict(fit.cox, reference = "sample")
  sim1$scores.additive = additive.score.sim(dat = sim1)
  
  
  count = atrisk =  Lambda =tmp.S = matrix(0, J, nJ); 
  S.cox = S.cox.additive.combine = S.additive.cox.combine =  S.naive = S.additive = matrix(0, J, length(time.list))
  S.cox[,1] = 1; S.additive[,1] = S.naive[,1]=   1
  S.cox.var = S.cox.additive.combine.var =   S.additive.cox.combine.var = S.naive.var = S.additive.var = matrix(0, J, length(time.list))
  
  # risk classes
  cox.class = c(quantile(sim1$scores.cox , seq(0.2, 1, 0.2))[1:4], max(sim1$scores.cox) + 0.1) 
  additive.class = c(quantile(sim1$scores.additive , seq(0.2, 1, 0.2))[1:4], max(sim1$scores.additive) + 0.1)
  
  # population level risk class probability
  all.cox.class = apply(as.matrix(sim1$scores.cox), 1, function(x) sum(x > cox.class)) + 1
  sim1$all.cox.class = all.cox.class
  all.additive.class = apply(as.matrix(sim1$scores.additive), 1, function(x) sum(x > additive.class)) + 1
  sim1$all.additive.class = all.additive.class
  
  # iterate across J centers:
  for(j in 1:J){
    tmp.center = sim1[sim1$center == j,]
    order.center = tmp.center[order(tmp.center$timeto),]
    ## prognostic score weights based on the stratified cox proportional hazard model
    table.mat = as.numeric(table(sim1$all.cox.class)) / nrow(sim1)
    center.observed.cox = c()
    for(q in 1:5){
      center.observed.cox[q] = sum(order.center$all.cox.class == q)
    }
    weights.cox = table.mat/center.observed.cox
    weights.cox = ifelse(is.infinite(weights.cox), 0, weights.cox)
    weights.cox = ifelse(is.na(weights.cox), 0, weights.cox)
    weights.cox = weights.cox*nJ/sum(table.mat[center.observed.cox > 0])
    
    ## prognostic score weights based on the stratified additive hazard model
    table.mat = as.numeric(table(sim1$all.additive.class)) / nrow(sim1)
    center.observed.additive = c()
    for(q in 1:5){
      center.observed.additive[q] = sum(order.center$all.additive.class == q)
    }
    weights.additive = table.mat/center.observed.additive
    weights.additive = ifelse(is.infinite(weights.additive), 0, weights.additive)
    weights.additive = ifelse(is.na(weights.additive), 0, weights.additive)
    weights.additive = weights.additive*nJ/sum(table.mat[center.observed.additive > 0])
    
    ## (i) using the weighted Nelson-Aalen estimator (prognostic scores from a cox model)
    order.center$weights = weights.cox[order.center$all.cox.class]
    for(t in 1:nrow(order.center)){
      if(order.center$Delta[t] == 1){
        count[j,t] = weights.cox[(order.center$all.cox.class[t])]
      }else{
        count[j,t] = 0
      }
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      atrisk[j,t] = sum(weights.cox[(order.center$all.cox.class[y.index])])
      
      # cumulative cause-specific hazards
      if(t==1){
        Lambda[j,t] = 0 + (count[j,t]/atrisk[j,t])
      }else{
        Lambda[j,t] = Lambda[j,(t-1)] + (count[j,t]/atrisk[j,t])
      }
      
      # survival function
      tmp.S[j,t] = exp(-Lambda[j,t])
    }
    
    for(k in 1:length(time.list)){
      S.cox[j,k] = tmp.S[j,which.min(abs(time.list[k] - order.center$timeto))]
      S.cox.var[j,k] = (S.cox[j,k])^2*var.Lambda(obs.Y = order.center$timeto, delta = order.center$Delta,
                                                 time.point = time.list[k], weights = order.center$weights)
    }
    
    
    
    ## (ii) using the weighted Nelson-Aalen estimator (prognostic scores from the additive hazard)
    order.center$weights = weights.additive[order.center$all.additive.class]
    for(t in 1:nrow(order.center)){ 
      if(order.center$Delta[t] == 1){
        count[j,t] = weights.additive[(order.center$all.additive.class[t])]
      }else{
        count[j,t] = 0
      }
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      atrisk[j,t] = sum(weights.additive[(order.center$all.additive.class[y.index])])
      
      # cumulative cause-specific hazards
      if(t==1){
        Lambda[j,t] = 0 + (count[j,t]/atrisk[j,t])
      }else{
        Lambda[j,t] = Lambda[j,(t-1)] + (count[j,t]/atrisk[j,t])
      }
      # survival function
      tmp.S[j,t] = exp(-Lambda[j,t])
    }
    for(k in 1:length(time.list)){
      S.additive[j,k] = tmp.S[j,which.min(abs(time.list[k] - order.center$timeto))]
      S.additive.var[j,k] = (S.additive[j,k])^2*var.Lambda(obs.Y = order.center$timeto, delta = order.center$Delta,
                                                           time.point = time.list[k], weights = order.center$weights)
    }
    
    ## (iii) using the Naive Nelson-Aalen estimator without any weighting
    order.center$weights = rep(1, nrow(order.center))
    for(t in 1:nrow(order.center)){ 
      if(order.center$Delta[t] == 1){
        count[j,t] = 1
      }else{
        count[j,t] = 0
      }
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      atrisk[j,t] = length(y.index)
      
      # cumulative cause-specific hazards
      if(t==1){
        Lambda[j,t] = 0 + (count[j,t]/atrisk[j,t])
      }else{
        Lambda[j,t] = Lambda[j,(t-1)] + (count[j,t]/atrisk[j,t])
      }
      # survival function
      tmp.S[j,t] = exp(-Lambda[j,t])
    }
    for(k in 1:length(time.list)){
      S.naive[j,k] = tmp.S[j,which.min(abs(time.list[k] - order.center$timeto))]
      S.naive.var[j,k] = (S.naive[j,k])^2*var.Lambda(obs.Y = order.center$timeto, delta = order.center$Delta,
                                                     time.point = time.list[k], weights = order.center$weights)
    }
  }
  
  additive_results[[ii]] =  list(S.cox = S.cox, S.naive = S.naive, S.additive = S.additive,
              S.cox.var = S.cox.var, S.naive.var = S.naive.var, S.additive.var = S.additive.var)
}

