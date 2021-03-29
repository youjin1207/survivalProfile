var.Lambda = function(obs.Y, delta, time.point, weights){
  ################################################
  ## a variance estimator (Hu and Huffer, 2020) ##
  ## obs.Y : observed failure times
  ## delta: failure indicators 
  ## time.point: a given time point
  ## weights: individual-specific weights
  ################################################
  dat = data.frame(obs.Y, delta, weights)
  dat = dat[order(dat$obs.Y),]
  obs.riskset = c()
  for(t in 1:nrow(dat)){
    obs.riskset[t] = sum(dat$obs.Y >= dat$obs.Y[t])
  }
  
  if(sum(dat$obs.Y < time.point) == 0){
    return(var.Lambda = 0)
  }
  
  k = max(which(dat$obs.Y-time.point <=0))  # k : select the closest timepoint
  var.Lambda = 0
  for(i in 1:k){
    var.Lambda = var.Lambda + sum((dat$weights)^2*(dat$obs.Y >= dat$obs.Y[i]))*
      sum(dat$weights*(dat$obs.Y >= dat$obs.Y[i])*(dat$obs.Y == dat$obs.Y[i] & dat$delta == 1)) / 
      sum(dat$weights*(dat$obs.Y >= dat$obs.Y[i]))^3
  }
  return(var.Lambda = var.Lambda)
}



additive.score.sim = function(dat){
  ########################################################################
  ## deriving prognostic scores from a stratified additive hazards model ##
  ## a data frame 'dat' contains the following columns:
  ## X1, X2, X3, X4, X5 : continuous baseline covariates 
  ## timeto: observed failure times
  ## delta: observed failure indicators
  ## center: center id
  ################################################
  
  
  dat = dat[order(dat$timeto),]
  dat$X1.mean  = rep(NA, nrow(dat))
  tmp1 = aggregate(X1 ~ center, FUN = mean, data = dat)$center
  tmp2 = aggregate(X1~ center, FUN = mean, data = dat)$X1
  for(j in 1:length(tmp1)){
    dat$X1.mean = ifelse(dat$center == tmp1[j], rep(tmp2[j], nrow(dat)),
                         dat$X1.mean) 
  }
  
  dat$X2.mean  = rep(NA, nrow(dat))
  tmp1 = aggregate(X2 ~ center, FUN = mean, data = dat)$center
  tmp2 = aggregate(X2~ center, FUN = mean, data = dat)$X2
  for(j in 1:length(tmp1)){
    dat$X2.mean = ifelse(dat$center == tmp1[j], rep(tmp2[j], nrow(dat)),
                         dat$X2.mean) 
  }
  
  dat$X3.mean  = rep(NA, nrow(dat))
  tmp1 = aggregate(X3 ~ center, FUN = mean, data = dat)$center
  tmp2 = aggregate(X3~ center, FUN = mean, data = dat)$X3
  for(j in 1:length(tmp1)){
    dat$X3.mean = ifelse(dat$center == tmp1[j], rep(tmp2[j], nrow(dat)),
                         dat$X3.mean) 
  }
  
  dat$X4.mean  = rep(NA, nrow(dat))
  tmp1 = aggregate(X4 ~ center, FUN = mean, data = dat)$center
  tmp2 = aggregate(X4~ center, FUN = mean, data = dat)$X4
  for(j in 1:length(tmp1)){
    dat$X4.mean = ifelse(dat$center == tmp1[j], rep(tmp2[j], nrow(dat)),
                         dat$X4.mean) 
  }
  
  dat$X5.mean  = rep(NA, nrow(dat))
  tmp1 = aggregate(X5 ~ center, FUN = mean, data = dat)$center
  tmp2 = aggregate(X5~ center, FUN = mean, data = dat)$X5
  for(j in 1:length(tmp1)){
    dat$X5.mean = ifelse(dat$center == tmp1[j], rep(tmp2[j], nrow(dat)),
                         dat$X5.mean) 
  }
  
  sum(is.na(dat))
  ## stratified additive hazards model
  Ahat = Uhat = 0
  for(i in 1:nrow(dat)){
    dummy1 = (dat$timeto >= dat$timeto[i])
    dummy2 = (dat$timeto == dat$timeto[i] & dat$Delta == 1)
    Zmat1 = cbind(dat$X1[dummy1]-dat$X1.mean[dummy1], dat$X2[dummy1]-dat$X2.mean[dummy1],
                  dat$X3[dummy1]-dat$X3.mean[dummy1], dat$X4[dummy1]-dat$X4.mean[dummy1],
                  dat$X5[dummy1]-dat$X5.mean[dummy1])
    if(sum(dummy2) == 0){
      Zmat2 = rep(0, 5)
    }else{
      Zmat2 = cbind(dat$X1[dummy2]-dat$X1.mean[dummy2], dat$X2[dummy2]-dat$X2.mean[dummy2],
                    dat$X3[dummy2]-dat$X3.mean[dummy2], dat$X4[dummy2]-dat$X4.mean[dummy2],
                    dat$X5[dummy2]-dat$X5.mean[dummy2])
      Zmat2 = colMeans(Zmat2)
    }
    Ahat = Ahat + t(Zmat1)%*%Zmat1
    Uhat = Uhat + Zmat2
  }
  theta_stratified_hat = solve(Ahat) %*% as.matrix(Uhat)
  
  ## return the prognostic scores
  additive.score = as.numeric(t(cbind(dat$X1, dat$X2, dat$X3, dat$X4, dat$X5) %*% (theta_stratified_hat)))
  return(additive.score)
}
