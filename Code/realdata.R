##################################################
library(survival)
library(gplots)
source("Code/aux.R")

## explore the data
dat = read.csv("Data/sample.csv", header = TRUE, sep = ",")
dim(dat) # 5000 x 27
names(dat)
summary(dat$folltime) # distribution of the observed survival times
table(dat$GF) # failure indicator
dat$REC_CTR_ID = as.character(dat$REC_CTR_ID)
table(dat$REC_CTR_ID)
# delete the centers with < 25 patients
toosmall = names(table(dat$REC_CTR_ID))[which(as.numeric(table(dat$REC_CTR_ID)) < 25)] 
subdat = dat[!(dat$REC_CTR_ID %in% toosmall),]
# check the center size 
summary(as.numeric(table(subdat$REC_CTR_ID)))
center.name = names(table(subdat$REC_CTR_ID))
subdat = subdat[order(subdat$folltime),]

## fit a stratified cox proportional model
fit.cox = coxph( Surv(folltime, GF) ~ age50 + year2010 + yrs_wl + log_KDRI +
                   rec_hcv + diag_poly + diag_hyper + diag_other + diabetes + BMI + 
                   yrs_dial + female + race_Black + race_Hispanic + race_Asian + 
                   race_Other + blood_a + blood_ab + blood_b + copd + hyperten + malig + 
                   ins_private + strata(REC_CTR_ID), data = subdat) 

## fit a multivariate propensity score model 
fit.multinom = multinom(REC_CTR_ID ~ age50 + year2010 + yrs_wl + log_KDRI +
                          rec_hcv + diag_poly + diag_hyper + diag_other + diabetes + BMI + 
                          yrs_dial + female + race_Black + race_Hispanic + race_Asian + 
                          race_Other + blood_a + blood_ab + blood_b + copd + hyperten + malig + 
                          ins_private, data = subdat, MaxNWts = 10000)

time.list = 365 # evaluate the survival after t=365 (days)
J = length(center.name)
subdat$scores.cox = predict(fit.cox, reference = "sample")
# population level risk class probability
all.cox.class = apply(as.matrix(subdat$scores.cox), 1, function(x) sum(x > cox.class)) + 1
subdat$all.cox.class = all.cox.class

S.cox = S.cox.var = rep(0, J)
for(j in 1:J){
  
  tmp.center = subdat[subdat$REC_CTR_ID == center.name[j],]
  order.center = tmp.center[order(tmp.center$folltime),]
  table.mat = as.numeric(table(subdat$all.cox.class)) / nrow(subdat)
  center.observed.cox = c()
  for(q in 1:5){
    center.observed.cox[q] = sum(order.center$all.cox.class == q)
  }
  weights.cox = table.mat/center.observed.cox
  weights.cox = ifelse(is.infinite(weights.cox), 0, weights.cox)
  weights.cox = ifelse(is.na(weights.cox), 0, weights.cox)
  weights.cox = weights.cox*nrow(tmp.center)/sum(table.mat[center.observed.cox > 0])
  
  count = atrisk =  Lambda =tmp.S = c()
  
  order.center$weights = weights.cox[order.center$all.cox.class]
  for(t in 1:nrow(order.center)){
    if(order.center$GF[t] == 1){
      count[t] = weights.cox[(order.center$all.cox.class[t])]
    }else{
      count[t] = 0
    }
    
    y.index = which((order.center$folltime >= order.center$folltime[t]))
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
  
  S.cox[j] = tmp.S[which.min(abs(time.list - order.center$folltime))]
  S.cox.var[j] = (S.cox[j])^2*var.Lambda(obs.Y = order.center$folltime, delta = order.center$GF,
                                         time.point = time.list, weights = order.center$weights)
  
}

## estimate the center effect (the excess survival)
S.cox.tau.var = S.cox.tau = c()
for(j in 1:length(center.name)){
  S.cox.tau[j] = S.cox[j] - mean(S.cox)
  S.cox.tau.var[j] = ((J-1)/J)^2*(S.cox.var[j]) + sum((1/J)^2*(S.cox.var[-j])) 
}



pdf("Figure/S_cox_tau_oneyear.pdf", width = 16, height = 8)
par(mfrow = c(1,1), cex.lab = 2, cex.axis = 2, 
    mar=c(4,6,3,1), tcl = 0.5, oma = c(2, 2, 3, 3), xpd = FALSE, cex.main = 2)
plot(1:length(center.name), 
     S.cox.tau[order(S.cox.tau)], col = "orangered", 
     cex = 0.8, ylim = c(-0.3, 0.2), xlim = c(0, J+1),
     ylab = expression(paste("Estimated ", hat(tau)[j], "(t) at t = 365 (one-year)")), 
     xlab = "Ordered Centers",
     main = expression(paste("Excess survival probability ", tau[j], "(t)")), type = "p", pch = 19)
abline(h = 0, col = "grey")
plotCI(S.cox.tau[order(S.cox.tau)], y=NULL, 
       ui = S.cox.tau[order(S.cox.tau)]+1.96*sqrt(S.cox.tau.var[order(S.cox.tau)]), 
       li = S.cox.tau[order(S.cox.tau)]-1.96*sqrt(S.cox.tau.var[order(S.cox.tau)]),
       err="y", pch=20, lty=1, col = "black", add= TRUE, cex = 0.1, lwd = 0.05)
dev.off()
