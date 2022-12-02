# This code explores the impact of varied number of risk classes on the real data application results and compares the estimated excess three-year survival probabilities to the one-year survival probabilities. 
# This code can be used to reproduce Figures S2-S4 and Table S6 in the Supporting Information.

library(sas7bdat)
library(survival)
library(gplots)
library(gtools)
source("Code/aux.R")
####
dat = read.sas7bdat("Data/postkt_13oct2020.sas7bdat")
head(dat); dim(dat) # 67552 x 31
names(dat)
summary(dat$folltime)
table(dat$GF) # more than 80% subjects are censored 
dat$REC_CTR_ID = as.character(dat$REC_CTR_ID)
table(dat$REC_CTR_ID)
summary(as.numeric(table(dat$REC_CTR_ID)))
length(as.numeric(table(dat$REC_CTR_ID))) # 251
toosmall = names(table(dat$REC_CTR_ID))[which(as.numeric(table(dat$REC_CTR_ID)) < 25)] # 44
subdat = dat[!(dat$REC_CTR_ID %in% toosmall),]
## check 
summary(as.numeric(table(subdat$REC_CTR_ID)))

fit.cox = coxph(Surv(folltime, GF) ~ age50 + year2010 + yrs_wl + log_KDRI +
                   rec_hcv + diag_poly + diag_hyper + diag_other + diabetes + BMI + 
                   yrs_dial + female + race_Black + race_Hispanic + race_Asian + 
                   race_Other + blood_a + blood_ab + blood_b + copd + hyperten + malig + 
                   ins_private + strata(REC_CTR_ID), data = subdat) 


subdat$scores.cox = predict(fit.cox, reference = "sample")

## risk classes
R = 20 # the number of risk sets, R = 5, 7, 10, 15, 20
cox.class = c(quantile(subdat$scores.cox , seq(1/R, 1, 1/R))[1:(R-1)], max(subdat$scores.cox) + 0.1)

## population level risk class probability
all.cox.class = apply(as.matrix(subdat$scores.cox), 1, function(x) sum(x > cox.class)) + 1
subdat$all.cox.class = all.cox.class

# iterate across J centers:
center.names = names(table(subdat$REC_CTR_ID))
J = length(center.names)
nj = c()
time.list = 365 # one-year survival

S.cox = S.cox.var = rep(NA, J)
for(j in 1:J){
  print(j)
  tmp.center = subdat[subdat$REC_CTR_ID == center.names[j],]
  nj[j] = nrow(tmp.center) # center size
  order.center = tmp.center[order(tmp.center$folltime),]
  table.mat = as.numeric(table(subdat$all.cox.class)) / nrow(subdat)
  center.observed.cox = c()
  for(q in 1:R){
    center.observed.cox[q] = sum(order.center$all.cox.class == q)
  }
  weights.cox = table.mat/center.observed.cox
  weights.cox = ifelse(is.infinite(weights.cox), 0, weights.cox)
  weights.cox = ifelse(is.na(weights.cox), 0, weights.cox)
  weights.cox = weights.cox*nrow(order.center)/sum(table.mat[center.observed.cox > 0])
  
  count = atrisk =  Lambda =tmp.S = c()
  # iterate across observed failure times within center:
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
    # approximate the time points to specified times points at time.list 
  }
  
  S.cox[j] = tmp.S[which.min(abs(time.list - order.center$folltime))]
  S.cox.var[j] = (S.cox[j])^2*var.Lambda(obs.Y = order.center$folltime, delta = order.center$GF,
                                         time.point = time.list, weights = order.center$weights)
}

### estimate tau
S.cox.tau.var = S.cox.tau = c()
for(k in 1:207){
  S.cox.tau[k] = S.cox[k] - sum(nj*S.cox)/sum(nj)
  S.cox.tau.var[k] = (1-nj[k]/sum(nj))^2*(S.cox.var[k]) +
    sum((nj[-k]/sum(nj))^2*(S.cox.var[-k])) 
}

#result.R5 = list(S.cox.tau, S.cox.tau.var) # when R = 5
#result.R7 = list(S.cox.tau, S.cox.tau.var) # when R = 7
#result.R10 = list(S.cox.tau, S.cox.tau.var) # when R = 10
#result.R15 = list(S.cox.tau, S.cox.tau.var) # when R = 15
result.R20 = list(S.cox.tau, S.cox.tau.var) # when R = 20

pdf("Figure/S_cox_tau_oneyear_R20.pdf", width = 16, height = 8)
par(mfrow = c(1,1), cex.lab = 2, cex.axis = 2, 
    mar=c(5,6,3,1), tcl = 0.5, oma = c(1, 2, 3, 3), xpd = FALSE, cex.main = 2)
plot(1:207, 
     S.cox.tau[order(result.R5[[1]])], col = "cornflowerblue", 
     cex = 0.8, ylim = c(-0.3, 0.2), xlim = c(0, 212),
     ylab = expression(paste("Estimated ", hat(tau)[j], "(t) at t = 365 (one-year)")), 
     xlab = "Ordered centers by the estimates based on R =5",
     main = expression(paste("Excess survival probability ", tau[j], "(t) using R = 20")), type = "p", pch = 19)
abline(h = 0, col = "grey")
plotCI(x = S.cox.tau[order(result.R5[[1]])], y = NULL, 
       ui = S.cox.tau[order(result.R5[[1]])]+1.96*sqrt(S.cox.tau.var[order(result.R5[[1]])]), 
       li = S.cox.tau[order(result.R5[[1]])]-1.96*sqrt(S.cox.tau.var[order(result.R5[[1]])]),
       err="y", pch=20, lty=1, col = "black", add= TRUE, cex = 0.1, lwd = 0.05,
       xlab = "", ylab = "", labels = FALSE)
dev.off()


### Table S6 (Spearman rank correlation between the results) ###
round(cor(result.R5[[1]], result.R7[[1]], method = "spearman"),3)
round(cor(result.R5[[1]], result.R10[[1]], method = "spearman"),3)
round(cor(result.R5[[1]], result.R15[[1]], method = "spearman"),3)
round(cor(result.R5[[1]], result.R20[[1]], method = "spearman"),3)

round(cor(result.R7[[1]], result.R10[[1]], method = "spearman"),3)
round(cor(result.R7[[1]], result.R15[[1]], method = "spearman"),3)
round(cor(result.R7[[1]], result.R20[[1]], method = "spearman"),3)

round(cor(result.R10[[1]], result.R15[[1]], method = "spearman"),3)
round(cor(result.R10[[1]], result.R20[[1]], method = "spearman"),3)

round(cor(result.R15[[1]], result.R20[[1]], method = "spearman"),3)


######## Supporting Information Figure S4 ########
time.list = 365*3 ## three year

S.cox = S.cox.var = rep(NA, J)
for(j in 1:J){
  print(j)
  tmp.center = subdat[subdat$REC_CTR_ID == center.names[j],]
  nj[j] = nrow(tmp.center) # center size
  order.center = tmp.center[order(tmp.center$folltime),]
  table.mat = as.numeric(table(subdat$all.cox.class)) / nrow(subdat)
  center.observed.cox = c()
  for(q in 1:R){
    center.observed.cox[q] = sum(order.center$all.cox.class == q)
  }
  weights.cox = table.mat/center.observed.cox
  weights.cox = ifelse(is.infinite(weights.cox), 0, weights.cox)
  weights.cox = ifelse(is.na(weights.cox), 0, weights.cox)
  weights.cox = weights.cox*nrow(order.center)/sum(table.mat[center.observed.cox > 0])
  
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

### estimate the excess survival probability (tau)
S.cox.tau.var = S.cox.tau = c()
for(k in 1:207){
  S.cox.tau[k] = S.cox[k] - sum(nj*S.cox)/sum(nj)
  S.cox.tau.var[k] = (1-nj[k]/sum(nj))^2*(S.cox.var[k]) +
    sum((nj[-k]/sum(nj))^2*(S.cox.var[-k])) 
}

### save the results
result.threeyear = list(S.cox.tau, S.cox.tau.var)

pdf("Figure/S_cox_tau_threeyear.pdf", width = 16, height = 8)
par(mfrow = c(1,1), cex.lab = 2, cex.axis = 2, 
    mar=c(5,6,3,1), tcl = 0.5, oma = c(1, 2, 3, 3), xpd = FALSE, cex.main = 2)
plot(1:207, 
     S.cox.tau[order(result.R5[[1]])], col = "orangered", 
     cex = 0.8, ylim = c(-0.3, 0.2), xlim = c(0, 212),
     ylab = expression(paste("Estimated ", hat(tau)[j], "(t) at t = 1095 (three-year)")), 
     xlab = "Ordered centers by the estimated one-year survival",
     main = expression(paste("Excess three-year survival probability ", tau[j], "(t)")), type = "p", pch = 19)
abline(h = 0, col = "grey")
plotCI(x = S.cox.tau[order(result.R5[[1]])], y = NULL, 
       ui = S.cox.tau[order(result.R5[[1]])]+1.96*sqrt(S.cox.tau.var[order(result.R5[[1]])]), 
       li = S.cox.tau[order(result.R5[[1]])]-1.96*sqrt(S.cox.tau.var[order(result.R5[[1]])]),
       err="y", pch=20, lty=1, col = "black", add= TRUE, cex = 0.1, lwd = 0.05,
       xlab = "", ylab = "", labels = FALSE)
dev.off()

pdf("Figure/S_cox_tau_3yr.pdf", width = 16, height = 8)
par(mfrow = c(1,1), cex.lab = 2, cex.axis = 2, 
    mar=c(5,6,3,1), tcl = 0.5, oma = c(1, 2, 3, 3), xpd = FALSE, cex.main = 2)
plot(1:207, 
     S.cox.tau[order(result.threeyear[[1]])], col = "orangered", 
     cex = 0.8, ylim = c(-0.3, 0.2), xlim = c(0, 212),
     ylab = expression(paste("Estimated ", hat(tau)[j], "(t) at t = 1095 (three-year)")), 
     xlab = "Ordered centers by the estimated three-year survival",
     main = expression(paste("Excess three-year survival probability ", tau[j], "(t)")), type = "p", pch = 19)
abline(h = 0, col = "grey")
plotCI(x = S.cox.tau[order(result.threeyear[[1]])], y = NULL, 
       ui = S.cox.tau[order(result.threeyear[[1]])]+1.96*sqrt(S.cox.tau.var[order(result.threeyear[[1]])]), 
       li = S.cox.tau[order(result.threeyear[[1]])]-1.96*sqrt(S.cox.tau.var[order(result.threeyear[[1]])]),
       err="y", pch=20, lty=1, col = "black", add= TRUE, cex = 0.1, lwd = 0.05,
       xlab = "", ylab = "", labels = FALSE)
dev.off()