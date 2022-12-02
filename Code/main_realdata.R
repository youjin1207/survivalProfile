# This code reads the real data example and estimates the center effects using the proposed methods. 
# This code can be used to reproduce Figures 2-3 in the main manuscript and Figure S1 in the Supporting Information.

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
center.name = names(table(subdat$REC_CTR_ID))

subdat = subdat[order(subdat$folltime),]

fit.cox = coxph( Surv(folltime, GF) ~ age50 + year2010 + yrs_wl + log_KDRI +
                   rec_hcv + diag_poly + diag_hyper + diag_other + diabetes + BMI + 
                   yrs_dial + female + race_Black + race_Hispanic + race_Asian + 
                   race_Other + blood_a + blood_ab + blood_b + copd + hyperten + malig + 
                   ins_private + strata(REC_CTR_ID), data = subdat) 


fit.multinom = multinom(REC_CTR_ID ~ age50 + year2010 + yrs_wl + log_KDRI +
                          rec_hcv + diag_poly + diag_hyper + diag_other + diabetes + BMI + 
                          yrs_dial + female + race_Black + race_Hispanic + race_Asian + 
                          race_Other + blood_a + blood_ab + blood_b + copd + hyperten + malig + 
                          ins_private, data = subdat, MaxNWts = 10000)


time.list = 365 # one-year survival probability
SMR.estimate =  matrix(0, length(center.name), length(time.list));

SMR.estimate = rep(NA, length(time.list))
J = 207
for(j in 1:J){
  print(j)
  main.center = subdat[subdat$REC_CTR_ID == center.name[j],]
  main.center$folltime = time.list
  main.center$GF = 1
  
  observed = mean(1-predict(fit.cox, type = "survival", newdata = main.center))
  expected = 0
  for(k in 1:J){
    print(k)
    tmp.center = main.center
    tmp.center$REC_CTR_ID = center.name[k]
    
    expected = expected + (1-predict(fit.cox, type = "survival", newdata = tmp.center))*
      predict(fit.multinom, newdata = tmp.center, "probs")[,k]
  }
  SMR.estimate[j] = observed / mean(expected)
}


real.SMR.365 = SMR.estimate
cor(real.tau.365, real.SMR.365)

# Save the bootstrap results in SMR.estimate.boot

### summarize the results ###
### ten lowest centers 
order.ind1 = order(real.tau.365)[1:10]
order.ind2 = order(real.tau.365, decreasing = TRUE)[1:10]

mat = matrix(NA, 2, 10)
rownames(mat) = c("Tau", "SMR")
colnames(mat) = c(1:10)
mat[1,] = S.cox.tau[order.ind1]
mat[2,] = real.SMR.365[order.ind1] - 1
lim = 1.2*max(mat)

error.bar = function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

sd.mat = matrix(NA, 2, 10)
rownames(sd.mat) = c("Tau", "SMR")
colnames(sd.mat) = c(1:10)
sd.mat[1,] = 1.96*sqrt(S.cox.tau.var[order.ind1])
sd.mat[2,] = 1.96*sqrt(apply(SMR.estimate.boot, 2, var)[order.ind1])

### Figure 2 (a) (main manuscript)
pdf("Figure/tencenters_oneyear_lowest.pdf", width = 10, height = 6)
par(mfrow = c(1,1), cex.lab = 1.5, cex.axis = 1.5, 
    mar=c(4,6,1,1), tcl = 0.5, oma = c(2, 2, 2, 3), xpd = FALSE, cex.main = 2)
ze_barplot <- barplot(mat , beside=T , 
                      ylim=c(-1,4) , ylab="Center effect estimates",
                      xlab = expression(paste("Centers with ten lowest ", hat(tau)[j], "(t)")),
                      col = c("dodgerblue", "gold"))
#points(c(1:10), mat[2,], pch = 19)
legend("topright", c(expression(paste(hat(tau)[j], "(t)")),
                     expression(paste("Estimated ", hat(SMR)[j], "(t) - 1"))),
       fill=c("dodgerblue" , "gold"), bty = "n", cex = 1.5)
abline(h = 0, lwd = 2)
error.bar(ze_barplot, mat, sd.mat)
dev.off()

### log scale
mat = matrix(NA, 2, 10)
rownames(mat) = c("Tau", "SMR")
colnames(mat) = c(1:10)
mat[1,] = S.cox.tau[order.ind1]
mat[2,] = log(real.SMR.365[order.ind1])
lim = 1.2*max(mat)
sd.mat = matrix(NA, 2, 10)
rownames(sd.mat) = c("Tau", "SMR")
colnames(sd.mat) = c(1:10)
sd.mat[1,] = 1.96*sqrt(S.cox.tau.var[order.ind1])
sd.mat[2,] = 1.96*sqrt(apply(SMR.estimate.boot, 2, var)[order.ind1]/(real.SMR.365[order.ind1])^2) # by Delta Method

### Figure S1 (a) (Supporting Information)
pdf("Figure/tencenters_oneyear_lowest_log.pdf", width = 10, height = 6)
par(mfrow = c(1,1), cex.lab = 1.5, cex.axis = 1.5, 
    mar=c(4,6,1,1), tcl = 0.5, oma = c(2, 2, 2, 3), xpd = FALSE, cex.main = 2)
ze_barplot <- barplot(mat , beside=T , 
                      ylim=c(-1,4) , ylab="Center effect estimates",
                      xlab = expression(paste("Centers with ten lowest ", hat(tau)[j], "(t)")),
                      col = c("dodgerblue", "gold"))
#points(c(1:10), mat[2,], pch = 19)
legend("topright", c(expression(paste(hat(tau)[j], "(t)")),
                     expression(paste("Estimated log(", hat(SMR)[j], "(t))"))),
       fill=c("dodgerblue" , "gold"), bty = "n", cex = 1.5)
abline(h = 0, lwd = 2)
error.bar(ze_barplot, mat, sd.mat)
dev.off()

### ten highest centers 
order.ind2 = order(real.tau.365, decreasing = TRUE)[3:12]
mat = matrix(NA, 2, 10)
rownames(mat) = c("Tau", "SMR")
colnames(mat) = c(1:10)
mat[1,] = S.cox.tau[order.ind2]
mat[2,] = real.SMR.365[order.ind2] - 1
sd.mat = matrix(NA, 2, 10)
rownames(sd.mat) = c("Tau", "SMR")
colnames(sd.mat) = c(1:10)
sd.mat[1,] = 1.96*sqrt(S.cox.tau.var[order.ind2])
sd.mat[2,] = 1.96*sqrt(apply(SMR.estimate.boot, 2, var)[order.ind2])

### Figure 3 (b) (main manuscript)
pdf("Figure/tencenters_oneyear_highest.pdf", width = 10, height = 6)
par(mfrow = c(1,1), cex.lab = 1.5, cex.axis = 1.5, 
    mar=c(4,6,1,1), tcl = 0.5, oma = c(2, 2, 2, 3), xpd = FALSE, cex.main = 2)
ze_barplot <- barplot(mat , beside=T , 
                      ylim=c(-1.5,0.5) , ylab="Center effect estimates",
                      xlab = expression(paste("Centers with ten highest ", hat(tau)[j], "(t)")),
                      col = c("dodgerblue", "gold"))
#points(c(1:10), mat[2,], pch = 19)
legend(x = 17, y = 0.6, c(expression(paste(hat(tau)[j], "(t)")),
                          expression(paste("Estimated ", hat(SMR)[j], "(t) - 1"))),
       fill=c("dodgerblue" , "gold"), bty = "n", cex = 1.5)
abline(h = 0, lwd = 2)
error.bar(ze_barplot, mat, sd.mat)
dev.off()

### log scale 
mat = matrix(NA, 2, 10)
rownames(mat) = c("Tau", "SMR")
colnames(mat) = c(1:10)
mat[1,] = S.cox.tau[order.ind2]
mat[2,] = log(real.SMR.365[order.ind2])
lim = 1.2*max(mat)
sd.mat = matrix(NA, 2, 10)
rownames(sd.mat) = c("Tau", "SMR")
colnames(sd.mat) = c(1:10)
sd.mat[1,] = 1.96*sqrt(S.cox.tau.var[order.ind2])
sd.mat[2,] = 1.96*sqrt(apply(SMR.estimate.boot, 2, var)[order.ind2]/(real.SMR.365[order.ind2])^2) # by Delta Method

### Figure S1 (b) (Supporting Information)
pdf("Figure/tencenters_oneyear_highest_log.pdf", width = 10, height = 6)
par(mfrow = c(1,1), cex.lab = 1.5, cex.axis = 1.5, 
    mar=c(4,6,1,1), tcl = 0.5, oma = c(2, 2, 2, 3), xpd = FALSE, cex.main = 2)
ze_barplot <- barplot(mat , beside=T , 
                      ylim=c(-4,2) , ylab="Center effect estimates",
                      xlab = expression(paste("Centers with ten highest ", hat(tau)[j], "(t)")),
                      col = c("dodgerblue", "gold"))
#points(c(1:10), mat[2,], pch = 19)
legend(x = 17, y = 2, c(expression(paste(hat(tau)[j], "(t)")),
                        expression(paste("Estimated log(", hat(SMR)[j], ")"))),
       fill=c("dodgerblue" , "gold"), bty = "n", cex = 1.5)
abline(h = 0, lwd = 2)
error.bar(ze_barplot, mat, sd.mat)
dev.off()
