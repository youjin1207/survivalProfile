set.seed(1234)
randoms = sample(1:100, 100)
time.list = seq(1, 30, 1)

prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1)
X = apply(rmultinom(500, 1, prob = prob), 2, function(x) which(x == 1)) / 10

### Figure 1 (a) (main manuscript)
pdf("Figure/diagram.pdf", width = 16, height = 8)
par(mfrow = c(1,2), cex.lab = 2.5, cex.axis = 2, 
    mar=c(4,8,3,1), tcl = 0.5, oma = c(2, 2, 3, 3), xpd = FALSE, cex.main = 2)
plot(seq(0.1, 1, 0.1), as.numeric(table(X))/500, bty = "n",
     xlab = "X", main = "", xaxt="n", ylab = "Density", yaxt = "n",
     xlim = c(0, 1.1), ylim = c(0, 0.3), cex = 1, pch  = 19, col = "royalblue")
probs = as.numeric(table(X))/500
for(q in 1:10){
  segments(x0 = 0.1*q, y0 = 0, x1 = 0.1*q, y1 = probs[q], col = "royalblue", lwd = 2)
}
axis(1, at = seq(0.1, 1, 0.1), label = seq(0.1, 1, 0.1))
axis(2, at = seq(0, 0.3, 0.1), label = seq(0, 0.3, 0.1))

dummy = runif(10, 0.5, 0.9)
plot(seq(0.1, 1, 0.1), dummy, bty = "n",
     xlab = "X", main = "", yaxt="n", ylab = expression(paste("Pr(", T^j, " > t | ", X, ")")), 
     ylim = c(0.4, 1), xlim = c(0, 1.1), cex = 1.5, pch = 19,  col = "orangered",
     xaxt = "n")
for(q in 1:10){
  segments(x0 = 0.1*q, y0 = 0.4, x1 = 0.1*q, y1 = dummy[q], col = "orangered", lwd = 2)
}
axis(1, at = seq(0.1, 1, 0.1), label = seq(0.1, 1, 0.1))
axis(2, at = seq(0.4, 1, 0.2), label = seq(0.4, 1, 0.2))

dev.off()


prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1) # scenario (ii)
prob = rep(prob, each = 50)
X = apply(rmultinom(500*50, 1, prob = prob), 2, function(x) which(x == 1)) / 500
dummy = c(runif(100, 0.5, 0.8), runif(100, 0.7, 1.0), runif(100, 0.6, 0.9),
          runif(100, 0.5, 0.7), runif(100, 0.4, 0.7))
cuts = c(as.numeric(quantile(X, seq(0.2, 1, 0.2))))

### Figure 1 (b) (main manuscript)
pdf("Figure/diagram_score.pdf", width = 16, height = 8)
par(mfrow = c(1,2), cex.lab = 2.5, cex.axis = 2, 
    mar=c(4,8,3,1), tcl = 0.5, oma = c(2, 2, 3, 3), xpd = FALSE, cex.main = 2)
plot(seq(1/500, 1, 1/500), as.numeric(table(X))/(500*50), bty = "n",
     xlab = expression(paste(eta[i], "(X)")), main = "", xaxt="n", ylab = "Density", yaxt = "n",
     xlim = c(0, 1.1), ylim = c(0, 0.005), cex = 0.5, pch  = 19, col = "royalblue")
probs = as.numeric(table(X))/(500*50)
for(q in 1:500){
  segments(x0 = (1/500)*q, y0 = 0, x1 = ((1/500))*q, y1 = probs[q], col = "royalblue", lwd = 0.3)
}
axis(1, at = seq(0.1, 1, 0.1), label = seq(0.1, 1, 0.1))
axis(2, at = seq(0, 0.005, 0.001), label = seq(0, 0.005, 0.001))
abline(v = cuts, lwd = 2, col = "black")
text(x = cuts[1]/2, y = 0.005, "Q=1", cex = 1.5)
text(x = 0.36, y = 0.005, "Q=2", cex = 1.5)
text(x = 0.5, y = 0.005, "Q=3", cex = 1.5)
text(x = 0.62, y = 0.005, "Q=4", cex = 1.5)
text(x = 0.85, y = 0.005, "Q=5", cex = 1.5)

plot(seq(1/500, 1, 1/500), dummy, bty = "n",
     xlab = expression(paste(eta[i], "(X)")), main = "", yaxt="n", ylab = expression(paste("Pr(", T^j, " > t | ", X, ")")), 
     ylim = c(0.4, 1), xlim = c(0, 1.1), cex = 0.5, pch = 19,  col = "orangered",
     xaxt = "n")
for(q in 1:500){
  segments(x0 = (1/500)*q, y0 = 0.4, x1 = (1/500)*q, y1 = dummy[q], col = "orangered", lwd = 0.3)
}
axis(1, at = seq(0.1, 1, 0.1), label = seq(0.1, 1, 0.1))
axis(2, at = seq(0.4, 1, 0.2), label = seq(0.4, 1, 0.2))
abline(v = cuts, lwd = 2, col = "black")
text(x = cuts[1]/2, y = 1, "Q=1", cex = 1.5)
text(x = 0.36, y = 1, "Q=2", cex = 1.5)
text(x = 0.5, y = 1, "Q=3", cex = 1.5)
text(x = 0.62, y = 1, "Q=4", cex = 1.5)
text(x = 0.85, y = 1, "Q=5", cex = 1.5)

dev.off()