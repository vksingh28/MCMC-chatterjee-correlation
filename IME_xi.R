library(XICOR)
library(httpgd)

reps = 1e4
k = 1e2

# AR(1) process
ar.x = numeric(reps)
sd = 1
eps = rnorm(reps, mean = 0, sd = sd)
ar.x[1] = rnorm(1)
rho = 0.9
for(t in 2:reps){
	ar.x[t] = rho*ar.x[t-1] + eps[t]
}
N = 1e3
acf_chatterjee = numeric(N)
for(i in 0:1e3){
	acf_chatterjee[i+1] = xicor(head(ar.x, reps-i), tail(ar.x, reps-i))
}
jpeg('acf_vs_chatterjee_ar_1.jpg')
par(mfrow=c(2, 1))
acf(ar.x)
plot(0:50, acf_chatterjee[1:51], pch=16, cex=0.1, ylim=c(-.05, 1), xlab="Lag", ylab="ACF_Chatterjee", title="AR(1)")
segments(0:50, 0, 0:50, lwd=1, acf_chatterjee[1:51])
abline(h = .01, col = "blue", lty = 2)
abline(h = 0, col = "black")
abline(h = -.01, col = "blue", lty = 2)
par(mfrow=c(1, 1))
dev.off()
head(acf_chatterjee, 20)
acf(ar.x, plot=FALSE)
