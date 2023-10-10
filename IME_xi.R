library(XICOR)
library(httpgd)

reps = 1e4
k = 1e2

# AR(1) process
ar.x = numeric(reps)
sd = 100
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
# jpeg('acf_vs_chatterjee_ar_1.jpg')
# par(mfrow=c(2, 1))
acf(ar.x)
plot(0:40, acf_chatterjee[1:41], pch=16, cex=0.1, ylim=c(-.05, 1), xlab="Lag", ylab="ACF_Chatterjee", title="AR(1)")
segments(0:40, 0, 0:40, lwd=1, acf_chatterjee[1:51])
abline(h = .01, col = "blue", lty = 2)
abline(h = 0, col = "black")
abline(h = -.01, col = "blue", lty = 2)
# par(mfrow=c(1, 1))
# dev.off()
head(acf_chatterjee, 20)
acf(ar.x, plot=FALSE)

# Initial Positive Sequence Estimate
first_negative_chatterjee = min(which(acf_chatterjee < 0))-1
first_negative_acf = min(which(acf(ar.x)$acf < 0))-1
acf_values = as.vector(acf(ar.x)$acf)
var_pearson_positive = -acf_values[1] + 2*sum(acf_values[1:first_negative_acf])
var_chatterjee_positive = -acf_values[1] + 2*sum(acf_values[1:first_negative_chatterjee])
cat(var_pearson_positive, var_chatterjee_positive, "\n")

# Initial Monotone Sequence Estimate
for (i in 2:length(acf_chatterjee)) {
  if (acf_chatterjee[i] >= acf_chatterjee[i - 1]) {
    break  # Exit the loop when it stops decreasing
  }
  first_monotone_chatterjee <- i
}
for (i in 2:length(acf_values)) {
  if (acf_values[i] >= acf_values[i - 1]) {
    break  # Exit the loop when it stops decreasing
  }
  first_monotone_acf <- i
}

var_pearson_monotone = -acf_values[1] + 2*sum(acf_values[1:first_monotone_acf])
var_chatterjee_monotone = -acf_values[1] + 2*sum(acf_values[1:first_monotone_chatterjee])
cat(var_pearson_monotone, var_chatterjee_monotone, "\n")
