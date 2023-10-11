library(XICOR)
library(httpgd)


variance_estimates = function(rho){
  reps = 500

  # AR(1) process
  ar.x = numeric(reps)
  sd = 5
  eps = rnorm(reps, mean = 0, sd = sd)
  ar.x[1] = rnorm(1)
  for(t in 2:reps){
    ar.x[t] = rho*ar.x[t-1] + eps[t]
  }
  N = 150
  acf_chatterjee = numeric(N)
  for(i in 0:N){
    acf_chatterjee[i+1] = xicor(head(ar.x, reps-i), tail(ar.x, reps-i))
  }
  # jpeg('acf_vs_chatterjee_ar_1.jpg')
  # par(mfrow=c(2, 1))
  # acf(ar.x)
  # plot(0:40, acf_chatterjee[1:41], pch=16, cex=0.1, ylim=c(-.05, 1), xlab="Lag", ylab="ACF_Chatterjee", title="AR(1)")
  # segments(0:40, 0, 0:40, lwd=1, acf_chatterjee[1:51])
  # abline(h = .01, col = "blue", lty = 2)
  # abline(h = 0, col = "black")
  # abline(h = -.01, col = "blue", lty = 2)
  # # par(mfrow=c(1, 1))
  # # dev.off()
  # head(acf_chatterjee, 20)
  # acf(ar.x, plot=FALSE)
  # Initial Positive Sequence Estimate
  acf_values = as.vector(acf(ar.x, lag.max = N, plot = FALSE)$acf)
  # acf_values = as.vector(acf(ar.x)$acf)
  first_negative_chatterjee = min(which(acf_chatterjee < 0), N+1)-1
  first_negative_acf = min(which(acf_values < 0), N+1)-1
  var_pearson_positive = -acf_values[1] + 2*sum(acf_values[1:min(first_negative_acf, N)])
  var_chatterjee_positive = -acf_values[1] + 2*sum(acf_values[1:min(first_negative_chatterjee, N)])
  return(list(var_pearson = var_pearson_positive, var_chatterjee = var_chatterjee_positive, first_negative_acf = first_negative_acf, first_negative_chatterjee = first_negative_chatterjee))
  # cat(first_negative_acf, first_negative_chatterjee, "\n")
  # cat(var_pearson_positive, var_chatterjee_positive, "\n")

  # # Initial Monotone Sequence Estimate
  # for (i in 2:length(acf_chatterjee)) {
  #   if (acf_chatterjee[i] >= acf_chatterjee[i - 1]) {
  #     break  # Exit the loop when it stops decreasing
  #   }
  #   first_monotone_chatterjee <- i
  # }
  # for (i in 2:length(acf_values)) {
  #   if (acf_values[i] >= acf_values[i - 1]) {
  #     break  # Exit the loop when it stops decreasing
  #   }
  #   first_monotone_acf <- i
  # }

  # var_pearson_monotone = -acf_values[1] + 2*sum(acf_values[1:first_monotone_acf])
  # var_chatterjee_monotone = -acf_values[1] + 2*sum(acf_values[1:first_monotone_chatterjee])
  # cat(first_monotone_acf, first_monotone_chatterjee, "\n")
  # cat(var_pearson_monotone, var_chatterjee_monotone, "\n")
}
rho_values = seq(0.9, 0.99, length.out = 10)

# chatterjee_variances = data.frame(matrix(0, ncol = 10, nrow = 500))
# names(chatterjee_variances) = rho_values
# pearson_variances = data.frame(matrix(0, ncol = 10, nrow = 500))
# names(pearson_variances) = rho_values
# chatterjee_first_index = data.frame(matrix(0, ncol = 10, nrow = 500))
# names(chatterjee_first_index) = rho_values
# acf_first_index = data.frame(matrix(0, ncol = 10, nrow = 500))
# names(acf_first_index) = rho_values

variances = data.frame(matrix(0, ncol = 20, nrow = 500))
first_index = data.frame(matrix(0, ncol = 20, nrow = 500))

for(j in 1:length(rho_values)){
  for(i in 1:250){
    rho = rho_values[j]
    variance_values = variance_estimates(rho)
    variances[i, 2*j-1] = variance_values$var_chatterjee
    variances[i, 2*j] = variance_values$var_pearson
    first_index[i, 2*j-1] = variance_values$first_negative_chatterjee
    first_index[i, 2*j] = variance_values$first_negative_acf
    # chatterjee_variances[i, j] = variance_values$var_chatterjee
    # pearson_variances[i, j] = variance_values$var_pearson
    # chatterjee_first_index[i, j] = variance_values$first_negative_chatterjee
    # acf_first_index[i, j] = variance_values$first_negative_chatterjee

  }
}
boxplot_title = c()
for(rho in rho_values){
  boxplot_title = c(boxplot_title, paste(as.character(rho), "xi"), "acf")
}
names(variances) = boxplot_title
jpeg("acf_vs_chatterjee_boxplot.jpeg")
boxplot(variances[1:250, ], boxwidth = 10)
dev.off()
boxplot(first_index[1:250, ])
