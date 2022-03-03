library(XICOR)
library(qbld)

reps = 1e5
k = 1e3
# AR(1)
ar.x = numeric(reps)
eps = rexp(reps, rate = .01);
ar.x[1] = 0
rho = 0.8
for(t in 2:reps) {
	ar.x[t] = rho*ar.x[t-1] + eps[t]
}

xicor(head(ar.x, reps-k), tail(ar.x, reps-k))
xicor(tail(ar.x, reps-k), head(ar.x, reps-k))
# Comments: xicor is symmetric. I believe the reason why its not coming out to be same is low value of rho and probably the use of correlated samples (not iid) to calculate xicor.

# MH
normals = rnorm(reps)
uniforms = runif(reps)

target <- function(x){
  return (dbeta(x, shape1 = 5, shape2 = 2))
}

# Symmetric MH
mh.x = numeric(reps)

mh.x[1] <- rexp(1, .01)
for(i in 2:reps){
  current_x <- mh.x[i-1]
  proposed_x <- current_x + normals[i]
  A <- min(1, target(proposed_x)/target(current_x))
  ifelse(uniforms[i] < A, mh.x[i] <- proposed_x, mh.x[i] <- current_x)
}

xicor(head(mh.x, reps-k), tail(mh.x, reps-k))
xicor(tail(mh.x, reps-k), head(mh.x, reps-k))
# Comment: it seems symmetric in this case

# Independence MH
mh_ind.x = numeric(reps)

mh_ind.x[1] <- rexp(1)
mean = 2
for(i in 2:reps){
  current_x <- mh_ind.x[i-1]
  proposed_x <- mean + normals[i]
  A <- min(1, (target(proposed_x)/target(current_x))*(dnorm(current_x, mean = mean, sd = 1)/dnorm(proposed_x, mean = mean, sd = 1)))
  ifelse(uniforms[i] < A, mh_ind.x[i] <- proposed_x, mh_ind.x[i] <- current_x)
}

xicor(head(mh_ind.x, reps-k), tail(mh_ind.x, reps-k))
xicor(tail(mh_ind.x, reps-k), head(mh_ind.x, reps-k))
# Comment: it seems symmetric in this case as well

#Gibbs sampler
gibbs = function (n, rho){    # a gibbs sampler implementation of a bivariate random number generator
    mat = matrix(ncol = 2, nrow = n)   # matrix for storing the random samples
    x = 0
    y = 0
    mat[1, ] = c(x, y)        # initialize the markov chain
    for (i in 2:n) {
            x = rnorm(1, rho * y, sqrt(1 - rho^2))        # sample from x conditional on y
            y = rnorm(1, rho * x, sqrt(1 - rho^2))        # sample from y conditional on x
            mat[i, ] = c(x, y)
    }
    mat
}
bvn = gibbs(reps,0.98)
calculateXI(head(bvn, reps-k), tail(bvn, reps-k))
calculateXI(tail(bvn, reps-k), head(bvn, reps-k))

# Time irreversible markov chain


# non stationary example to test symmetricity
n = 1e6
x1 = rnorm(n)
x2 = rnorm (n)
x3 = rnorm(n)
y = x1*x2 + sin(x1*x3)
calculateXI(x1, y)
calculateXI(y, x1)
# Not symmetric

?model.qbld
