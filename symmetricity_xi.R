library(XICOR)

reps = 1e6

# AR(1)
ar.x = numeric(reps)
ar.y = numeric(reps-1)
eps = rexp(reps, rate = .01);
ar.x[1] = 0
rho = 0.25
for(t in 2:reps) {
	ar.x[t] = rho*ar.x[t-1] + eps[t]
	ar.y[t-1] = ar.x[t]
}
xicor(head(ar.x, reps-1), ar.y)
xicor(ar.y, head(ar.x, reps-1))
# Comments: Either xi is not symmetric, or there is some condition that is missing in the exp example

# MH
normals = rnorm(reps)
uniforms = runif(reps)

target <- function(x){
  return (dnorm(x, mean = 0, sd = 1))
}

# Symmetric MH
mh.x = numeric(reps)
mh.y = numeric(reps-1)

mh.x[1] <- rexp(1, .01)
for(i in 2:reps){
  current_x <- mh.x[i-1]
  proposed_x <- current_x + normals[i]
  A <- min(1, target(proposed_x)/target(current_x))
  ifelse(uniforms[i] < A, mh.x[i] <- proposed_x, mh.x[i] <- current_x)
}
for(i in 1:reps-1){
	mh.y[i] = mh.x[i+1]
}

xicor(head(mh.x, reps-1), mh.y)
xicor(mh.y, head(mh.x, reps-1))
# Comment: it seems symmetric in this case

# Independence MH
mh_ind.x = numeric(reps)
mh_ind.y = numeric(reps-1)

mh_ind.x[1] <- rexp(1)
mean = 2
for(i in 2:reps){
  current_x <- mh_ind.x[i-1]
  proposed_x <- mean + normals[i]
  A <- min(1, (target(proposed_x)/target(current_x))*(dnorm(current_x, mean = mean, sd = 1)/dnorm(proposed_x, mean = mean, sd = 1)))
  ifelse(uniforms[i] < A, mh_ind.x[i] <- proposed_x, mh_ind.x[i] <- current_x)
}
for(i in 1:reps-1){
	mh_ind.y[i] = mh_ind.x[i+1]
}
xicor(head(mh_ind.x, reps-1), mh_ind.y)
xicor(mh_ind.y, head(mh_ind.x, reps-1))
# Comment: it seems symmetric in this case as well


# Time irreversible markov chain
