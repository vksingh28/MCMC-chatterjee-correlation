library(XICOR)

reps = 1e5

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

# MH

# MH independent

# Time irreversible MC
