GaussApprox <- function(f,
                        n,
                        Q,
                        k = 5,
                        h = 0.005,
                        x = NULL) {
    dsmall <- .Machine$double.eps^0.5
    m <- ncol(Q)
    if (is.null(x))
        x <- double(m)
    for (j in 1:k) {
        fx <- f(x[1:n])
        fa  <- f(x[1:n] - h)
        fb  <- f(x[1:n] + h)
        cc <- (2*fx -fa -fb)/(h^2)
        cc[cc<dsmall] <- dsmall
        bb <- (fb - fa)/(2*h) + x[1:n]*cc
        Qn <- Q + Diagonal(m, c(cc, rep(0, m-n)))
        L <- chol(Qn) ### "aparentlty" L is not use, but IS
        x <- drop(solve(Qn, c(bb, rep(0, m-n))))
    }
    return(list(x=x, bb=bb, cc=cc, Q=Qn, L=L))
}

n <- 300
x <- sin(3*pi*1:n/n) 

plot(x, type = 'l')

mu0 <- 0
y <- rpois(n, exp(mu0 + x))
par(mfrow=c(1,1), mar=c(2,3,0.5,0.5), mgp=c(2,0.5,0))
plot(y, pch=19)

## define the likeihood function
llp <- function(llam)
    dpois(y, exp(llam), log=TRUE)

## define the precision matrix structure
library(Matrix)
rworder <- 2
R <- crossprod(diff(Diagonal(n), differences = rworder))
attr(R, "rank") <- n-rworder

R[1:10, 1:10]

## the precision Q = tau * R
## try different tau
ga1 <- GaussApprox(llp, n, R*300)
ga2 <- GaussApprox(llp, n, R*3000)
ga3 <- GaussApprox(llp, n, R*30000)
gas <- list(ga1, ga2, ga3)

for(i in 1:length(gas))
    gas[[i]]$sd <- sqrt(diag(solve(gas[[i]]$Q)))

par(mfrow = c(2, 1), mar = c(3,3,0.5,0.5),
    mgp = c(1.5,0.5,0), bty = 'n')
plot(mu0 + x, type = 'l', lwd = 2,
     ylim = range(mu0+gas[[1]]$x -2*gas[[1]]$sd,
                  mu0+gas[[1]]$x +2*gas[[1]]$sd))
for(i in 1:3) {
    lines(gas[[i]]$x, col = i+1, lty = 2)
    lines(gas[[i]]$x - 2*gas[[i]]$sd, lty = 2, col = i+1)
    lines(gas[[i]]$x + 2*gas[[i]]$sd, lty = 2, col = i+1)
}
legend("topright", bty = "n",
       legend = c("TRUE",
                  paste("Posterior", 1:3, "and CI")), 
       lty = c(1,2:4), col = 1:4)
plot(y, pch = 19)
for(i in 1:3) {
    lines(exp(gas[[i]]$x), col = i+1, lty = 2)
    lines(exp(gas[[i]]$x - 2*gas[[i]]$sd), lty = 2, col = i+1)
    lines(exp(gas[[i]]$x + 2*gas[[i]]$sd), lty = 2, col = i+1)
}
legend("topright", bty = "n",
       legend = c("TRUE",
                  paste("Posterior", 1:3, "and CI")), 
       lty = c(1,2:4), col = 1:4)

## remains to improve it by doing either 
## 1. Laplace approximation for x
## 2. low rank VB correction
## and the inference for tau

library(INLA)

## Gaussian approximation with INLA
fga <- y ~ 0 +
    f(i, model='rw2', constr = FALSE,
      hyper = list(prec = list(param = c(1, 0.01))),
      vb.correct = FALSE)
fit.ga <- inla(fga, family = 'poisson', data = list(y = y, i = 1:n))

## GA+VBC: https://www.jmlr.org/papers/v25/21-1405.html
fvbc <- y ~ 0 +
    f(i, model='rw2', constr = FALSE,
      hyper = list(prec = list(param = c(1, 0.01))),
      vb.correct = TRUE)
fit.vbc <- inla(fvbc, family = 'poisson', data = list(y = y, i = 1:n))

### take the tau (fitted by INLA) and to the GA with that
ga.otau <- GaussApprox(llp, n, R * exp(fit.ga$mode$theta))
ga.otau$sd <- diag(solve(ga.otau$Q))

par(mfrow = c(1,1))
plot(mu0 + x, type = 'l', lwd = 2)
lines(ga.otau$x, col = 2, lty = 2)
lines(fit.ga$summary.random$i$mean, col = 3, lty = 2)
lines(fit.vbc$summary.random$i$mean, col = 4, lty = 2)
