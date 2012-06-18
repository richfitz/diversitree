library(diversitree)
dyn.load("diversitree/src/expm.so")
use.expm <- FALSE

set.seed(1)
n <- 64L
Q <- diversitree:::mkn.Q(runif(n*(n-1)))
v <- runif(n)
v <- v / sum(v)
t <- 0.02

## 1. Plain matrix exponential
E.f <- matrix(.Fortran("dfexpmv", Q, n, t, out=numeric(n*n),
                       iflag=numeric(1))$out, n, n)
if ( use.expm ) {
  require(expm)
  E.e <- expm(Q*t)
  all.equal(E.f, E.e)
}

## 2. exp(Qt)v
ans.f <- .Fortran("ddexpmv", Q, n, v, t, out=numeric(n),
                  iflag=numeric(1))$out
if ( use.expm ) {
  require(expm)
  ans.e <- c(expm(Q*t) %*% v)
  all.equal(ans.f, ans.e)
}

## 3. sparse exp(Qt)v

set.seed(1)
i <- sort(sample(n*(n-1), 2*n))
p <- rep(0.0, n*(n-1))
p[i] <- runif(2*n)
Q2 <- diversitree:::mkn.Q(p)

ans.f <- .Fortran("ddexpmv", Q2, n, v, t, out=numeric(n),
                  iflag=numeric(1))$out
if ( use.expm ) {
  require(expm)
  ans.e <- c(expm(Q2*t) %*% v)
  all.equal(ans.e, ans.f)
}

## Pre-processing
idx <- which(Q2 != 0, TRUE)
nz <- as.integer(nrow(idx))
qnorm <- max(abs(Q2[i]))
tol <- 1e-8

## This takes two runs for some reason?
ans.s <- .Fortran("dsexpmv", Q2[idx], n, idx[,1], idx[,2], nz,
                  qnorm, v, t, tol, out=numeric(n), iflag=integer(1))$out

if ( use.expm )
  all.equal(ans.e, tmp$out)

tt <- seq(0, t, length=11)[-1]
lt <- as.integer(length(tt))

ans.t <- .Fortran("dsexpmvi", Q2[idx], n, idx[,1], idx[,2], nz,
                  qnorm, v, tt, lt, tol, out=numeric(n*lt),
                  iflag=integer(1))$out
ans.t <- matrix(ans.t, n, lt)

target <- sapply(tt, function(t) c(expm(Q2*t) %*% v))
all.equal(ans.t, target)

f.pars <- function(Q) {
  idx <- which(Q2 != 0, TRUE)
  nz <- as.integer(nrow(idx))
  qq <- Q2[i]
  list(Q=qq, ia=idx[,1], ja=idx[,2], nz=nrow(idx), qnorm=max(abs(qq)))
}

