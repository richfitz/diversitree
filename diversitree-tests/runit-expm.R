test.expm <- function() {
  library(diversitree)
  library(expm) # slow
  library(RUnit)

  set.seed(1)
  n <- 64L
  Q <- diversitree:::mkn.Q(runif(n*(n-1)))
  v <- runif(n)
  v <- v / sum(v)
  t <- 0.02

  ## 1. Plain matrix exponential
  E.f <- matrix(.Fortran("dexpmf", Q, n, t, out=numeric(n*n),
                         iflag=numeric(1))$out, n, n)
  E.e <- expm(Q*t)
  all.equal(E.f, E.e)

  ## 2. exp(Qt)v
  ans.f <- .Fortran("ddexpmv", Q, n, v, t, out=numeric(n),
                    iflag=numeric(1))$out
  ans.e <- c(expm(Q*t) %*% v)
  all.equal(ans.f, ans.e)

  ## 3. sparse exp(Qt)v
  set.seed(1)
  i <- sort(sample(n*(n-1), 2*n))
  p <- rep(0.0, n*(n-1))
  p[i] <- runif(2*n)
  Q2 <- diversitree:::mkn.Q(p)

  ans.e <- c(expm(Q2*t) %*% v)

  ans.d <- diversitree:::expmv.expokit.dense(Q2, t, v)
  all.equal(ans.e, ans.d)

  Q2.sparse <- diversitree:::expm.expokit.sparse.pars(Q2)
  ans.s <- c(diversitree:::expmv.expokit.sparse(Q2, t, v, tol=1e-10))
  all.equal(ans.s, ans.d)

  ## Then, with a vector of times:
  tt <- seq(0, t, length=11)[-1]

  ans.ee <- sapply(tt, function(t) c(expm(Q2*t) %*% v))
  ans.dd <- diversitree:::expmv.expokit.dense(Q2, tt, v)
  ans.ss <- diversitree:::expmv.expokit.sparse(Q2.sparse, tt, v, tol=1e-10)

  all.equal(ans.ee, ans.dd)
  all.equal(ans.ee, ans.ss)
}
