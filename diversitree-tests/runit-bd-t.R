test.bd.t2 <- function() {
  library(diversitree)
  library(RUnit)
  
  set.seed(1)
  ## BUG: tree.bd chokes with wrong length parameters...
  phy <- tree.bd(c(.1, 0), max.taxa=30)

  p0 <- c(.1, 0)
  p1 <- c(.1, .05)
  p0.t <- c(p0[1], 0, p0[2])
  p1.t <- c(p1[1], 0, p1[2])
  p2.t <- c(p1[1], .005, p1[2])

  lik1 <- make.bd(phy, control=list(method="ode"))
  ll0 <- lik1(p0)
  ll1 <- lik1(p1)

  lik2 <- make.bd.t.old(phy, constant.t)
  checkEquals(lik2(p0), ll0)
  checkEquals(lik2(p1), ll1)

  ## old version
  lik3 <- make.bd.t.old(phy, list(linear.t, constant.t))
  checkEquals(lik3(p0.t), ll0)
  checkEquals(lik3(p1.t), ll1)
  ll2 <- -20.332980072
  checkEquals(lik3(p2.t), ll2)

  ## new version
  lik4 <- make.bd.t(phy, c("linear.t", "constant.t"))
  checkEquals(lik4(p0.t), ll0)
  checkEquals(lik4(p1.t), ll1)
  checkEquals(lik4(p2.t), ll2)

  ## Also with CVODES:
  lik5 <- make.bd.t(phy, c("linear.t", "constant.t"),
                    control=list(backend="CVODES"))
  checkEquals(lik5(p0.t), ll0, tol=2e-7)
  checkEquals(lik5(p1.t), ll1, tol=2e-7)
  checkEquals(lik5(p2.t), ll2, tol=2e-7)

  ## Now, a spline fit:
  t.max <- max(branching.times(phy)) * 1.001

  ## First, the full function (a sin wave down to the root of the
  ## tree).
  sin.t <- function(t, y0, y1)
    y0 + (y1 - y0) * (sin(t/t.max*2*pi)+1)/2
  lik6 <- make.bd.t.old(phy, list(sin.t, constant.t))

  p3.t <- c(p0[1], p0[1], p0[2])
  p4.t <- c(p0[1], 2*p0[1], p0[2])

  checkEquals(lik6(p3.t), ll0)
  ll4 <- -23.0593106966584
  checkEquals(lik6(p4.t), ll4)

  ## Then, with a simple fit through this.
  x <- seq(0, t.max, length.out=101)
  y <- sin(x/t.max*2*pi)
  spline.data <- list(t=x, y=y)

  lik7 <- make.bd.t(phy, c("spline.t", "constant.t"),
                    spline.data=spline.data)
  checkEquals(lik7(p3.t), ll0)
  checkEquals(lik7(p4.t), ll4)
}
