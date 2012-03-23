## ODE interface to MK models.

initial.tip.mkn.ode <- function(cache) {
  k <- cache$info$k
  y <- matrix(rep(c(rep(0, k)), k + 1), k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,]) <- 1

  y <- matrix.to.list(y)

  y.i <- cache$states
  y.i[is.na(y.i)] <- k + 1

  tips <- cache$tips
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}

###########################################################################
## Additional functions
## For historical and debugging purposes, not used directly in the
## calculations, but branches function is generated this way
## internally.
make.branches.mkn <- function(cache, control)
  make.branches.dtlik(cache$info, control)
