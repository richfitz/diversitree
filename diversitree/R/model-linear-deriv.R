## This uses a matrix exponentiation approach, rather than the usual
## ODE approach.  As such, the entire Pij matrix needs to be
## computed.  An ODE version might be faster, but this is just for
## playing really.
make.mkn.deriv <- function(tree, states, k, strict=FALSE) {
  cache <- make.cache.mkn.deriv(tree, states, k, strict)

  nd <- as.integer(k)
  np <- as.integer(k*(k-1))

  t         <- cache$len
  state.C    <- cache$state.C
  order.C    <- cache$order.C
  children.C <- cache$children.C
  Q.d       <- unlist(cache$Q.d)

  f.pars <- make.pars.mkn(k)

  ll <- function(pars, with.gr=TRUE) {
    qmat <- f.pars(pars)
    es <- eigen(qmat, FALSE)
    d <- es$values
    Amat <- es$vectors
    Ainv <- solve(Amat)

    browser()

    .Call("linear_deriv", nd, np,
          d, Amat, Ainv, Q.d, t,
          state.C, order.C, children.C, as.logical(with.gr),
          PACKAGE="diversitree")
  }

  class(ll) <- c("mkn.deriv", "mkn", "dtlik", "function")
  ll
}

make.info.mkn.deriv <- function(k, phy) {
  info <- make.info.mkn(k, phy)
  info$name <- "mkn.deriv"
  info$name.pretty <- "Mk(n) (with derivatives)"
  info$ny <- NA # special
  info$reference <- NULL # not sure
  info
}

make.cache.mkn.deriv <- function(tree, states, k, strict) {
  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=strict,
                         strict.vals=1:k)
  if ( any(is.na(states)) )
    stop("Unknown states not yet allowed")

  cache <- make.cache(tree) 
  cache$np <- np <- k * (k - 1)
  cache$states  <- states

  cache$Q.d <- make.deriv.Q.mkn(k)

  cache$state.C    <- toC.int(  cache$states)
  cache$children.C <- toC.int(t(cache$children))
  cache$order.C    <- toC.int(  cache$order)

  cache
}

make.deriv.Q.mkn <- function(k) {
  out <- rep(list(matrix(0.0, k, k)), k * (k-1))
  idx <- 1
  for ( i in 1:k ) {
    for ( j in 1:k ) {
      if ( i != j ) {
        out[[idx]][i,j] <- 1
        out[[idx]][i,i] <- -1
        idx <- idx + 1
      }
    }
  }
  out
}
