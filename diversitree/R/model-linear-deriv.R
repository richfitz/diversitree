## TODO: Keep an R version handy?
## TODO: mkn or linear-deriv?  Or have mkn.deriv be a special case of
##   linear.deriv? (yes)
make.mkn.deriv <- function(tree, states, k, strict=FALSE) {
  cache <- make.cache.mkn.deriv(tree, states, k, strict)

  nd <- as.integer(k)
  np <- as.integer(k*(k-1))
  fail <- structure(-Inf, gr=rep(NA, np))

  t <- cache$len
  state0 <- cache$state0
  order0 <- cache$order0
  children0 <- cache$children0

  Q <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  Q.d <- unlist(cache$Q.d)

  ll.mkn.deriv <- function(pars, with.gr=TRUE) {
    if ( length(pars) != np )
      stop(sprintf("Invalid length parameters (expected %d)", np))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(fail)

    Q[idx] <- pars
    diag(Q) <- -rowSums(Q)

    es <- eigen(Q, FALSE)
    d <- es$values
    Amat <- es$vectors
    Ainv <- solve(Amat)

    .Call("linear_deriv", nd, np,
          d, Amat, Ainv, Q.d, t,
          state0, order0, children0, as.logical(with.gr),
          package="diversitree")
  }

  class(ll.mkn.deriv) <- c("mkn.deriv", "mkn", "function")
  attr(ll.mkn.deriv, "k") <- k
  ll.mkn.deriv
}

## 2: print
print.mkn.deriv <- function(x, ...) {
  cat("Mk-n likelihood function [with gradients]:\n")
  print(unclass(x))
}

## 3. argnames / argnames<- (inherited from mkn)
## 4. find.mle (inherited from mkn, for now)

## 5: make.cache
make.cache.mkn.deriv <- function(tree, states, k, strict) {
  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=strict,
                         strict.vals=1:k)
  if ( any(is.na(states)) )
    stop("Unknown states not yet allowed")

  cache <- make.cache(tree) 
  cache$k <- k
  cache$np <- np <- k * (k - 1)
  cache$tip.state  <- states

  cache$Q.d <- make.deriv.Q.mkn(k)

  cache$state0 <- as.integer(cache$tip.state - 1)
  cache$children0 <- as.integer(t(cache$children-1))
  cache$order0 <- as.integer(cache$order-1)

  cache
}

## 6: ll
## make.ll.linear.deriv <- function(cache, nd, np) {
##   nd <- as.integer(cache$nd)
##   np <- as.integer(cache$np)
##   fail <- structure(-Inf, gr=rep(NA, np))
  
##   t <- cache$len
##   state0 <- cache$state0
##   order0 <- cache$order0
##   children0 <- cache$children0

##   f.Q   <- cache$f.Q
##   f.Q.d <- cache$f.Q.d

##   ll.linear.deriv <- function(pars, with.gr=TRUE) {
##     if ( length(pars) != np )
##       stop(sprintf("Invalid length parameters (expected %d)", np))
##     if ( any(!is.finite(pars)) || any(pars < 0) )
##       return(fail)

##     Q <- f.Q(pars)
##     Q.d <- f.Q.d(pars)

##     es <- eigen(Q, FALSE)
##     d <- es$values
##     Amat <- es$vectors
##     Ainv <- solve(Amat)

##     .Call("linear_deriv", nd, np,
##           d, Amat, Ainv, Q.d, t,
##           state0, order0, children0, as.logical(with.gr),
##           package="diversitree")
##   }
## }


## Additional functions:
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

