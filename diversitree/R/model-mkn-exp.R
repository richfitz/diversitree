## TODO: This is somewhat misleading, as I don't actually use matrix
## exponentiation; I'm computing Pij = exp(Qt) through integrating a
## series of k^2 ODEs.

## Models should provide:
##   1. make
##   2. print
##   3. argnames / argnames<-
##   4. find.mle
## Generally, make will require:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches
make.mkn.exp <- function(tree, states, k, strict=TRUE, control=list()) {
  ## (eps is ignored for this model)
  control <- check.control.ode(control)
  use.mk2 <- control$use.mk2
  cache <- make.cache.mkn.exp(tree, states, k, use.mk2, strict)

  cache$f.pij <- if ( k == 2 && use.mk2 )
    pij.mk2 else make.pij.mkn(k, control$safe, control$tol)

  f.pars <- make.mkn.pars(k)

  ll.mkn <- function(pars, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE) {
    if (!is.null(root.p) && root != ROOT.GIVEN) 
      warning("Ignoring specified root state")
    check.pars.mkn(pars, k)
    
    ans <- all.branches.mkn.exp(f.pars(pars), cache)

    d.root <- ans$init[,cache$root]
    root.p <- root.p.mkn(d.root, pars, root, root.p)
    loglik <- root.mkn(d.root, ans$lq, root.p)

    if ( intermediates ) {
      ans$init[,cache$tips] <- cache$y$y[cache$y$i,]
      ans$root.p <- root.p
      attr(loglik, "intermediates") <- ans
      attr(loglik, "vals") <- d.root
    }
    loglik
  }

  class(ll.mkn) <- c("mkn.exp", "mkn", "function")
  attr(ll.mkn, "k") <- k
  ll.mkn
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.mkn.exp <- function(tree, states, k, use.mk2, strict) {
  tree <- check.tree(tree)
  if ( !is.null(states) ) # for multitrait
    states <- check.states(tree, states, strict=strict, strict.vals=1:k)

  cache <- make.cache(tree)
  cache$k <- k
  cache$tip.state  <- states
  ## TODO: deal with unknown states; tip calculations for unknown tip
  ## states just return the sum of nonzero ys from the matrix
  ## multiplication, so they can just be done separately.
  cache$map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
  cache$idx.tip <- cbind(c(cache$map[cache$tip.state,]),
                         rep(seq_len(cache$n.tip), k))
  cache$len.uniq <- sort(unique(cache$len))
  cache$len.idx <- match(cache$len, cache$len.uniq)
  if ( !is.null(states) ) # multitrait
    cache$y <- initial.tip.mkn.exp(cache)

  ## Alter things to make it more speedy.  The '0' denotes base zero
  ## indices.
  cache$children0 <- as.integer(t(cache$children-1))
  cache$order0 <- as.integer(cache$order-1)
  cache$use.mk2 <- use.mk2
  
  cache
}
initial.tip.mkn.exp <- function(cache) {
  k <- cache$k
  tip.state <- cache$tip.state
  if ( any(tip.state < 1 | tip.state > k, na.rm=TRUE) )
    stop(sprintf("tip states must be in the range [1, %d]", k))
  
  y <- matrix(rep(rep(0, k), k + 1), k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,1:k]) <- 1
  i <- cache$tip.state
  i[is.na(i)] <- k + 1
  list(y=y, i=i, types=sort(unique(i)))
}

## 6: ll (done within make.mkn)

## 7: initial.conditions (done in underlying C code)

## 8: branches (separate for mk2 and mkn)
pij.mk2 <- function(len, pars) {
  ## The 3,2 indices here are because pars is a Q matrix by the time
  ## this gets called.
  q01 <- pars[3]
  q10 <- pars[2]
  x <- exp(-(q01+q10)*len)
  rbind((x*q01 + q10),
        (1 - x)*q10,
        (1 - x)*q01,
        (x*q10 + q01)) / (q01 + q10)
}

make.pij.mkn <- function(k, safe=FALSE, tol=1e-8) {
  pij.ode <- make.ode("derivs_mkn_pij", "diversitree",
                      "initmod_mkn", k*k, safe)
  ATOL <- RTOL <- tol
  yi <- diag(k)
  
  function(len, pars)
    pij.ode(yi, c(0, len), pars,
            rtol=RTOL, atol=ATOL)[-1,-1,drop=FALSE]
}

## The new integrator, especially for Mkn.
all.branches.mkn.exp <- function(pars, cache) {
  ## At this point, the parameters are assumed to be a Q matrix
  k <- cache$k

  idx.tip <- cache$idx.tip
  n.tip <- cache$n.tip
  n <- length(cache$len)
  order0 <- cache$order0

  len.uniq <- cache$len.uniq
  len.idx <- cache$len.idx
  pij <- cache$f.pij(len.uniq, pars)[,len.idx]
  
  lq <- numeric(n)
  branch.init <- branch.base <- matrix(NA, k, n)
  storage.mode(branch.init) <- "numeric"

  ## tips
  ans <- matrix(pij[idx.tip], n.tip, k)
  q <- rowSums(ans)
  branch.base[,seq_len(n.tip)] <- t.default(ans/q)
  lq[seq_len(n.tip)] <- log(q)

  ## branches
  ans <- .C("r_mkn_core",
            k        = as.integer(k),
            n        = length(order0) - 1L,
            order    = order0,
            children = cache$children0,
            pij      = pij,
            init     = branch.init,
            base     = branch.base,
            lq       = lq,
            NAOK=TRUE, DUP=FALSE)

  list(init=ans$init,
       base=ans$base,
       lq=ans$lq,
       pij=pij)
}
