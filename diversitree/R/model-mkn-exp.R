## TODO: This is somewhat misleading, as I don't actually use matrix
## exponentiation; I'm computing Pij = exp(Qt) through integrating a
## series of k^2 ODEs.
make.all.branches.mkn.exp <- function(cache, control) {
  k <- cache$info$k
  
  f.pij <- make.pij.mkn(cache$info, control)

  idx.tip <- cache$idx.tip
  n.tip <- cache$n.tip
  n <- length(cache$len)

  ## TODO: deal with unknown states; tip calculations for unknown tip
  ## states just return the sum of nonzero ys from the matrix
  ## multiplication, so they can just be done separately.
  map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
  idx.tip <- cbind(c(map[cache$states,]),
                   rep(seq_len(n.tip), k))
  len.uniq <- sort(unique(cache$len))
  len.idx <- match(cache$len, len.uniq)

  ## Alter things to make it more speedy.  The '.C' denotes C-style
  ## base-0 indices.
  children.C <- toC.int(t(cache$children))
  order.C    <- toC.int(cache$order)
  
  ## At this point, the parameters are assumed to be a Q matrix
  function(pars, intermediates, preset=NULL) {
    if ( !is.null(preset) )
      stop("Preset values not allowed")
    pij <- f.pij(len.uniq, pars)[,len.idx]

    lq <- numeric(n)
    branch.init <- branch.base <- matrix(NA, k, n)
    storage.mode(branch.init) <- "numeric"

    ## tips
    ans <- matrix(pij[idx.tip], n.tip, k)
    q <- rowSums(ans)
    branch.base[,seq_len(n.tip)] <- t.default(ans/q)
    lq[seq_len(n.tip)] <- log(q)

    ans <- .C("r_mkn_core",
              k        = as.integer(k),
              n        = length(order.C) - 1L,
              order    = order.C,
              children = children.C,
              pij      = pij,
              init     = branch.init,
              base     = branch.base,
              lq       = lq,
              NAOK=TRUE, DUP=FALSE)
    list(init=ans$init,
         base=ans$base,
         lq=ans$lq,
         vals=ans$init[,cache$root],
         pij=pij)
  }
}

######################################################################
## Mkn-special stuff.
## The calculations here are quite different to the rest of the
## package, so thre is a lot below.

## Compute Pij matrices:
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
make.pij.mkn <- function(info, control) {
  if ( control$method == "mk2" )
    return(pij.mk2)
  control <- check.control.ode(control)
  ## TODO/NEW: cvodes should be easy too.
  if ( control$backend != "deSolve" )
    stop("Only deSolve backend available")

  k <- info$k
  info <- list(name="mkn_pij",
               ny=k*k, np=k*k, idx.d=integer(0))
  info <- check.info.ode(info, control)
  pij.ode <- make.ode.deSolve(info, control)

  yi <- diag(k) # initial conditions always same.

  function(len, pars)
    pij.ode(yi, c(0, len), pars)
}
