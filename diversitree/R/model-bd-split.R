## Split BD model.  This is basically MEDUSA, but a more concrete
## likelihood-function based version, rather than the highly optimised
## search version that Alfaro et al. describe.

## In the future, I hope to use the update() function to efficiently
## push this to match up with MEDUSA's capabilities.  As it stands,
## this will not be an efficient way of looping over nodes and running
## optimisations in the same way as MEDUSA.  But it should provide a
## decent reference implementation of the calculations.

## 1: make
make.bd.split <- function(tree, nodes, split.t, sampling.f=NULL,
                          unresolved=NULL) {
  if ( missing(split.t) )
    split.t <- rep(Inf, length(nodes))
  cache <- make.cache.bd.split(tree, nodes, split.t, sampling.f,
                               unresolved)
  n.part <- cache$n.part
  ll.part <- lapply(seq_len(n.part), make.bd.split.part, cache=cache)

  ll <- function(pars, condition.surv=TRUE, ...) {
    if ( length(pars) != 2 * n.part )
      stop(sprintf("Expected %d parameters, but got %d",
                   2 * n.part, length(pars)))
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)
    
    pars <- matrix(pars, 2, n.part)

    res <- numeric(n.part)
    for ( i in seq_len(n.part) )
      res[i] <- ll.part[[i]](pars[,i], condition.surv)
    sum(res)
  }

  class(ll) <- c("bd.split", "function")
  attr(ll, "n.part") <- cache$n.part
  
  ll
}

## 2: print
print.bd.split <- function(x, ...) {
  cat("BD(split) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.bd.split <- function(x, ...) {
  obj <- attr(x, "argnames")
  n <- attr(x, "n.part")
  if ( is.null(obj) )
    obj <- list(base=c("lambda", "mu"),
                levels=seq_len(n))

  paste(obj$base, rep(obj$levels, each=2), sep=".")
}
`argnames<-.bd.split` <- function(x, value) {
  n <- attr(x, "n")
  if ( !is.list(value) || length(value) != 2 )
    stop("'value' must be a list of length 2")
  if ( length(value[[1]]) != 2 || length(value[[2]]) != n )
    stop(sprintf("value's elements must be of length 2, %d", n))

  names(value) <- c("base", "levels")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.bd.split <- function(func, x.init, method, fail.value=NA,
                              ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method,
             class.append="fit.mle.bd.split")
}

## 5: make.cache
make.cache.bd.split <- function(tree, nodes, split.t=Inf,
                                sampling.f=NULL, unresolved=NULL) {
  tree <- check.tree(tree)

  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")

  if ( !isTRUE(all.equal(unique(split.t), Inf, check.attr=FALSE)) )
    stop("split.t cannot yet be changed")
  nodes <- check.split(tree, nodes, rep(Inf, length(nodes)))$nodes
  n <- length(nodes)

  sampling.f <- check.sampling.f(sampling.f, n)

  ## Unresolved here is a different format to that expected by bisse.
  if ( is.null(unresolved) ) {
    n.taxa <- 1
  } else {
    ## TODO: This might not be correct.
    n.taxa <- unresolved$n.taxa[match(tree$tip.label, unresolved$tip.label)]
    n.taxa[is.na(n.taxa)] <- 1
  }

  edge <- tree$edge
  n.tip <- length(tree$tip.label)
  bt <- as.numeric(branching.times(tree))

  i <- seq_len(max(edge))
  j <- match(i, edge[,2])
  z <- cbind(anc=edge[j,1], dec=i,
             t.0=NA,
             t.1=bt[match(edge[j,1], (n.tip+1):max(edge))],
             t.len=tree$edge.length[j],
             n0=1, nt=NA,
             group=NA)
  z[,"t.0"] <- z[,"t.1"] - z[,"t.len"]

  z[match(seq_len(n.tip), z[,2]),"nt"] <- n.taxa
  z[n.tip + 1,1] <- n.tip + 1

  group <- make.split.phylo.vec(tree, nodes)[j]
  group[n.tip + 1] <- 1
  z[,"group"] <- group

  obj <- list(z=z, n.taxa=n.tip, n.node=tree$Nnode,
              sampling.f=sampling.f, t.root=max(bt),
              n.part=n)
}

make.bd.split.part <- function(cache, i) {
  z <- cache$z[cache$z[,"group"] == i,]
  f <- cache$sampling.f[i]
  n.node <- sum(is.na(z[,"nt"]))    

  ## Determine if this is the root group, and if so remove the root
  ## (needs to be done after the node counting, as the root node still
  ## contributes a lambda term).
  is.root <- is.na(z[,"t.1"])
  if ( any(is.root) ) {
    z <- z[!is.root,]
    is.root <- TRUE
    t.root <- cache$t.root

    ng <- length(cache$sampling.f)
    n <- tabulate(cache$z[!is.na(cache$z[,"nt"]),"group"], ng)

    ## Where all the sampling.f are the same, this will be
    ## lfactorial(n.taxa) + n.taxa * log(f)
    root.constant <- 
      lfactorial(cache$n.taxa - 1) + sum(n*log(cache$sampling.f))
  } else {
    is.root <- FALSE
  }

  t0 <- z[,"t.0"]
  t1 <- z[,"t.1"]
  dt <- z[,"t.len"]
  
  function(pars, condition.surv=TRUE) {
    lambda <- pars[1]
    mu <- pars[2]
    r <- lambda - mu

    ## The abs() here is justified by this being
    ##   log(x^2) -> 2 log(abs(x))
    d <- r * dt +
      2*(log(abs((f * exp(r * t0) + 1-f) * lambda - mu)) -
         log(abs((f * exp(r * t1) + 1-f) * lambda - mu)))

    log.lik <- sum(d) + n.node * log(lambda)

    if ( is.root ) {
      log.lik <- log.lik + root.constant

      if ( condition.surv )
        log.lik <- log.lik -
          log(f * f * r * (1 - mu/lambda)) +
            2*log(abs(exp(-r * t.root)*(mu/lambda - 1 + f) - f))
    }

    log.lik
  }
}

