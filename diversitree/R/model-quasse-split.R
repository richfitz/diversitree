## 1: make
make.quasse.split <- function(tree, states, states.sd, lambda, mu,
                              nodes, split.t, control=NULL,
                              sampling.f=NULL) {
  cache <- make.cache.quasse.split(tree, states, states.sd,
                                   lambda, mu, nodes, split.t,
                                   control, sampling.f)
  control <- cache$control

  ## TODO: This is an annoying limit, but requires shuffling around
  ## aux.i in an unfortunate way.
  if ( min(branching.times(tree)[nodes]) < control$tc )
    stop("Sorry - all nodes must root below tc")
  
  if ( control$method == "fftC" ) {
    branches <- make.branches.quasse.fftC(control)
    branches.aux <- make.branches.aux.quasse.fftC(control,
                                                  cache$sampling.f)
    tips <- vector("list", length(cache$cache))
    if ( control$tips.combined ) {
      for ( i in seq_along(tips) ) {
        x <- cache$cache[[i]]
        tips[[i]] <- make.tips.quasse.fftC(control, x$len[x$tips],
                                           x$tips)
      }
    }
  } else {
    ## I am not sure what the technical reason for this not being
    ## ready is, but it can't be that bad.
    stop("Alternative methods not yet implemented")
  }
    
  initial.conditions <- make.initial.conditions.quasse(control)

  ll <- function(pars, ...)
    ll.quasse.split(cache, pars, branches, branches.aux,
                    initial.conditions, tips, ...)
  
  attr(ll, "n.part") <- cache$n.part
  attr(ll, "f.lambda") <- lambda
  attr(ll, "f.mu") <- mu
  class(ll) <- c("quasse.split", "quasse", "function")
  ll
}

## 2: print
print.quasse.split <- function(x, ...) {
  cat("QuaSSE/Split likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
## Argnames for QuaSSE/split is complicated compared with BiSSE/split,
## because there are infinite possible functions that may be used for
## the different functions.  Also, each different component may have
## different functions, so the basic "twopart" idea does not work here.
argnames.quasse.split <- function(x, ...) {
  ## TODO: update to use the argsnames system
  args.names <- environment(x)$cache$args.names
  n <- length(args.names)
  if ( n != attr(x, "n.part") )
    stop("Corrupt QuaSSE/split function")
  m <- sapply(args.names, length)
  sprintf("%s.%d", unlist(args.names), rep(seq_len(n), m))
}
`argnames<-.quasse.split` <- function(x, value) {
  .NotYetImplemented()
}



## 4: find.mle (inherited from quasse), though this will be changed
## soon to better deal with cases where single parameter moves are
## done.

## 5: make.cache:
make.cache.quasse.split <- function(tree, states, states.sd,
                                    lambda, mu, nodes, split.t,
                                    control, sampling.f) {
  ## 1: tree
  tree <- check.tree(tree, node.labels=TRUE)

  ## 2: states & errors
  tmp <- check.states.quasse(tree, states, states.sd)
  states <- tmp$states
  states.sd <- tmp$states.sd

  ## 3: Control structure
  control <- check.control.quasse(control, tree, states)

  n <- length(nodes) + 1 # +1 for base group

  ## 4: Check sampling.f
  sampling.f <- as.list(check.sampling.f(sampling.f, n))

  ## 5: Speciation/extinction functions
  ## Most of the processing of these is currently done below.
  tmp <- check.f.quasse.split(lambda, n)
  lambda <- tmp$f
  n.lambda <- tmp$n
  names.lambda <- tmp$names

  tmp <- check.f.quasse.split(mu, n)
  mu <- tmp$f
  n.mu <- tmp$n
  names.mu <- tmp$names
  
  ## 6: Generic cache splitting:
  cache <- make.cache.split(tree, nodes, split.t)

  for ( i in seq_len(n) ) {
    x <- cache$cache[[i]]
    x$states  <- states[x$tip.label]
    x$states.sd <- states.sd[x$tip.label]
    x$sampling.f <- sampling.f[[i]]
    x$lambda <- lambda[[i]]
    x$mu <- mu[[i]]

    cache$cache[[i]] <- x    
  }

  ## TODO: This causes a problem here: the auxillary information will
  ## not go in the right place through aux.i...
  if ( any(split.t < control$tc) )
    stop("split.t < control$tc not yet handled")
  cache$sampling.f <- sampling.f
  cache$aux.i <- seq_len(control$nx)
  cache$control <- control

  ## 6: Arguments and corresponding parameter names:
  i <- rbind(lambda=n.lambda, mu=n.mu, drift=1, diffusion=1)
  j <- split(seq_len(sum(i)), rep(seq_along(i), i))
  dim(j) <- dim(i)
  rownames(j) <- c("lambda", "mu", "drift", "diffusion")
  
  cache$args <- t(j)
  cache$n.args <- sum(i)

  cache$args.names <-
    mapply(c, lapply(names.lambda, sprintf, fmt="l.%s"),
           lapply(names.mu, sprintf, fmt="m.%s"),
           "drift", "diffusion", SIMPLIFY=FALSE)

  cache
}

## 6: ll
ll.quasse.split <- function(cache, pars, branches, branches.aux,
                            initial.conditions, tips,
                            condition.surv=TRUE, root=ROOT.OBS,
                            root.p=NULL, intermediates=FALSE) {
  n.part <- cache$n.part
  args <- cache$args

  names(pars) <- NULL # Because of use of do.call, strip names
  if ( length(pars) != cache$n.args )
    stop(sprintf("Incorrect number of arguments (expected %d, got %d)",
                 cache$n.args, length(pars)))

  ## Checking and extent construction:
  drift <- pars[unlist(args[,3])]
  diffusion <- pars[unlist(args[,4])]
  if ( any(drift != 0) )
    stop("Non-zero drift not yet implemented.") # TODO
  if ( any(diffusion <= 0) )
    stop("Strictly positive diffusion required")
  ext <- quasse.extent(cache$control, 0, max(diffusion))

  pars.l <- lapply(seq_len(n.part), function(i)
                   expand.pars.quasse(cache$cache[[i]]$lambda,
                                      cache$cache[[i]]$mu,
                                      args[i,], ext, pars))

  lambda.x <- unlist(lapply(pars.l, function(x)
                            c(x[[1]]$lambda, x[[2]]$lambda)))
  mu.x <- unlist(lapply(pars.l, function(x) c(x[[1]]$mu, x[[2]]$mu)))

  if ( any(lambda.x < 0) || any(mu.x < 0) )
    stop("Illegal negative parameters")
  if ( !any(lambda.x > 0) )
    stop("No positive lambda; cannot compute likelihood")

  init <- lapply(cache$cache, initial.tip.quasse,
                 cache$control, ext$x[[1]])

  if ( cache$control$tips.combined ) {
    for ( i in seq_len(n.part) ) {
      ## TODO: This is where the tip calculations would go, following
      ## along with the pattern in plain QuaSSE....
      cache$cache[[i]]$preset <- tips[[i]](init[[i]]$y, pars.l[[i]])
      ## supress tip calculation now:
      cache$cache[[i]]$tips <- integer(0) 
      cache$cache[[i]]$y <- NULL
      ## cache$cache[[i]]$y <- init[[i]]
    }
  } else {
    for ( i in seq_len(n.part) )
      cache$cache[[i]]$y <- init[[i]]
  }

  ## Normal:
  ans <- all.branches.split(pars.l, cache, initial.conditions,
                           branches, branches.aux)

  vars <- matrix(ans[[1]]$base, cache$control$nx, 2)
  lq <- unlist(lapply(ans, "[[", "lq"))

  ## TODO: this assumes that the root node is in the low-condition.  I
  ## think that this is enforced by the checking, but I cannot
  ## remember.
  loglik <- root.quasse(vars[seq_len(ext$ndat[2]),], lq,
                       cache$control$dx, pars.l[[1]]$lo$lambda,
                       condition.surv)

  if ( intermediates )
    attr(loglik, "intermediates") <- ans
  loglik
}

## 7: initial.conditions: from quasse

## 8: branches: from quasse.  However the 'branches.aux' function is
## required to compute the E0, E1 values after a partition.

## Other:
check.f.quasse.split <- function(f, rep) {
  if ( is.function(f) ) {
    n <- rep.int(check.f.quasse(f), rep)
    names <- rep.int(list(names(formals(f))[-1]), rep)
    f <- rep.int(list(f), rep)
  } else {
    if ( length(f) != rep )
      stop("Invalid length for speciation/extinction function")
    n <- sapply(f, check.f.quasse)
    names <- lapply(f, function(x) names(formals(x))[-1])
  }
  list(n=n, f=f, names=names)
}
