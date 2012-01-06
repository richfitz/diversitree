## 1: make
make.quasse.split <- function(tree, states, states.sd, lambda, mu,
                              nodes, split.t, control=NULL,
                              sampling.f=NULL, defaults=NULL) {
  cache <- make.cache.quasse.split(tree, states, states.sd,
                                   lambda, mu, nodes, split.t,
                                   control, sampling.f)
  control <- cache$control # TODO: hmmm? Different to BiSSE

  if ( control$tips.combined )
    ## TODO: Will not be too hard.
    stop("tips.combined not yet implemented for split QuaSSE")

  ## TODO (1): Should use sanitised nodes?
  ## TODO (2): Remove this restriction.  This is an issue as the
  ## auxilliary variables will go in the wrong place...
  ## TODO (3): This might be moveable now.
  nodes.i <- cache$nodes - length(tree$tip.label)
  if ( min(branching.times(tree)[nodes.i]) < control$tc )
    stop("Sorry - all nodes must root below tc")

  initial.conditions.quasse <- make.initial.conditions.quasse(control) 

  if ( control$method == "fftC" ) {
    ## TODO: Should take arguments cache, control, and return the
    ## correct branches functions switching on control.
    branches.main <- make.branches.quasse.fftC(control)
    branches.aux <- make.branches.aux.quasse.fftC(control,
                                                  cache$sampling.f)
    branches <- make.branches.split(cache, branches.main, branches.aux)
    initial.conditions <-
      make.initial.conditions.split(cache, initial.conditions.quasse)
  } else {
    ## This is because the .aux functions are not done for these
    ## methods...
    stop("Alternative methods not yet implemented")
  }

  n.part <- cache$n.part
  args <- cache$args

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.f=NULL, intermediates=FALSE) {

    ## 1: Check parameters (move this to simplify...)
    names(pars) <- NULL # Because of use of do.call, strip names
    if ( length(pars) != cache$n.args )
      stop(sprintf("Incorrect number of arguments (expected %d, got %d)",
                   cache$n.args, length(pars)))
    ## Checking and extent construction:
    drift <- pars[unlist(args[,3])]
    diffusion <- pars[unlist(args[,4])]

    ext <- quasse.extent(cache$control, drift, diffusion)

    ## expand the parameters, with our current extent.
    pars.l <- lapply(seq_len(n.part), function(i)
                     expand.pars.quasse(cache$lambda[[i]],
                                        cache$mu[[i]],
                                        args[i,], ext, pars))

    if ( control$caching.branches )
      caching.branches.set.pars(pars.l, branches)

    ## *More* parameter checking:
    lambda.x <- unlist(lapply(pars.l, function(x)
                              c(x[[1]]$lambda, x[[2]]$lambda)))
    mu.x <- unlist(lapply(pars.l, function(x) c(x[[1]]$mu, x[[2]]$mu)))
    if ( any(lambda.x < 0) || any(mu.x < 0) || any(diffusion <= 0) )
      stop("Illegal negative parameters")
    if ( !any(lambda.x > 0) )
      stop("No positive lambda; cannot compute likelihood")

    cache$y <- initial.tip.quasse.split(cache, cache$control, ext$x[[1]])

    ans <- all.branches.list(pars.l, cache, initial.conditions, branches)
    pars.root <- pars.l[[cache$group.nodes[cache$root]]] # always 1
    vals <- matrix(ans$init[[cache$root]],
                   cache$control$nx, 2)[seq_len(ext$ndat[2]),] 
    root.p <- root.p.quasse(vals, ext, root, root.f)
    
    loglik <- root.quasse(vals, pars.root$lo$lambda, ans$lq, condition.surv,
                          root.p, cache$control$dx)

    if ( intermediates )
      attr(loglik, "intermediates") <- ans

    loglik  
  }

  if ( !is.null(defaults) )
    ll <- set.defaults(ll, defaults)

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
  cache <- make.cache.quasse(tree, states, states.sd, lambda[[1]], mu,
                             control, sampling.f, TRUE)
  cache <- make.cache.split(tree, cache, nodes, split.t)

  n.part <- cache$n.part
  control <- cache$control <- check.control.split(cache$control)

  cache$sampling.f <- check.sampling.f.split(sampling.f, 1, n.part)
  cache$aux.i <- seq_len(control$nx)
 
  ## 5: Speciation/extinction functions
  ## Most of the processing of these is currently done below.
  tmp <- check.f.quasse.split(lambda, n.part)
  cache$lambda <- tmp$f
  n.lambda <- tmp$n
  names.lambda <- tmp$names

  tmp <- check.f.quasse.split(mu, n.part)
  cache$mu <- tmp$f
  n.mu <- tmp$n
  names.mu <- tmp$names

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

## This is actually a bad time sink; can take 0.01s from a 0.024s
## minimal evaluation.  Really, we don't need to compute this each
## time, but just subset by excluding nkl and nkr points.
initial.tip.quasse.split <- function(cache, control, x) {
  nx <- control$nx * control$r
  npad <- nx - length(x)
  tips <- cache$tips
  e0 <- lapply(cache$sampling.f, function(f)
               rep(1 - f, nx))[cache$group.nodes[tips]]
  y <- mapply(function(e0, mean, sd)
              c(e0, dnorm(x, mean, sd), rep(0, npad)),
              e0, cache$states, cache$states.sd, SIMPLIFY=FALSE)
  
  dt.tips.ordered(y, tips, cache$len[tips])
}
