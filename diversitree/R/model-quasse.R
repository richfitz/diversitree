
## 1: make
make.quasse <- function(tree, states, states.sd, lambda, mu,
                        control=NULL, sampling.f=NULL) {
  cache <- make.cache.quasse(tree, states, states.sd, lambda, mu,
                             control, sampling.f)

  ## TODO: branches here changes to allow backend switching, basd on
  ## cache$control$method...
  control <- cache$control
  tips <- NULL
  
  if ( control$method == "fftC" ) {
    branches <- make.branches.quasse.fftC(control)
    if ( control$tips.combined )
      tips <- make.tips.quasse.fftC(control, cache$len[cache$tips],
                                   cache$tips)
  } else if ( control$method == "fftR" ) {
    branches <- make.branches.quasse.fftR(control)
  }

  initial.conditions <- make.initial.conditions.quasse(control)

  ll <- function(pars, ...)
    ll.quasse(cache, pars, branches, initial.conditions, tips, ...)
  attr(ll, "f.lambda") <- lambda
  attr(ll, "f.mu") <- mu
  class(ll) <- c("quasse", "function")
  ll
}

## 2: print
print.quasse <- function(x, ...) {
  cat("QuaSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.quasse <- function(x, ...) {
  f.lambda <- attr(x, "f.lambda")
  f.mu <- attr(x, "f.mu")
  if ( is.null(f.lambda) || is.null(f.mu) )
    stop("Corrupt QuaSSE function")
  c(sprintf("l.%s", names(formals(f.lambda))[-1]),
    sprintf("m.%s", names(formals(f.mu))[-1]),
    "drift", "diffusion")
}
`argnames<-.quasse` <- function(x, value) {
  ## It might be easier to get at cache?
  ## TODO: fix this if the cache structure changes.
  if ( length(value) != environment(x)$cache$n.args )
    stop("Invalid names length")
  else
    attr(x, "argnames") <- value
  x
}

## 4: find.mle
find.mle.quasse <- function(func, x.init, method, fail.value=NA,
                            verbose=TRUE, ...) {
  if (missing(method)) 
    method <- "subplex"
  find.mle.default(func, x.init, method,
                   fail.value, "fit.mle.quasse", verbose=verbose, ...)
}

## 5: make.cache:
make.cache.quasse <- function(tree, states, states.sd, lambda, mu,
                              control, sampling.f) {
  ## 1: tree
  tree <- check.tree(tree)  

  ## 2: states & errors
  tmp <- check.states.quasse(tree, states, states.sd)
  states <- tmp$states
  states.sd <- tmp$states.sd

  ## 3: Control structure (lots of checking!)
  control <- check.control.quasse(control, tree, states)

  ## 4: Speciation/extinction functions
  n.lambda <- check.f.quasse(lambda)
  n.mu     <- check.f.quasse(mu)
  n.args   <- n.lambda + n.mu + 2
  args <- list(lambda=seq_len(n.lambda),
               mu=seq_len(n.mu) + n.lambda,
               drift=n.lambda + n.mu + 1,
               diffusion=n.lambda + n.mu + 2)

  sampling.f <- check.sampling.f(sampling.f, 1)

  ## The raw cache just provides the traversal information:
  cache <- make.cache(tree)
  cache$states  <- states
  cache$states.sd <- states.sd
  cache$sampling.f <- sampling.f
  cache$lambda <- lambda
  cache$mu <- mu

  cache$control <- control
  
  cache$args <- args
  cache$n.args <- n.args

  cache
}

initial.tip.quasse <- function(cache, control, x) {
  nx <- control$nx * control$r
  npad <- nx - length(x)
  e0 <- 1 - cache$sampling.f

  if ( control$tips.combined ) {
    tips <- cache$tips
    t <- cache$len[tips]
    i <- order(t)
    target <- tips[i]

    states <- cache$states[i]
    states.sd <- cache$states.sd[i]

    y <- mapply(function(mean, sd)
                c(dnorm(x, mean, sd), rep(0, npad)),
                states, states.sd, SIMPLIFY=FALSE)
    y <- matrix(c(rep(e0, nx), unlist(y)), nx, length(target)+1)

    list(target=target, y=y, t=t[i])
  } else {
    y <- mapply(function(mean, sd)
                c(rep(e0, nx), dnorm(x, mean, sd), rep(0, npad)),
                cache$states, cache$states.sd, SIMPLIFY=FALSE)
    dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
  }
}

## 6: ll
ll.quasse <- function(cache, pars, branches, initial.conditions, tips,
                      condition.surv=TRUE, root=ROOT.OBS,
                      root.p=NULL, intermediates=FALSE) {
  names(pars) <- NULL # Because of use of do.call, strip names
  args <- cache$args

  drift <- pars[args$drift]
  diffusion <- pars[args$diffusion]

  ext <- quasse.extent(cache$control, drift, diffusion)
  ## This would confirm the translation:
  ##   all.equal(ext$x[[1]][ext$tr], ext$x[[2]])

  ## Parameters, expanded onto the extent:
  pars <- expand.pars.quasse(cache$lambda, cache$mu, args, ext, pars)

  lambda.x <- c(pars$hi$lambda, pars$lo$lambda)
  mu.x <- c(pars$hi$mu, pars$lo$mu)

  if ( any(lambda.x < 0) || any(mu.x < 0) || diffusion <= 0 )
    stop("Illegal negative parameters")
  if ( !any(lambda.x > 0) )
    stop("No positive lambda; cannot compute likelihood")

  init <- initial.tip.quasse(cache, cache$control, ext$x[[1]])

  if ( cache$control$tips.combined ) {
    cache$preset <- tips(init$y, pars)

    ## supress tip calculation now:
    cache$tips <- integer(0) 
    cache$y <- NULL
  } else {
    cache$y <- init
  }

  ans <- all.branches(pars, cache, initial.conditions, branches)
  
  vals <- matrix(ans$init[[cache$root]], cache$control$nx, 2)

  ## TODO: this assumes that the root node is in the low-condition.  I
  ## think that this is enforced by the checking, but I cannot
  ## remember.
  loglik <- root.quasse(vals[seq_len(ext$ndat[2]),], ans$lq,
                        cache$control$dx, pars$lo$lambda,
                        condition.surv)
  if ( intermediates )
    attr(loglik, "intermediates") <- ans

  loglik
}

## 7: initial.conditions
make.initial.conditions.quasse <- function(control) {
  tc <- control$tc
  r <- control$r
  nx.lo <- control$nx
  nx.hi <- nx.lo * r

  ## There is the chance that we could be slightly off on the depth
  ## by rounding error.  Because of this, I've done the testing
  ## against the *length* of the data, and then checked that the time
  ## is appropriate (to within eps of the correct value).  It is
  ## possible that two different branches with different numbers of
  ## nodes that finish right at the critical interval might have
  ## conflicting lengths.
  eps <- 1e-8
  function(init, pars, t, is.root=FALSE) {
    if ( length(init[[1]]) != length(init[[2]]) )
      stop("Data have incompatible length")

    if ( t < tc ) {
      ## if ( length(init[[1]]) / 2 == nx.hi ) { # t < tc
      ## if ( !((t - eps) < tc) )
      ##   stop("Wrong data size")
      nx <- nx.hi
      lambda <- pars[[1]]$lambda
    } else {
      ## if ( !((t + eps) > tc) )
      ##   stop("Wrong data size")
      nx <- nx.lo
      lambda <- pars[[2]]$lambda
    }
    
    ndat <- length(lambda)
    i <- seq_len(nx)
    j <- seq.int(nx+1, nx + ndat)

    c(init[[1]][i],
      init[[1]][j] * init[[2]][j] * lambda,
      rep.int(0.0, nx - ndat))
  }
}

## 8: Branches - see backend files for the bits that make up f.hi and
## f.lo
make.branches.quasse <- function(f.hi, f.lo, control) {
  nx <- control$nx
  dx <- control$dx
  tc <- control$tc
  r <- control$r
  
  function(y, len, pars, t0) {
    if ( t0 >= tc ) {
      ans <- f.lo(y, len, pars$lo, t0)
    } else if ( t0 + len < tc ) {
      ans <- f.hi(y, len, pars$hi, t0)
      dx <- dx / r
    } else {
      len.hi <- tc - t0
      ans.hi <- f.hi(y, len.hi, pars$hi, t0)
      y.lo <- ans.hi[pars$tr,]
      if ( nrow(y.lo) < nx )
        y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr), 2))
      ans <- f.lo(y.lo, len - len.hi, pars$lo, t0)
    }

    ## This should probably be togglable.
    ans[,2][ans[,2] < 0] <- 0

    q <- sum(ans[,2]) * dx
    ans[,2] <- ans[,2] / q
    c(log(q), ans)
  }
}


## TODO: A bunch more options to come in here for different root
## styles, etc.  Then drop the 'old.style' bit (hence the inefficiency
## there).
## Still to do: ROOT.FLAT, ROOT.GIVEN (ROOT.GIVEN requires computing
## things against the x extent...)

## Normally, the probability of speciation at the root and subsequent
## survival is computed as
##   lambda * (1 - E)^2
## Here I will compute
##   \int p(x) \lambda(x) (1 - E(x))^2 dx
root.quasse <- function(vars, log.comp, dx, lambda, condition.surv) {
  e.root <- vars[,1]
  d.root <- vars[,2]

  p.root <- d.root / (sum(d.root) * dx)

  if ( condition.surv ) {
    ## d.root <- d.root / (lambda * (1-e.root)^2)
    p.surv <- sum(p.root * lambda * (1 - e.root)^2) * dx
    d.root <- d.root / p.surv
  }

  log(sum(p.root * d.root) * dx) + sum(log.comp)
}

