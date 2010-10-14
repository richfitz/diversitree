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

## 1: make
make.bisse <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, strict=TRUE, safe=FALSE) {
  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra,
                            strict=strict)
  branches <- make.branches.bisse(safe)
  ll <- function(pars, ...) ll.bisse(cache, pars, branches, ...)
  class(ll) <- c("bisse", "function")
  ll
}

## 2: print
print.bisse <- function(x, ...) {
  cat("BiSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.bisse <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
  else
    ret
}
`argnames<-.bisse` <- function(x, value) {
  if ( length(value) != 6 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.bisse <- function(func, x.init, method,
                           fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.bisse")
}

mcmc.bisse <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
make.cache.bisse <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10,
                             strict=TRUE) {
  ## TODO: There is a potential issue here with states, as
  ## 'unresolved' may contain one of the states.  For now I am
  ## disabling the check, but this is not great.
  if ( strict && !is.null(unresolved) ) {
    strict <- FALSE
  }
  
  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=strict, strict.vals=0:1)

  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- make.unresolved.bisse(tree$clades, states)
    states <- states[tree$tip.label]
  }

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, 2)

  ## Check 'unresolved' (there is certainly room to streamline this in
  ## the future).

  cache <- make.cache(tree)
  cache$tip.state  <- states
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f

  unresolved <- check.unresolved(cache, unresolved, nt.extra)
  cache$unresolved <- unresolved

  if ( !is.null(unresolved) ) {
    cache$tips <- cache$tips[-unresolved$i]
    cache$tip.state <- cache$tip.state[-unresolved$i]
  }

  ## This avoids a warning in the case where all tips are unresolved.
  if ( length(cache$tips) > 0 )
    cache$y <- initial.tip.bisse(cache)

  cache
}

## By the time this hits, unresolved clades and any other non-standard
## tips have been removed.  We have an index "tips" (equal to 1:n.tip
## for plain bisse) that is the "index" (in phy$edge numbering) of the
## tips, and a state vector cache$tip.state, both of the same length.
## The length of the terminal branches is cache$len[cache$tips].
##
## Allowing for unknown state tips, there are three possible states
##   (0, 1, NA -> 1, 2, 3)
## Initial conditions at the tips are given by their tip states:
## There are three types of initial condition in bisse:
##   state0: c(f_0, 0,   1-f_0, 1-f_1)
##   state1: c(0,   f_1, 1-f_0, 1-f_1)
##   state?: c(f_0, f_1, 1-f_0, 1-f_1)

initial.tip.bisse <- function(cache) {
  f <- cache$sampling.f
  y <- list(c(1-f, f[1], 0),
            c(1-f, 0, f[2]),
            c(1-f, f))
  y.i <- cache$tip.state + 1
  y.i[is.na(y.i)] <- 3

  tips <- cache$tips

  ## This would return a data structure appropriate for the more basic
  ## tip treatment:
  ##   dt.tips.ordered(y[y.i], tips, cache$len[tips])
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}


## 6: ll
ll.bisse <- function(cache, pars, branches, prior=NULL,
                     condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE,
                     root.p0=NA, root.p1=NA) {
  if ( !is.null(prior) )
    stop("'prior' argument to likelihood function no longer accepted")
  if ( length(pars) != 6 )
    stop("Invalid parameter length (expected 6)")
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)

  if ( !is.na(root.p0) ) {
    warning("root.p0 is deprecated: please use root.p instead")
    root.p <- c(root.p0, 1-root.p0)
  } else if ( !is.na(root.p1) ) {
    warning("root.p1 is deprecated: please use root.p instead")
    root.p <- c(1-root.p1, root.p1)
  }
  if ( !is.null(root.p) &&  root != ROOT.GIVEN )
    warning("Ignoring specified root state")

  if ( !is.null(cache$unresolved) )
    cache$preset <- 
      branches.unresolved.bisse(pars, cache$unresolved)

  ll.xxsse(pars, cache, initial.conditions.bisse, branches,
           condition.surv, root, root.p, intermediates)
}

## 7: initial.conditions:
initial.conditions.bisse <- function(init, pars, t, is.root=FALSE)
  c(init[[1]][c(1,2)],
    init[[1]][c(3,4)] * init[[2]][c(3,4)] * pars[c(1,2)])

## 8: branches
make.branches.bisse <- function(safe=FALSE) {
  RTOL <- ATOL <- 1e-8
  bisse.ode <- make.ode("derivs_bisse", "diversitree",
                        "initmod_bisse", 4, safe)
  branches <- function(y, len, pars, t0)
    t(bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  
  make.branches(branches, 3:4)
}

## 9: branches.unresolved
branches.unresolved.bisse <- function(pars, unresolved) {
  Nc <- unresolved$Nc
  k <- unresolved$k
  nsc <- unresolved$nsc
  t <- unresolved$len
  nt <- max(Nc) + unresolved$nt.extra
  
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  base <- bucexpl(nt, mu0, mu1, lambda0, lambda1, q01, q10, t,
                  Nc, nsc, k)[,c(3,4,1,2),drop=FALSE]

  q <- rowSums(base[,3:4,drop=FALSE])
  base[,3:4] <- base[,3:4] / q

  list(target=unresolved$target,
       lq=log(q),
       base=matrix.to.list(base))
}

## Additional functions
bisse.stationary.freq <- function(pars) {
  .Deprecated("stationary.freq.bisse")
  stationary.freq.bisse(pars)
}
stationary.freq.bisse <- function(pars) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  g <- (lambda0 - mu0) - (lambda1 - mu1)
  eps <- (lambda0 + mu0 + lambda1 + mu1)*1e-14
  if ( abs(g) < eps ) {
    if ( q01 + q10 == 0 )
      0.5
    else
      q10/(q01 + q10)
  } else {
    roots <- quadratic.roots(g, q10+q01-g, -q10)
    roots <- roots[roots >= 0 & roots <= 1]
    if ( length(roots) > 1 )
      NA
    else
      roots
  }
}

starting.point.bisse <- function(tree, q.div=5, yule=FALSE) {
  ## TODO: Use qs estimated from Mk2?  Can be slow is the only reason
  ## I have not set this up by default.
  ## find.mle(constrain(make.mk2(phy, phy$tip.state), q10 ~ q01), .1)$par
  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  if  ( pars.bd[1] > pars.bd[2] )
    p <- rep(c(pars.bd, (pars.bd[1] - pars.bd[2]) / q.div), each=2)
  else
    p <- rep(c(pars.bd, pars.bd[1] / q.div), each=2)
  names(p) <- argnames.bisse(NULL)
  p
}
starting.point <- bisse.starting.point <- function(tree, q.div=5) {
  .Deprecated("starting.point.bisse")
  starting.point.bisse(tree, q.div)
}

## This is here for reference, but not exported yet.  It should be
## tweaked in several ways
##   1. Starting parameter guessing should be done internally, at
##      least as an option.
##   2. Better listing of arguments
##   3. Automatic parsing of results into some sort of table; this
##      proabably requires classing this.
all.models.bisse <- function(f, p, ...) {
  f3 <- constrain(f, lambda1 ~ lambda0, mu1 ~ mu0, q01 ~ q10)
  f4.lm <- constrain(f, lambda1 ~ lambda0, mu1 ~ mu0)
  f4.lq <- constrain(f, lambda1 ~ lambda0, q01 ~ q10)
  f4.mq <- constrain(f, mu1 ~ mu0, q01 ~ q10)
  f5.l <- constrain(f, lambda1 ~ lambda0)
  f5.m <- constrain(f, mu1 ~ mu0)
  f5.q <- constrain(f, q01 ~ q10)

  ## Fit six and three parameter models
  if ( length(p) != 3 )
    stop("Starting point must be of length 3 (lambda, mu, q)")
  ans3 <- find.mle(f3, p, ...)

  ## Using the values from the 3p model, fit the 4p and 5p models:
  l <- ans3$par[1]
  m <- ans3$par[2]
  q <- ans3$par[3]

  ## Start the searches from the best model of the previous type.
  ans4.lm <- find.mle(f4.lm, c(l, m, q, q), ...)
  ans4.lq <- find.mle(f4.lq, c(l, m, m, q), ...)
  ans4.mq <- find.mle(f4.mq, c(l, l, m, q), ...)

  p.l <- if ( ans4.lm$lnLik > ans4.lq$lnLik )
    ans4.lm$par[c(1:2,2:4)] else ans4.lq$par[c(1:4,4)]
  p.m <- if ( ans4.lm$lnLik > ans4.mq$lnLik )
    ans4.lm$par[c(1,1:4)] else ans4.mq$par[c(1:4,4)]
  p.q <- if ( ans4.lq$lnLik > ans4.mq$lnLik )
    ans4.lq$par[c(1,1:4)] else ans4.mq$par[c(1:3,3:4)]
  ans5.l  <- find.mle(f5.l, p.l, ...)
  ans5.m  <- find.mle(f5.m, p.m, ...)
  ans5.q  <- find.mle(f5.q, p.q, ...)

  tmp <- list(ans5.l, ans5.m, ans5.q)
  i <- which.max(sapply(tmp, "[[", "lnLik"))
  p6 <- tmp[[i]]$par
  j <- list(c(1, 1:5), c(1:2, 2:5), c(1:5, 5))
  ans6 <- find.mle(f, p6[j[[i]]], ...)

  list(ans6=ans6,
       ans5.l =ans5.l,  ans5.m =ans5.m,  ans5.q =ans5.q,
       ans4.lm=ans4.lm, ans4.lq=ans4.lq, ans4.mq=ans4.mq,
       ans3=ans3)
}

check.unresolved <- function(cache, unresolved, nt.extra) {
  if ( is.null(unresolved) ) {
    return(NULL)
  } else if ( nrow(unresolved) == 0 ) {
    warning("Ignoring empty 'unresolved' argument")
    return(NULL)
  }

  required <- c("tip.label", "Nc", "n0", "n1")
  if ( !all(required %in% names(unresolved)) )
    stop("Required columns missing from unresolved clades")

  unresolved$tip.label <- as.character(unresolved$tip.label)
  if ( !all(unresolved$tip.label %in% cache$tip.label) )
    stop("Unknown tip species in 'unresolved'")

  if ( max(unresolved$Nc + nt.extra) > 200 )
    stop("The largest unresolved clade supported has %d species",
         200 - nt.extra)

  unresolved$i <- match(unresolved$tip.label, cache$tip.label)
  unresolved$target <- cache$tips[unresolved$i]

  unresolved$k   <- unresolved$n1
  unresolved$nsc <- unresolved$n0 + unresolved$n1
  unresolved <- as.list(unresolved)
  unresolved$nt.extra <- nt.extra

  unresolved$len <- cache$len[unresolved$target]

  unresolved
}
