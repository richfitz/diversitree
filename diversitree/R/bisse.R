
ROOT.FLAT  <- 1
ROOT.EQUI  <- 2
ROOT.OBS   <- 3
ROOT.GIVEN <- 4
ROOT.BOTH  <- 5

## TODO: Pass sampling.f through with the bissecache
## (probably?) or through with bisse.ll?  bissecache seems to make
## more sense, but it continues to cause problems with the wrong
## levels of abstraction.
make.bisse <- function(tree, states, unresolved=NULL,
                       sampling.f=NULL, ...) {
  cache <- bissecache(tree, states, unresolved, sampling.f, ...)

  f <- function(pars, ..., fail.value=NULL) {
    if ( !is.null(fail.value) )
      protect(bisse.ll, fail.value)(cache, pars, ...)
    else
      bisse.ll(cache, pars, ...)
  }

  class(f) <- c(class(f), "bisse")
  f
}

bisse.ll <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     condition.surv=TRUE, intermediates=FALSE,
                     root.p0=NA, root.p1=NA) {
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)
  if ( !missing(root.p0) || !missing(root.p1) )
    if ( missing(root) ) root <- ROOT.GIVEN
    else warning("Ignoring specified root state")
  if ( !missing(root.p1) )
    root.p0 <- 1 - root.p1
  ans <- bisse.branches(pars, cache)
  e.root <- ans$init[cache$root,1:2]
  d.root <- ans$init[cache$root,3:4]

  ## Just the stationary frequency for now
  if ( root == ROOT.FLAT )
    p <- 0.5
  else if ( root == ROOT.EQUI )
    p <- stationary.freq(pars)
  else if ( root == ROOT.OBS )
    p <- d.root[1]/sum(d.root)
  else if ( root == ROOT.GIVEN )
    p <- root.p0
  else if ( root != ROOT.BOTH )
    stop("Invalid root mode")

  ## Condition on survival
  if ( condition.surv )
    d.root <- d.root / (1-e.root)^2
  logcomp <- sum(ans$lq)
  if ( root == ROOT.BOTH )
    loglik <- log(c(sum(c(1,0) * d.root),
                    sum(c(0,1) * d.root)))
  else
    loglik <- log(sum(c(p, 1-p) * d.root)) + logcomp

  if ( !is.null(prior) )
    loglik <- loglik + bisse.prior(pars, prior)
  
  if ( intermediates ) {
    attr(loglik, "cache") <- cache
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- ans$init[cache$root,]
    attr(loglik, "logComp") <- logcomp
  }

  loglik
}

## bissecache does some caching that can be done for *every* set of
## parameters.  Mostly, I am caching tree topology information.
bissecache <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, ...) {
  if ( !is.null(unresolved) && nrow(unresolved) == 0 ) {
    unresolved <- NULL
    warning("Ignoring empty 'unresolved' argument")
  }

  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- make.unresolved(tree$clades, states)
    states <- states[tree$tip.label]
    names(states) <- tree$tip.label
  }

  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else if ( is.null(sampling.f) )
    sampling.f <- c(1, 1)
  else if ( length(sampling.f) != 2 )
    stop("sampling.f must be of length 2 (or NULL)")
  else if ( max(sampling.f) > 1 || min(sampling.f) < 0 )
    stop("sampling.f must be on range [0,1]")

  if ( !is.null(unresolved) ) {
    required <- c("tip.label", "Nc", "n0", "n1")
    if ( !all(required %in% names(unresolved)) )
      stop("Required columns missing from unresolved clades")
    if ( !all(unresolved$tip.label %in% tree$tip.label) )
      stop("Unknown tip species in 'unresolved'")
    unresolved$k   <- unresolved$n1
    unresolved$nsc <- unresolved$n0 + unresolved$n1
    unresolved$i   <- match(unresolved$tip.label, tree$tip.label)
  }

  if ( is.null(names(states)) ) {
    stop("The states vector must contain names")
  } else {
    known <- names(states)
    if ( !is.null(unresolved) )
      known <- unique(c(known, as.character(unresolved$tip.label)))
    if ( !all(tree$tip.label %in% known) )
      stop("Not all species have state information")
  }

  edge <- tree$edge
  idx <- seq_len(max(edge))
  ntip <- length(tree$tip.label)
  root <- ntip + 1
  
  is.tip <- idx <= ntip

  tip.state <- rep(NA, ntip)
  tip.state[] <- states[match(tree$tip.label, names(states))]

  children <- lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
  if ( !all(sapply(children, length)==2) )
    stop("Multifircations/unbranched nodes in tree - best get rid of them")
  children <- rbind(matrix(NA, ntip, 2), t(matrix(unlist(children), 2)))

  ans <- list(len=tree$edge.len[match(idx, edge[,2])],
              is.tip=is.tip,
              tip.state=tip.state,
              tip.label=tree$tip.label,
              children=children,
              order=get.ordering(children, is.tip, root),
              root=root,
              unresolved=unresolved,
              sampling.f=sampling.f,
              nt.extra=nt.extra)
  class(ans) <- "bissecache"
  ans
}

## Utility, used by bissecache
get.ordering <- function(children, is.tip, root) {
  todo <- list(root)
  i <- root
  repeat {
    kids <- children[i,]
    i <- kids[!is.tip[kids]]
    if ( length(i) > 0 )
      todo <- c(todo, list(i))
    else
      break
  }
  as.vector(unlist(rev(todo)))
}

## Within tips():
## The which() calls here should probably be omitted for speed
## The tapply() call should be replaced by the required bits of code
##   from tapply()
## sort(uniq()) fails do do the right thing for some tips, where the
##   difference is ~1e-15.  I'm not sure what the best approach there
##   is though.
## If this is a time sink, though, the best option will probably be
## just to cache everything.
##
## 2009/03/9: Changed to interpret NA states as validly missing data,
## and set D_0(0) = D_1(0) = f[i]
bisse.branches <- function(pars, cache) {
  tips <- function(state) {
    if ( is.na(state) )
      i <- setdiff(which(is.na(tip.state)), unresolved$i)
    else
      i <- setdiff(which(tip.state == state), unresolved$i)

    if ( length(i) > 0 ) {
      t <- len[i]
      j <- tapply(t, t)
      t.uniq <- sort(unique(t))
      
      yi <- c(1-sampling.f, 0, 0)
      if ( is.na(state) )
        yi[3:4] <- sampling.f
      else
        yi[state + 3] <- sampling.f[state + 1]
      
      ans <- t(solve(yi, t.uniq, pars))
      n <- length(i)
      branch.init[i,] <<- matrix(rep(yi, n), n, 4, TRUE)
      branch.base[i,] <<- ans[j,-1]
    }
  }
  unresolved.tips <- function() {
    i <- unresolved$i
    Nc <- unresolved$Nc
    k <- unresolved$k
    nsc <- unresolved$nsc
    t <- len[i]
    nt.extra <- cache$nt.extra
    branch.base[i,] <<- solve.unresolved(Nc, k, t, pars, nsc, nt.extra)
  }

  len <- cache$len
  tip.state <- cache$tip.state
  children <- cache$children
  order <- cache$order
  root <- cache$root
  unresolved <- cache$unresolved
  nnode <- length(order)
  sampling.f <- cache$sampling.f

  ## Set up required space
  n <- length(len)
  branch.init <- branch.base <- matrix(NA, n, 4)
  lq <- rep(0, n)

  if ( !is.null(unresolved) )
    unresolved.tips()
  tips(0)
  tips(1)
  tips(NA)

  ## This does all the branches down which we integrate
  for ( i in order[-nnode] ) {
    y.in <- initial.conditions(branch.base[children[i,],], pars)
    lq[i] <- y.in[5]
    y.in <- y.in[-5]
    if ( any(is.na(y.in) | y.in < 0) )
      stop("Invalid conditions")
    branch.init[i,] <- y.in
    branch.base[i,] <- solve(y.in, len[i], pars)[-1]
  }

  ## The final node join is best done without the speciation event
  y.in <- initial.conditions(branch.base[children[root,],], pars, TRUE)
  lq[root] <- y.in[5]
  branch.init[root,] <- y.in[-5]
  list(init=branch.init, base=branch.base, lq=lq)
}

## Given the conditions at the base of a nodes two daughter branches,
## construct the initial conditions for a new node (as well as
## computing the underflow compensation)
initial.conditions <- function(init, pars, is.root=FALSE) {
  e <- init[1,c(1,2)]
  d <- init[1,c(3,4)] * init[2,c(3,4)]
  if ( !is.root )
    d <- d * pars[c(1,2)]
  q <- min(d)
  if ( q != 0 )
    c(e, d/q, log(q))
  else
    c(e, d, 0)
}

stationary.freq <- function(pars) {
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

quadratic.roots <- function(a, b, c)
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)

bisse.prior <- function(pars, r)
  - sum(pars * r)

solve.unresolved <- function(Nc, k, t, pars, nsc=Nc, nt.extra=10,
                             nt=max(Nc)+nt.extra) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  bucexpl(nt, mu0, mu1, lambda0, lambda1, q01, q10, t,
          Nc, nsc, k)[,c(3:4,1:2)]
}
