## My simple-minded EB calculator.

## This is a bit of a trick, because we need to augment everything
## with the depth.

## 1: make
make.eb <- function(tree, states, states.sd=0, control=list()) {
  control <- check.control.continuous(control)
  cache <- make.cache.eb(tree, states, states.sd, control)

  if (control$method == "vcv") {
    all.branches <- make.all.branches.eb.vcv(cache, control)
    rootfunc <- rootfunc.bm.vcv
  } else {
    all.branches <- make.all.branches.eb.pruning(cache, control)
    rootfunc <- rootfunc.bm.pruning
  }

  ## Next, find the root depth:
  max.t <- max(cache$depth)

  ll <- function(pars, root=ROOT.MAX, root.x=NULL,
                 intermediates=FALSE) {
    check.pars.eb(pars)
    ans <- all.branches(c(pars, max.t), intermediates)
    rootfunc(ans, pars, root, root.x, intermediates)
  }
  class(ll) <- c("eb", "dtlik", "function")
  ll
}

## 2: info
make.info.eb <- function(phy) {
  list(name="eb",
       name.pretty="Early Burst (AC/DC)",
       ## Parameters:
       np=2L,
       argnames=default.argnames.eb(),
       ## Variables:
       ny=3L,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=FALSE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Blomberg et al. 2003. Evolution 57: 717"))
}
default.argnames.eb <- function()
  c("s2", "a")

## 3: make.cache
make.cache.eb <- function(tree, states, states.sd, control) {
  cache <- make.cache.bm(tree, states, states.sd, control)
  cache$info <- make.info.eb(tree)
  cache
}

###########################################################################
## Additional functions

## Checking
check.pars.eb <- function(pars) {
  if (length(pars) != 2)
    stop("Incorrect parameter length")
  check.nonnegative(pars[1])
  check.nonpositive(pars[2])
  TRUE
}

## VCV approach:
make.all.branches.eb.vcv <- function(cache, control) {
  n.tip <- cache$n.tip
  states <- cache$states
  states.sd <- cache$states.sd

  eb.rescale <- make.eb.rescale(cache$info$phy)

  function(pars, intermediates) {
    s2 <- pars[1]
    a <- pars[2]

    vcv <- vcv.phylo(eb.rescale(a))
    vv <- s2 * vcv
    # Below here is identical to make.all.branches.ou.vcv
    diag(vv) <- diag(vv) + states.sd^2
    mu <- phylogMean(vv, states)
    dmvnorm2(states, rep(mu, n.tip), vv, solve(vv), log=TRUE)
  }
}

## Q&D reimplementation of rescale.phylo(EB)
make.eb.rescale <- function(phy) {
  ht <- heights.phylo(phy)
  N <- Ntip(phy)
  Tmax <- ht$start[N + 1]
  mm <- match(seq_len(nrow(ht)), phy$edge[, 2])
  ht$t1 <- Tmax - ht$end[phy$edge[mm, 1]]
  ht$t2 <- ht$start - ht$end + ht$t1

  ## NOTE: I think that this is in terms of *log* a.  This might be an
  ## issue with some of the other parameters, too.  In the Blomberg
  ## formulation, '1' is the switch point between accelerate and
  ## decelerate, but here it's log(1).  This is also consistent with
  ## the rewrite that I found from last week, where the definite
  ## integral of the Blomberg function (and end up with something
  ## like:
  ##   (g^t2 - g^t1)/log(g); a = log(g)
  function(a) {
    if (a == 0) 
      return(phy)
    ## OK, this little bit here is what we need to actually get.  This
    ## should be really easy as we have depth (t0) and len.
    bl <- (exp(a * ht$t2) - exp(a * ht$t1))/a
    phy$edge.length <- bl[phy$edge[, 2]]
    phy
  }
}

heights.phylo <- function(phy) {
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)
  n.node <- phy$Nnode
  xx <- numeric(n + n.node)
  for (i in nrow(phy$edge):1)
    xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
  root <- ifelse(is.null(phy$root.edge), 0, phy$root.edge)
  labs <- c(phy$tip.label, phy$node.label)
  depth <- max(xx)
  tt <- depth - xx
  idx <- 1:length(tt)
  dd <- phy$edge.length[idx]
  mm <- match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
  dd <- c(phy$edge.length, root)[mm]
  ss <- tt + dd
  res <- cbind(ss, tt)
  rownames(res) <- idx
  colnames(res) <- c("start", "end")
  res <- data.frame(res)
  res
}

make.all.branches.eb.pruning <- function(cache, control) {
  if (control$backend == "R")
    function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache,
                          initial.conditions.bm.pruning,
                          branches.eb, preset)
  else {
    ## NOTE: This is a hack, but allow here for the extra paramter:
    cache$info$np <- 3L
    make.all.branches.continuous(cache, control)
  }
}


## The issue that I have here is that time is computed against the
## root, so we'll need to know that.  t0 is the *depth* of the
## branch *tip* so t0+len is the *depth* of the branch base:
##    tr    t0+len         t0      0
##    |-----|--------------|-------|
##    0     s0             s1      st
## So let, st = tr (time of root, time of tip), then
##    s1 = tr - t0
##    s0 = tr - (t0 + len) = s1 - len
branches.eb <- function(y, len, pars, t0, idx) {
  m <- y[[1]]
  v <- y[[2]]
  z <- y[[3]]

  sigma2 <- pars[[1]]
  a      <- pars[[2]]
  tr     <- pars[[3]]

  if (a != 0) {
    s1 <- tr - t0
    s0 <- s1 - len
    len <- (exp(a * s1) - exp(a * s0))/a
  }

  list(z, c(m, v + sigma2 * len, 0))
}
