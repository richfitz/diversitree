## I'm going to take code here from arbutus, rather than from geiger,
## because we're tracking it more carefully and it's already tested.
## Once arbutus is released, perhaps we can either depend on arbutus
## or we can pull some of these utilities out.

## The difference here is that SE rescaling is handled separately.  We
## avoid scaling by s^2 (that's for the BM calculations to deal
## with).  Parameter vectors are handled differently (and assumed
## correct).
branching.times.edge <- function(phy) {
  ht <- unname(branching.heights(phy)) # same as node.depth.edgelength()
  list(start=ht[phy$edge[,1]],
       end=ht[phy$edge[,2]],
       max=max(ht))
}

make.rescale.phylo.ou <- function(phy) {
  t <- branching.times.edge(phy)
  t.start <- t$start
  t.end   <- t$end
  t.max   <- t$max

  ## This should also be easy to port over to a theta-less pruning
  ## approach.
  function(alpha) {
    if (alpha > 0) {
      tau.max   <- 2 * alpha * t.max
      tau.start <- 2 * alpha * t.start
      tau.end   <- 2 * alpha * t.end
      phy$edge.length <-
        exp(-tau.max) / (2 * alpha) * (
          exp(tau.start) * expm1(-tau.start) -
          exp(tau.end)   * expm1(-tau.end))
    }
    phy
  }
}

make.rescale.phylo.eb <- function(phy, pars) {
  t <- branching.times.edge(phy)
  t.start <- t$start
  t.end   <- t$end
  t.max   <- t$max

  function(a) {
    if (a != 0) {
      phy$edge.length <-
        (exp(a * t.end) - exp(a * t.start)) / a
    }
    phy
  }
}

make.rescale.phylo.lambda <- function(phy) {
  t <- branching.times.edge(phy)
  t.start <- t$start
  t.end   <- t$end
  t.max   <- t$max

  n.tip <- length(phy$tip.label)
  interns <- which(phy$edge[, 2] >  n.tip)
  externs <- which(phy$edge[, 2] <= n.tip)

  t <- branching.times.edge(phy)
  height <- t$end[externs]

  function(lambda) {
    phy$edge.length[interns] <-
      phy$edge.length[interns] * lambda
    # This *might* be correct:
    phy$edge.length[externs] <-
      phy$edge.length[externs] * lambda + (1 - lambda) * height
    phy
  }
}

rescale.phylo.se <- function(phy, se) {
  tips <- phy$edge[,2] <= length(phy$tip.label)
  phy$edge.length[tips] <- phy$edge.length[tips] + se
  phy
}
