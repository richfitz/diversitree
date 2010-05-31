## Checking that is specific to splitting:
check.split <- function(phy, nodes, split.t) {
  if ( length(split.t) != length(nodes) )
    stop("'nodes' and 'split.t' must be the same length")
    
  ## Check that all nodes are ok
  n.tip <- length(phy$tip.label)
  if ( is.character(nodes) )
    nodes <- match(nodes, phy$node.label) + n.tip
  if ( any(is.na(nodes) | nodes <= n.tip+1 | nodes > n.tip + phy$Nnode) )
    stop("Invalid node specification")
  
  ## Check that split times are OK
  t0 <- branching.times(phy)[nodes - n.tip]
  t1 <- t0 + phy$edge.length[match(nodes, phy$edge[,2])]
  split.t[split.t == 0] <- t0
  split.t[split.t == Inf] <- t1
  if ( any(split.t < t0 | split.t > t1) )
    stop("Invalid split time")

  list(nodes=c(NA, nodes), split.t=c(Inf, split.t))
}

## This uses the split generic: x is a phylogeny, f is a vector of
## nodes to split at, drop is ignored (there only for compatibility),
## and split.t is a vector along 'f' with the time that the split
## occurs at.
split.phylo <- function(x, f, drop=FALSE, split.t, ...) {
  tmp <- check.split(x, f, split.t)

  phy <- x
  nodes <- tmp$nodes
  split.t <- tmp$split.t
  
  n.tip <- length(phy$tip.label)
  bt <- branching.times(phy)
  tt <- structure(rep(0, n.tip), names=phy$tip.label) # tip times

  group <- make.split.phylo.vec(phy, nodes)

  ## Find the *parent* of the different groups, so that we can
  ## establish an order to work in:
  parent.node <- phy$edge[match(nodes, phy$edge[,2]),1]
  parent.group <- group[match(parent.node, phy$edge[,2])]

  ## 'extra.tips' is the index for the tips that represent daughter
  ## groups, the index being against 'phy' (it will change in the sub
  ## sub phylogeny)
  extra.tips <- c(NA, sapply(nodes[-1], function(x)
                             min(descendants(x, phy$edge))))
  names(extra.tips) <- phy$tip.label[extra.tips]

  split1 <- function(idx) {
    i <- phy$edge[group == idx,2]
    daughters <- which(parent.group == idx)

    to.keep <- c(extra.tips[daughters], i[i <= n.tip])
    to.drop <- seq_len(n.tip)[-to.keep]
    phy2 <- drop.tip.fixed(phy, to.drop)
    phy2$bt <- bt[phy2$node.label]
    phy2$tt <- tt[phy2$tip.label]

    phy2$parent <- parent.group[idx]
    
    if ( length(daughters) > 0 ) {
      ## Trim the dummy branches
      spp <- names(extra.tips[daughters])
      j <- match(match(spp, phy2$tip.label), phy2$edge[,2])
      phy2$edge.length[j] <- phy2$edge.length[j] - split.t[daughters]
      names(daughters) <- names(extra.tips[daughters])
      phy2$tt[spp] <- split.t[daughters]
      phy2$daughters <- daughters[order(split.t[daughters])]
    }

    if ( !is.na(nodes[idx]) )
      ## Calculate how much additional space needs computing at the end of
      ## the base node.
      phy2$trailing <- split.t[idx] - bt[nodes[idx] - n.tip]

    phy2
  }

  lapply(seq_along(nodes), split1)
}

## Generate a grouping vector from a phylogeny.  This takes a vector
## 'nodes', and classifies the tree into different groups.
make.split.phylo.vec <- function(phy, nodes, group=NULL) {
  edge <- phy$edge

  descendants.idx <- function(node)
    which(edge[,1] == node | edge[,2] %in% descendants(node, edge))
  
  f <- function(node, phy, group) {
    i <- which(edge[,1] == node | edge[,2] %in% descendants(node, edge))
    base <- group[edge[,2] == node]
    group[i[group[i] == base]] <- max(group) + 1
    group
  }

  if ( is.null(group) )
    group <- rep(1, nrow(edge))
  else if ( length(group) != nrow(edge) )
    stop("Invalid length grouping vector")

  for ( i in nodes )
    group <- f(i, phy, group)

  group
}

## This is slightly dodgy, but appears to work.
dt.split.order <- function(daughters, parents) {
  pending <- rep(TRUE, length(daughters))
  order <- integer(0)
  
  while ( any(pending) ) {
    i <- which(pending & sapply(daughters, length) == 0)
    j <- parents[i]
    order <- c(order, i)
    daughters[j] <- lapply(daughters[j], setdiff, i)
    pending[i] <- FALSE
  }
  
  order
}
