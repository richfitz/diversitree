## TODO: document this!
## much of this file is to split functions what diversitree-branches.R
## is to the core models; code for generating caches, and then runing
## integrations given those caches are provided.

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

  ## This is a hack:
  ## Because the root node does not terminate in an edge in the tree,
  ## this needs catching.  There may be a clearer way of doing this.
  ## The root node is always group 1.  I hope.
  if ( length(i <- which(parent.node == n.tip + 1)) > 0 )
    parent.group[i] <- 1

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

all.branches.split <- function(pars, cache, initial.conditions, branches,
                               branches.aux, intermediates=FALSE) {
  ## TODO: 'n' is a stupid name: change to n.part or n.split
  n.part <- cache$n.part
  obj <- res <- vector("list", n.part)
  aux.i <- cache$aux.i # 1, 2 (E0, E1) for BiSSE

  ## TODO: 'order' is confusing with the *other* cache$order
  ## change to order.part
  for ( i in cache$order ) {
    x <- cache$cache[[i]]

    ## It might be the case that auxilliary variables need computing
    ## at splits.  For example, in BiSSE (and friends), when switching
    ## parameter space, rather than take the E values from the
    ## daughter lineage, I need to recompute them to use the new
    ## parameters at this point (a lineage arising just below the
    ## split would have a totally different probability of extinction
    ## than one after the split).  In this case I have to run a new
    ## branch down with the parent set and parent sampling.f.  The 'D'
    ## values from doing this would be ignored.  For mk2 (and friends)
    ## this is not an issue.
    if ( length(x$daughters) > 0 ) {
      y <- res[x$daughters]
      ## TODO:
      ##   x$depth[names(x$daughters)]
      ## might be replaced with
      ##   x$depth[x$daughters.i]
      ## I think?
      x$depth[names(x$daughters)]
      
      aux <- branches.aux(i, x$depth[names(x$daughters)], pars[[i]])
      if ( is.matrix(aux) )
        aux <- matrix.to.list(aux)

      preset <- list(target=x$daughters.i,
                     base=vector("list", length(x$daughters)),
                     lq=numeric(length(x$daughters)))
      
      for ( j in seq_along(x$daughters) ) {
        yj <- res[[x$daughters[j]]]$base
        yj[aux.i] <- aux[[j]]
        idx <- x$daughters.i[j]
        tmp <- branches(yj, x$len[idx], pars[[i]], x$depth[idx])
        preset$lq[j] <- tmp[1]
        preset$base[[j]] <- tmp[-1]
      }

      if ( is.null(x$preset) ) {
        x$preset <- preset
      } else { # TODO: This should be done with something like modifyList
        x$preset$target <- c(x$preset$target, preset$target)
        x$preset$lq     <- c(x$preset$lq,     preset$lq)
        x$preset$base   <- c(x$preset$base,   preset$base)
      }
    }

    obj[[i]] <-
      all.branches(pars[[i]], x, initial.conditions, branches)

    base <- obj[[i]]$init[[x$root]]
    if ( !is.na(cache$parent[i]) ) { # trailing to deal with...
      tmp <- branches(base, x$trailing.len, pars[[i]], x$trailing.t0)
      res[[i]] <- list(lq=tmp[1] + sum(obj[[i]]$lq), base=tmp[-1])
    } else {
      res[[i]] <- list(lq=sum(obj[[i]]$lq), base=base)
    }

    if ( intermediates )
      res[[i]]$intermediates <- obj[[i]]
  }

  res
}

make.cache.split <- function(tree, nodes, split.t) {
  ## TODO: This has the poor effect of not working correctly to create
  ## single branch partitions.  I should be careful about that.  This
  ## is also quite slow, so that making split functions is not very quick.
  subtrees <- split.phylo(tree, nodes, split.t=split.t)

  n.part <- length(subtrees)

  cache <- list()
  cache$cache <- vector("list", n.part)
  cache$tip.tr <- structure(rep(NA, length(tree$tip.label)),
                            names=tree$tip.label)
  cache$parents <- integer(n.part)
  cache$daughters <- vector("list", n.part)

  for ( i in seq_len(n.part) ) {
    tree.sub <- subtrees[[i]]

    res <- make.cache(tree.sub)
    res$tips <- seq_len(res$n.tip) # TODO: remove after update
    res$tip.label <- tree.sub$tip.label
        res$trailing.len <- tree.sub$trailing
    res$trailing.t0 <- max(res$depth)
    res$daughters <- tree.sub$daughters
    n <- length(res$daughters)
    if ( n > 0 ) {
      res$daughters.i <- match(names(res$daughters), tree.sub$tip.label)
      res$n.tip <- res$n.tip - n
      res$tips <- res$tips[-res$daughters.i]
      res$tip.label <- res$tip.label[-res$daughters.i]
    }

    cache$cache[[i]] <- res
    
    cache$parents[i] <- tree.sub$parent
    cache$daughters[i] <- list(res$daughters)
    if ( length(cache$daughters) != n.part )
      browser()
    cache$tip.tr[tree.sub$tip.label[res$tips]] <- i
  }

  cache$n.parts <- length(subtrees)
  cache$order.parts <- dt.split.order(cache$daughters, cache$parents)

  cache
}
