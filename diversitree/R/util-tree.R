## This is the same as the function in quasse2
descendants <- function(node, edge) {
  ans <- node
  repeat {
    node <- edge[edge[,1] %in% node,2]
    if ( length(node) > 0 )
      ans <- c(ans, node)
    else
      break
  }

  unlist(ans)
}

ancestors <- function(phy, i=seq_along(phy$tip.label)) {
  anc <- i
  edge <- phy$edge
  while ( any(!is.na(i)) ) {
    i <- edge[match(i, edge[,2]),1]
    anc <- cbind(anc, i, deparse.level=0)
  }

  apply(anc, 1, function(x)
        c(rev(x[!is.na(x)]), rep(NA, sum(is.na(x)))))
}

## Compute the MRCA of tips with indices in 'tips'
mrca.tipset <- function(phy, tips) {
  if ( is.character(tips) )
    tips <- match(tips, phy$tip.label)
  
  if ( length(tips) == 1 )
    tips
  else {
    anc <- ancestors(phy, tips)
    j <- which(apply(anc, 1, function(x) length(unique(x))) > 1)[1]
    anc[j-1,1]
  }
}

## Similar to ape's branching.times(), but returns the height above
## the root node, even for non-ultrametric trees.  Includes tip times.
branching.heights <- function(phy) {
  if (!inherits(phy, "phylo"))
    stop('object "phy" is not of class "phylo"')

  edge <- phy$edge
  n.node <- phy$Nnode
  n.tip <- length(phy$tip.label)

  ht <- numeric(n.node + n.tip) # zero'd
  for (i in seq_len(nrow(edge)) )
    ht[edge[i, 2]] <- ht[edge[i, 1]] + phy$edge.length[i]

  ## Ugly, but fairly compatible with branching.times()
  names.node <- phy$node.label
  if (is.null(names.node))
    names.node <- (n.tip + 1):(n.tip + n.node)
  names(ht) <- c(phy$tip.label, names.node)

  ht
}

## This only works for ultrametric trees:
branching.depth <- function(len, children, order, tips) {
  depth <- numeric(nrow(children))
  depth[tips] <- 0
  for ( i in order )
    depth[i] <- depth[children[i,1]] + len[children[i,1]]
  depth
}
