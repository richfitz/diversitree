## Support for trees with unresolved clades

make.clade.tree <- function(tree, clades) {
  if ( !identical(class(tree), "phylo") )
    stop("tree must be a plain 'phylo' tree")
  if ( !all(names(clades) %in% tree$tip.label) )
    stop("Unknown clade representatives")
  if ( !all(sapply(clades, is.character)) )
    stop("'clades' must be a list of character vectors")
  if ( any(duplicated(unlist(clades))) )
    stop("Duplicated species names")
  
  tree$clades <- clades
  class(tree) <- c("clade.tree", class(tree))
  tree
}

## This takes a state vector and a tree that has been split into
## unresolved clades and produces an 'unresolved' data structure.
make.unresolved <- function(clades, states) {
  if ( !all(unlist(clades) %in% names(states)) )
    stop("Species in 'clades' do not have states information")

  count.true  <- function(x) sum(!is.na(x) & x)
  count.false <- function(x) sum(!is.na(x) & !x)
  
  states.clades <- lapply(clades, function(x) states[x])
  unresolved <-
    data.frame(tip.label=names(clades),
               Nc=sapply(clades, length),
               n0=sapply(states.clades, count.false),
               n1=sapply(states.clades, count.true))
  rownames(unresolved) <- NULL
  unresolved
}

polytomies.to.clades <- function(tree) {
  from <- tree$edge[,1]
  to   <- tree$edge[,2]
  n.taxa <- length(tree$tip.label)

  edge.counts <- tapply(to, from, length)
  poly.nodes <- as.integer(names(edge.counts[edge.counts > 2]))

  ## Find who the ancestors of a node are:
  ans <- lapply(poly.nodes, ancestors, tree)

  is.node <- seq_len(max(tree$edge)) %in% from

  ## Now, some of these are nested within one another:
  ans1 <- mapply(setdiff, ans, poly.nodes)
  clades <-  lapply(ans1[!(poly.nodes %in% unlist(ans1))],
                    function(x) x[x <= n.taxa])

  ## One representative species (arbitrarily the first species in the
  ## clade) will be left in the phylogeny.
  clades.repr <- tree$tip.label[sapply(clades, "[",  1)]
  names(clades) <- clades.repr
  clades.spp <- lapply(clades, function(x) tree$tip.label[x])

  ## Drop all but the representative species from each clade.
  clades.drop <- sort(unlist(lapply(clades, "[", -1)))
  tree2 <- drop.tip(tree, clades.drop)

  make.clade.tree(tree2, clades.spp)
}

## All ancestors of node 'x' in 'tree'
## ancestors <- function(x, tree, tips.only=FALSE) {
##   from <- tree$edge[,1]
##   to   <- tree$edge[,2]
##   n.taxa <- length(tree$tip.label)
##   is.node <- seq_len(n.taxa) %in% from
##   anc <- list(x)
##   n <- 1
##   while ( length(x) > 0 ) {
##     kids <- to[from %in% x]
##     anc[[n <- n + 1]] <- kids
##     x <- kids[is.node[kids]]
##   }
##   anc <- sort(unlist(anc))
##   if ( tips.only )
##     anc[anc <= n.taxa]
##   else
##     anc
## }

ancestors <- function(x, tree, tips.only=FALSE) {
  from <- tree$edge[,1]
  to   <- tree$edge[,2]

  is.node <- seq_len(max(tree$edge)) %in% from

  anc <- list(x)
  n <- 1
  while ( length(x) > 0 ) {
    kids <- to[from %in% x]
    anc[[n <- n + 1]] <- kids
    x <- kids[is.node[kids]]
  }

  anc <- sort(unlist(anc))
  if ( tips.only )
    anc[anc <= length(tree$tip.label)]
  else
    anc
}
