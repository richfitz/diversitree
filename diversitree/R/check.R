## Checking utilities.  These are things that happen over and over
## again, and are tedious to have to write into each function.
check.tree <- function(tree, ultrametric=TRUE, bifurcating=TRUE,
                       node.labels=TRUE) {
  if ( !inherits(tree, "phylo") )
    stop("'tree' must be a valid phylo tree")
  if ( ultrametric && !is.ultrametric(tree) )
    stop("'tree' must be ultrametric")
  ## ape's is.binary.tree() can let a few nasties through - for
  ## e.g. each tritomy, an unbranched node and this gets through.
  ## This expression is a little stricter, even if a touch slower.
  if ( bifurcating && (!is.binary.tree(tree) ||
                       any(tabulate(tree$edge[, 1]) == 1)) )
    stop("'tree must be bifurcating (no polytomies or unbranched nodes)'")

  if ( any(duplicated(tree$tip.label)) )
    stop("Tree contains duplicated tip labels")
  
  if ( node.labels ) {
    if ( is.null(tree$node.label) )
      tree$node.label <- sprintf("nd%d", seq_len(tree$Nnode))
    else if ( any(duplicated(tree$node.label)) )
      stop("Tree contains duplicated node labels")
  }

  tree
}

check.states <- function(tree, states, allow.unnamed=FALSE) {
  if ( is.null(names(states)) ) {
    if ( allow.unnamed ) {
      if ( length(states) == length(tree$tip.label) ) {
        names(states) <- tree$tip.label
        warning("Assuming states are in tree$tip.label order")
      } else {
        stop(sprintf("Invalid states length (expected %d)",
                     length(tree$tip.label)))
      }
    } else {
      stop("The states vector must contain names")
    }
  }
  
  if ( !all(tree$tip.label %in% names(states)) )
    stop("Not all species have state information")

  states[tree$tip.label]
}

check.sampling.f <- function(sampling.f, n) {
  if ( is.null(sampling.f) )
    sampling.f <- rep(1, n)
  else if ( length(sampling.f) != n )
    stop(sprintf("sampling.f must be of length %d (or NULL)", n))
  else if ( max(sampling.f) > 1 || min(sampling.f) <= 0 )
    stop("sampling.f must be on range (0,1]")
  sampling.f
}

check.unresolved <- function(tree, unresolved, nt.extra) {
  if ( !is.null(unresolved) && nrow(unresolved) == 0 ) {
    unresolved <- NULL
    warning("Ignoring empty 'unresolved' argument")
  } else if ( !is.null(unresolved) ) {
    required <- c("tip.label", "Nc", "n0", "n1")
    if ( !all(required %in% names(unresolved)) )
      stop("Required columns missing from unresolved clades")
    unresolved$tip.label <- as.character(unresolved$tip.label)
    if ( !all(unresolved$tip.label %in% tree$tip.label) )
      stop("Unknown tip species in 'unresolved'")
    unresolved$k   <- unresolved$n1
    unresolved$nsc <- unresolved$n0 + unresolved$n1
    unresolved$i   <- match(unresolved$tip.label, tree$tip.label)
    unresolved <- as.list(unresolved)
    unresolved$nt.extra <- nt.extra

    if ( max(unresolved$Nc + nt.extra) > 200 )
      stop("The largest unresolved clade supported has %d species",
           200 - nt.extra)
  }

  unresolved
}
