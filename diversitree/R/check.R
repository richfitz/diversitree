## Checking utilities.  These are things that happen over and over
## again, and are tedious to have to write into each function.
check.tree <- function(tree, ultrametric=TRUE, bifurcating=TRUE,
                       node.labels=FALSE) {
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

check.states <- function(tree, states, allow.unnamed=FALSE,
                         strict=FALSE, strict.vals=NULL) {
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

  if ( strict && !is.null(strict.vals) )
    if ( !isTRUE(all.equal(sort(strict.vals),
                           sort(unique(na.omit(states))))) )
      stop("Because strict state checking requested, all (and only) ",
           sprintf("states in %s are allowed",
                   paste(strict.vals, collapse=", ")))

  if ( inherits(tree, "clade.tree") ) {
    spp.clades <- unlist(tree$clades)
    if ( !all(spp.clades %in% names(states)) )
      stop("Species in 'clades' do not have states information")
    states[union(tree$tip.label, spp.clades)]
  } else {
    states[tree$tip.label]
  }
}

check.par.length <- function(x, length) {
  if ( length(x) == 1 )
    rep(x, length)
  else if ( length(x) == length )
    x
  else
    stop(sprintf("'%s' of incorrect length",
                 deparse(substitute(x))))
}

check.sampling.f <- function(sampling.f, n) {
  if ( is.null(sampling.f) )
    sampling.f <- rep(1, n)
  else
    sampling.f <- check.par.length(sampling.f, n)

  if ( max(sampling.f) > 1 || min(sampling.f) <= 0 )
    stop("sampling.f must be on range (0,1]")
  sampling.f
}

check.sampling.f.split <- function(sampling.f, n, n.part) {
  if ( is.null(sampling.f) )
    rep(list(rep(1, n)), n.part)
  else if ( is.numeric(sampling.f) )
    rep(list(check.sampling.f(sampling.f, n)), n.part)
  else if ( is.list(sampling.f) )
    lapply(sampling.f, check.sampling.f, n)
  else
    stop("Invalid sampling.f")
}

check.bounds <- function(lower, upper, x0=NULL) {
  if ( !is.null(x0) && (any(x0 < lower) || any(x0 > upper)) )
    stop("Starting parameter falls outside of problems bounds")
  if ( any(lower >= upper) )
    stop("'upper' must be strictly greater than 'lower'")
}

check.par.multipart <- function(pars, n.part, n.per) {
  if ( is.matrix(pars) ) {
    if ( nrow(pars) != n.part )
      stop(sprintf("Expected %d parameter sets", n.part))
    if ( ncol(pars) != n.per )
      stop(sprintf("Expected %d parameters in each set", n.per))
    pars <- matrix.to.list(pars)
  } else if ( is.list(pars) ) {
    if ( length(pars) != n.part )
      stop(sprintf("Expected %d parameter sets", n.part))
    if ( !all(unlist(lapply(pars, length)) == n.per) )
      stop(sprintf("Expected %d parameters in each set", n.per))      
  } else {
    if ( length(pars) != n.part * n.per )
      stop(sprintf("Expected %d parameters", n.part * n.per))
    pars <- matrix.to.list(matrix(pars, n.part, n.per, TRUE))
  }
  pars
}

## Check that a number can reasonably be considered an integer.
check.integer <- function(x) {
  if ( max(abs(x - round(x))) > 1e-8 )
    stop("Non-integer argument for ", deparse(substitute(x)))
  as.integer(x)
}
