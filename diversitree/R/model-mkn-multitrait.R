## This differs from normal mkn in the interpretation of the arguments
## only; the actual calculation of the values is passed off to mkn, as
## usual...

## 1: make
make.mkn.multitrait <- function(tree, states, depth=NULL,
                                allow.multistep=FALSE,
                                strict=TRUE,
                                control=list()) {
  ## TODO: This cannot be hard to do, if I can do it in MuSSE!
  if ( any(is.na(states)) )
    stop("Missing data not yet allowed")

  n.trait <- ncol(states)
  k <- 2^n.trait

  if ( is.null(control$method) )
    control$method <- if (n.trait > 3) "ode" else "mkn"
  control <- check.control.mkn(control, k)
  
  states <- check.states.musse.multitrait(tree, states, strict=strict,
                                          strict.vals=0:1)

  code <- as.matrix(states)
  storage.mode(code) <- "character"
  code[is.na(code)] <- "."
  code <- apply(code, 1, paste, collapse="")
  types <- unique(code)
  
  ## This is the "key"; the mapping from a series of binary traits
  ## onto the Mk mapping
  key <- apply(do.call(expand.grid, rep(list(0:1), n.trait)),
               1, paste, collapse="")

  states.mkn <- match(code, types)
  names(states.mkn) <- rownames(states)

  tr <- mkn.multitrait.translate(n.trait, depth,
                                 colnames(states),
                                 allow.multistep)

  lik <- make.mkn(tree, states.mkn, k, FALSE, control)

  ll.mkn <- function(pars, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE, pars.only=FALSE) {
    pars2 <- drop(tr %*% pars)
    if ( pars.only )
      pars2
    else
      lik(pars2, root=root, root.p=root.p,
          intermediates=intermediates)
  }

  class(ll.mkn) <- c("mkn.multitrait", "mkn", "function")
  attr(ll.mkn, "k") <- k
  attr(ll.mkn, "n.trait") <- n.trait
  ll.mkn
}

## 2: print
print.mkn.multitrait <- function(x, ...) {
  cat("Mkn (Multitrait) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.mkn.multitrait <- function(x, k=attr(x, "k"), ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    colnames(environment(x)$tr)
  } else {
    ret
  }
}
`argnames<-.mkn.multitrait` <- function(x, value) {
  stop("Not yet implemented")
}

## 4: find.mle
## I think that nlminb is not doing a good job with complicated models
## here.
find.mle.mkn.multitrait <- function(func, x.init, method,
                                    fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.mkn")
}
## Zero is no longer an appropriate lower bound, as the effects can
## just as easily be negative.  We will just have to find bounds by
## trial and error; the function won't allow movement of the
## underlying basis into the negative region.
mcmc.mkn.multitrait <- mcmc.default

## 5: make.cache

## 6: ll is done within make.musse.multitrait

## 7: initial.conditions (unchanged from mkn)

## 8: branches (unchanged from musse)

## Parameter translation
mkn.multitrait.translate <- function(n.trait, depth=NULL,
                                     names=NULL,
                                     allow.multistep=FALSE) {
  if ( is.null(names) )
    names <- LETTERS[seq_len(n.trait)]

  if ( is.null(depth) )
    depth <- n.trait - 1
  else if ( length(depth) != 1 )
    stop("Depth must be 1")
  else if ( depth < 0 )
    stop("'depth' must be nonnegative")
  else if ( depth > n.trait - 1 )
    stop("requested depth too large")

  reparam.q(n.trait, depth, names, allow.multistep)
}
