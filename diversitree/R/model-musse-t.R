## Models should provide:
##   1. make
##   2. print
##   3. argnames / argnames<-
##   4. find.mle
## Generally, make will require:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches

## MuSSE/t is a bunch harder than the basic BiSSE/t or BD/t models, as
## the diagonal of the Q matrix needs to be specified; this is
## generated from the given parameters for each parameter update.
##
## This is easiest to deal with for models where all the Q elements do
## not change (and it might be worth having a special case for that).
## However, I think that the best bet is to intercept the functions
## and generate the Q matrix ourselves.  This might violate the
## assumptions of one element per function (and then require unlist()
## again).
make.musse.t <- function(tree, states, k, functions, sampling.f=NULL,
                         strict=TRUE, control=list()) {
  control <- modifyList(list(safe=FALSE, tol=1e-8, eps=0), control)
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)

  if ( is.null(names(functions)) && length(functions) == k*(k+1) )
    names(functions) <- argnames.musse(NULL, k)

  pars.t <- make.pars.t.musse(functions, k)
  n.args <- attr(pars.t, "n.args")
  is.constant.arg <- attr(pars.t, "is.constant.arg")

  branches <- make.branches.musse.t(k, control$safe, control$tol,
                                    control$eps)
  initial.conditions <-
    make.initial.conditions.t(initial.conditions.musse)

  ll.musse.t <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                         root.p=NULL, intermediates=FALSE) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid length parameters (expected %d)", n.args))
    pars.const <- pars[is.constant.arg]
    if ( any(pars.const < 0) || any(!is.finite(pars.const)) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    f.pars <- pars.t(pars)

    ll.xxsse.t(f.pars, cache, initial.conditions, branches,
               condition.surv, root, root.p, intermediates)
  }

  class(ll.musse.t) <- c("musse.t", "musse", "function")
  attr(ll.musse.t, "argnames") <- attr(pars.t, "argnames")
  ll.musse.t
}

make.pars.t.musse <- function(functions, k) {
  if ( length(functions) != k * (k + 1) )
    stop("Wrong number of functions")

  obj <- check.functions.t(functions)

  n.args <- obj$n.args
  idx <- obj$idx
  is.constant <- obj$is.constant
  idx.constant <- obj$idx.constant
  i.var <- which(!is.constant)

  ## This includes the diagonal elements of the Q matrix.
  out <- numeric(k * (k + 2))

  is.constant.f <- obj$is.constant

  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.q <- (2*k+1):(k*(1+k))

  if ( !all(is.constant.f[idx.q]) ) {
    stop("Time varying Q matrix not yet dealt with")
  }

  idx.constant.lm <- unlist(idx[intersect(1:(2*k), which(is.constant.f))])
  idx.constant.q <- unlist(idx[intersect(idx.q, which(is.constant.f))])
  idx.constant.target <- c(which(is.constant[1:(2*k)]), (2*k+1):(k*(k+2)))
  
  ret <- function(pars) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid argument length: expected %d", n.args))
    names(pars) <- NULL # because of do.call

    qmat[idx.qmat] <- pars[idx.constant.q]
    diag(qmat) <- -rowSums(qmat)
    out[idx.constant.target] <- c(pars[idx.constant.lm], qmat)

    function(t) {
      ## Surprisingly, this for loop was faster than lapply.
      for ( i in i.var )
        out[[i]] <- do.call(functions[[i]],
                            c(list(t), pars[idx[[i]]]))
      out
    }
  }

  attr(ret, "n.args") <- n.args
  attr(ret, "argnames") <- obj$argnames
  attr(ret, "is.constant.f") <- is.constant
  attr(ret, "is.constant.arg") <- obj$is.constant.arg

  ret
}

`argnames<-.musse.t` <- function(x, value) {
  .NotYetImplemented()
}

## 8: branches
## Note that this is quite different to the normal
## make.musse.branches, as the assembly of the Q matrix is done before
## this gets here.  Normally this is done on the way in (TODO:
## arguably, it should be done in the likelihood funtion, and this may
## change in the future)
make.branches.musse.t <- function(k, safe=FALSE, tol=1e-8, eps=0) {
  RTOL <- ATOL <- tol
  e <- new.env()
  
  musse.t.ode <- make.ode("derivs_musse_t", "diversitree",
                          "initmod_musse_t", 2*k, safe)
  branches <- function(y, len, pars, t0)
    t(musse.t.ode(y, c(t0, t0+len), list(pars, e),
                  rtol=RTOL, atol=ATOL)[-1,-1])
  make.branches(branches, (k+1):(2*k), eps)
}
