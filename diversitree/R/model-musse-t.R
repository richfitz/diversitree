make.musse.t <- function(tree, states, k, functions, sampling.f=NULL,
                         strict=TRUE, control=list()) {
  cache <- make.cache.musse.t(tree, states, k, functions,
                              sampling.f, strict)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.musse)
  rootfunc <- make.rootfunc.t(cache, rootfunc.musse)
  f.pars <- make.pars.t.musse(cache$functions, cache)
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("musse.t", "musse", "dtlik", "function")
  ll
}

make.cache.musse.t <- function(tree, states, k, functions,
                               sampling.f, strict) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  update.cache.t(cache, functions)  
}

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
make.pars.t.musse <- function(functions, cache=NULL,
                        check.negative.const=TRUE,
                        check.negative.var=TRUE) {
  if ( is.null(cache) )
    stop("Cache must be provided for make.pars.musse.t")
  k <- cache$info$k
  if ( length(functions) != k * (k + 1) )
    stop("Wrong number of functions")

  obj <- cache$functions.info

  ## Unpack the function information:
  n.args <- obj$n.args
  idx <- obj$idx
  is.constant <- obj$is.constant
  is.constant.arg <- obj$is.constant.arg # extra to MuSSE
  idx.constant <- obj$idx.constant
  i.var <- which(!is.constant)

  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.q <- (2*k+1):(k*(1+k))
  if ( !all(is.constant[idx.q]) )
    stop("Time varying Q matrix not yet dealt with")
  
  ## This includes the diagonal elements of the Q matrix.
  out <- numeric(k * (k + 2))

  ## Index in the *input* vector of constant lambda/mu parameters
  idx.constant.lm <- unlist(idx[intersect(1:(2*k), which(is.constant))])
  ## Input indices of the q parameters.
  idx.constant.q <- unlist(idx[intersect(idx.q, which(is.constant))])
  ## And corresponding indices of where both go in the *output*
  idx.constant.target <-
    c(which(is.constant[1:(2*k)]), (2*k+1):(k*(k+2)))
  
  function(pars) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid argument length: expected %d", n.args))
    names(pars) <- NULL # because of do.call
    check.nonnegative(pars[is.constant.arg])

    ## Fill in constant times.
    qmat[idx.qmat] <- pars[idx.constant.q]
    diag(qmat) <- -rowSums(qmat)
    out[idx.constant.target] <- c(pars[idx.constant.lm], qmat)

    function(t) {
      for ( i in i.var ) # Loop faster than lapply
        out[[i]] <- do.call(functions[[i]],
                            c(list(t), pars[idx[[i]]]))
      check.nonnegative(out[i.var])      
      out
    }
  }
}

make.branches.musse.t <- function(cache, control)
  make.branches.dtlik.t(cache$info, control)
