## TODO:
## One complication here is that without evaluating the functions
## (which requires an initial position) we have no way of knowing how
## many parameters should/will be created.  This will need to be added
## to the C code.

## There is a chance here that it will end up being just as efficient
## (and a bit clearer) to have the calculations done purely in R.  I
## should check this at some point.

## Time dependent code.  It will be useful to have a couple of
## primatives here.
constant.t <- function(t, c) rep.int(c, length(t))
linear.t <- function(t, c, m) {
  out <- c + m * t
  out[out < 0] <- 0
  out
}
sigmoid.t <- function(t, y0, y1, tmid, r)
  y0 + (y1 - y0)/(1 + exp(r * (tmid - t)))
exponential.t <- function(t, y0, y1, r)
  y0 + (y1 - y0) * exp(- r * t)
stepf.t <- function(t, y0, y1, tc)
  if ( t <= tc ) y0 else y1

check.f.t <- function(f) {
  args <- names(formals(f))
  if ( args[1] != "t" )
    stop("First argument of time-dependent function must be t")
  length(args) - 1
}

## TODO: I should assert somewhere that the return values of the
## functions are scalars.  This is enforced by the out[[i]], but this
## will lead to cryptic error messages.
## make.pars.t <- function(...) {
##   functions <- list(...)
##   n <- sapply(functions, check.f.t)
##   idx <- split(seq_len(sum(n)), rep(seq_along(n), n))
##   n.args <- sum(n)
##   n.types <- length(idx)
##   idx.seq <- seq_along(idx)

##   is.constant <- sapply(functions, identical, constant.t)
##   idx.constant <- unlist(idx[is.constant])

##   is.var <- !is.constant
##   idx.var <- idx[is.var]
##   i.var <- which(is.var)

##   out <- numeric(length(functions))

##   ret <- function(pars) {
##     if ( length(pars) != n.args )
##       stop(sprintf("Invalid argument length: expected %d", n.args))
##     names(pars) <- NULL # because of do.call
##     out[is.constant] <- pars[idx.constant]
    
##     function(t) {
##       ## Surprisingly, this for loop was faster than lapply.
##       for ( i in i.var )
##         out[[i]] <- do.call(functions[[i]],
##                             c(list(t), pars[idx[[i]]]))
##       out
##     }
##   }

##   ## Build argument names.
##   if ( is.null(names(functions)) )
##     stop("functions list must be named")
##   argnames <- vector("list", length(functions))
##   argnames[is.constant] <- names(functions)[is.constant]
##   argnames[is.var] <-
##     lapply(names(functions)[is.var], function(i)
##            paste(i, names(formals(functions[[i]]))[-1], sep="."))
##   argnames <- unlist(argnames)
##   if ( any(duplicated(argnames)) )
##     stop("Duplicate argument names: consider different prefixes?")
  
##   attr(ret, "n.args") <- n.args
##   attr(ret, "argnames") <- argnames
##   attr(ret, "is.constant.f") <- is.constant
##   attr(ret, "is.constant.arg") <- rep(is.constant, n)
       
##   ret
## }

make.initial.conditions.t <- function(initial.conditions) {
  function(init, pars, t, is.root=FALSE)
    initial.conditions(init, pars(t), t, is.root)
}

## This is identical to the version in diversitree-branches.R, except
## for the root parameter treatment.
ll.xxsse.t <- function(pars, cache, initial.conditions,
                       branches, condition.surv, root, root.p,
                       intermediates) {
  pars.root <- pars(cache$depth[cache$root])

  ans <- all.branches(pars, cache, initial.conditions, branches)
  vals <- ans$init[[cache$root]]
  root.p <- root.p.xxsse(vals, pars.root, root, root.p)
  loglik <- root.xxsse(vals, pars.root, ans$lq, condition.surv,
                       root.p)
  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
}

check.functions.t <- function(functions) {
  if ( is.null(names(functions)) )
    stop("functions list must be named")
  
  n <- sapply(functions, check.f.t)
  idx <- split(seq_len(sum(n)), rep(seq_along(n), n))
  n.args <- sum(n)
  n.types <- length(idx)
  idx.seq <- seq_along(idx)

  is.constant <- sapply(functions, identical, constant.t)
  idx.constant <- unlist(idx[is.constant])

  ## Build argument names.
  argnames <- vector("list", length(functions))
  argnames[is.constant] <- names(functions)[is.constant]
  argnames[!is.constant] <-
    lapply(names(functions)[!is.constant], function(i)
           paste(i, names(formals(functions[[i]]))[-1], sep="."))
  argnames <- unlist(argnames)
  if ( any(duplicated(argnames)) )
    stop("Duplicate argument names: consider different prefixes?")

  list(is.constant=is.constant, idx=idx,
       idx.constant=idx.constant,
       n.args=n.args, argnames=argnames,
       is.constant.arg=rep(is.constant, n))
}

make.pars.t <- function(functions) {
  obj <- check.functions.t(functions)

  out <- numeric(length(functions))

  n.args <- obj$n.args
  idx <- obj$idx
  is.constant <- obj$is.constant
  idx.constant <- obj$idx.constant
  i.var <- which(!is.constant)

  ret <- function(pars) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid argument length: expected %d", n.args))
    names(pars) <- NULL # because of do.call
    out[is.constant] <- pars[idx.constant]
    
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

