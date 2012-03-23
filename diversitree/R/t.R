make.all.branches.t.dtlik <- function(cache, control,
                                      initial.conditions) {
  control <- check.control.ode(control)

  if ( control$backend == "CVODES" ) {
    stop("CVODES not yet available for time-varying models")
  } else {
    branches.t <- make.branches.dtlik.t(cache$info, control)    
    initial.conditions.t <-
      make.initial.conditions.t(initial.conditions)
    function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache, initial.conditions.t,
                          branches.t, preset)
  }
}

make.branches.dtlik.t <- function(info, control) {
  info <- check.info.ode(info, control)
  br <- make.branches.dtlik(info, control)
  dll <- info$dll
  e <- new.env()

  if ( control$backend == "deSolve" ) {
    function(y, len, pars, t0, idx)
      br(y, len, list(pars, e), t0, idx)
  } else {
    ## I can't quite remember why this works...
    setfunc <- sprintf("r_set_tfunc_%s", info$name)
    dummy <- rep(0.0, info$np)
    function(y, len, pars, t0, idx) {
      .Call(setfunc, pars, e, PACKAGE=dll)
      br(y, len, dummy, t0, idx)
    }
  }
}

make.initial.conditions.t <- function(initial.conditions) {
  function(init, pars, t, idx)
    initial.conditions(init, pars(t), t, idx)
}

## TODO/NEW This is somewhat tedious, as we really should check for
## 'root' not being ROOT.EQUI, as that can't be done in time-dependent
## models.  However, because this is used by different functions with
## different argument lists, that is hard to do.  But then, that check
## is duplicated in too many functions.  For now, I'm skipping this.
make.rootfunc.t <- function(cache, rootfunc) {
  t.root <- cache$depth[cache$root]
  function(ans, pars, ...)
    rootfunc(ans, pars(t.root), ...)
}

## Compared with update.*.td, this is cache, not info...
update.cache.t <- function(cache, functions) {
  info <- cache$info
  n.args <- length(info$argnames)

  if ( is.function(functions) )
    functions <- rep(list(functions), n.args)
  if ( is.null(names(functions)) && length(functions) == n.args )
    names(functions) <- info$argnames
  cache$functions.info <- check.functions.t(functions) 
  cache$functions <- functions

  info$time.varying <- TRUE
  info$argnames <- cache$functions.info$argnames
  info$name.ode <- sprintf("%s_t", cache$info$name)
  info$name.pretty <- sprintf("%s (time-varying)", info$name.pretty)
  info$name <- sprintf("%s.t", cache$info$name)

  cache$info <- info
  cache
}

######################################################################
## Sample time-dependent functions:
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

######################################################################
## Parameter generator function:
make.pars.t <- function(functions, cache=NULL,
                        check.negative.const=TRUE,
                        check.negative.var=TRUE) {
  if ( is.null(cache) )
    obj <- check.functions.t(functions)
  else
    obj <- cache$functions.info

  ## Unpack the function information:
  n.args <- obj$n.args
  idx <- obj$idx
  is.constant <- obj$is.constant
  idx.constant <- obj$idx.constant
  i.var <- which(!is.constant)

  ## TODO: Possibly I can get musse working if I change this to
  ## info$np? (requires cache)
  out <- numeric(length(functions))

  function(pars) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid argument length: expected %d", n.args))
    names(pars) <- NULL # because of do.call
    out[is.constant] <- pars[idx.constant]

    if ( check.negative.const )
      check.nonnegative(out[is.constant])
    
    function(t) {
      for ( i in i.var ) # Loop faster than lapply.
        out[[i]] <- do.call(functions[[i]],
                            c(list(t), pars[idx[[i]]]))
      check.nonnegative(out[i.var])
      out
    }
  }
}

######################################################################
## Checking
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

check.f.t <- function(f) {
  args <- names(formals(f))
  if ( args[1] != "t" )
    stop("First argument of time-dependent function must be t")
  length(args) - 1
}
