make.time.machine <- function(functions, t.range, np.out, nonnegative,
                              spline.data=NULL) {
  info.t <- time.machine.types()
  types <- match(functions, names(info.t))
  if ( any(is.na(types)) )
    stop("Unknown time-varying function")

  if ( any(functions == "spline.t") )
    spline.data <- check.spline.data(spline.data, t.range)
  else
    spline.data <- NULL
  
  n <- sapply(info.t[types], length)
  ## Get the starting position
  start <- cumsum(n) - n + 1

  argnames <- mapply(paste, names(functions), info.t[types],
                     sep=".", SIMPLIFY=FALSE)
  is.constant <- functions == "constant.t"
  argnames[is.constant] <- names(functions)[is.constant]
  argnames <- unlist(argnames)
  names(argnames) <- NULL
  if ( any(duplicated(argnames)) )
    stop("Duplicate argument names: consider different prefixes?")

  if ( is.null(nonnegative) )
    nonnegative <- TRUE
  else if ( !is.logical(nonnegative) )
    stop("nonnegative must be NULL or logical")
  ## Then check/expand the values:
  if ( length(nonnegative) == 1 )
    nonnegative <- rep(nonnegative, np.out)
  else if ( length(nonnegative) != np.out )
    stop("Invalid length for nonnegative")

  np.in <- length(argnames)  
  ret <- list(np.in=np.in,
              np.out=as.integer(np.out),
              types=types,
              start=start,
              argnames=argnames,
              nonnegative=nonnegative,
              t.range=t.range,
              spline.data=spline.data)
  ret.C <- ret
  ret.C$types <- toC.int(ret$types)
  ret.C$start <- toC.int(ret$start)
  ptr <- ret$ptr <-
    .Call("r_make_time_machine", ret.C, PACKAGE="diversitree")

  ret$set <- function(pars) {
    if ( length(pars) != np.in )
      stop(sprintf("Invalid argument length: expected %d", np.in))
    .Call("r_init_time_machine", ptr, pars, PACKAGE="diversitree")
    invisible(TRUE)
  }

  ret$get <- function(t)
    .Call("r_run_time_machine", ptr, t, PACKAGE="diversitree")
  ret$getv <- function(t) t(sapply(t, ret$get))
  
  ret
}

## A special function will still be needed for cvodes and CVODES, as
## the normal branches function possibly won't work?
make.all.branches.t2.dtlik <- function(cache, control,
                                       initial.conditions) {
  control <- check.control.ode(control)
  time.machine <- cache$time.machine
  tm.ptr <- time.machine$ptr

  if ( control$backend == "CVODES" ) {
    stop("CVODES not yet available for time-varying models") 
  } else if ( control$backend == "cvodes" ) {
    stop("cvodes not yet available for time-varying models")
  } else {
    branches <- make.branches.dtlik(cache$info, control)
    initial.conditions.t <-
      make.initial.conditions.t2(cache, initial.conditions)
    function(pars, intermediates, preset=NULL) {
      time.machine$set(pars)
      all.branches.matrix(tm.ptr, cache, initial.conditions.t,
                          branches, preset)
    }
  }
}

update.cache.t2 <- function(cache, functions, spline.data) {
  info <- cache$info
  n.args <- length(info$argnames)

  ## This could be check.functions.t2():
  if ( !is.character(functions) )
    stop("'functions' must be characters [will relax soon]")
  if ( length(functions) == 1L )
    functions <- rep(list(functions), n.args)
  if ( is.null(names(functions)) && length(functions) == n.args )
    names(functions) <- info$argnames

  t.range <- range(0, cache$depth[cache$root])
  cache$time.machine <- make.time.machine(functions, t.range, n.args,
                                          cache$nonnegative,
                                          spline.data)

  info$time.varying <- TRUE
  info$argnames <- cache$functions.info$argnames
  info$name.ode <- sprintf("%s_t2", cache$info$name)
  info$name.pretty <- sprintf("%s (time-varying[v2])", info$name.pretty)
  info$name <- sprintf("%s.t2", cache$info$name)

  cache$info <- info
  cache
}

make.initial.conditions.t2 <- function(cache, initial.conditions) {
  pars.t <- cache$time.machine$get
  function(init, pars, t, idx)
    initial.conditions(init, pars.t(t), t, idx)
}

## TODO/NEW This is somewhat tedious, as we really should check for
## 'root' not being ROOT.EQUI, as that can't be done in time-dependent
## models.  However, because this is used by different functions with
## different argument lists, that is hard to do.  But then, that check
## is duplicated in too many functions.  For now, I'm skipping this.
make.rootfunc.t2 <- function(cache, rootfunc) {
  pars.t <- cache$time.machine$get
  t.root <- cache$depth[cache$root]
  function(ans, pars, ...)             # pars here is ignored...
    rootfunc(ans, pars.t(t.root), ...) # ...because tm version used.
}

## Information about possible types.
time.machine.types <- function() {
  info.t <- list(constant.t="c",
                 linear.t=c("c", "m"),
                 stepf.t=c("y0", "y1", "tc"),
                 sigmoid.t=c("y0", "y1", "tmid", "r"),
                 spline.t=c("y0", "y1"))
  info.t[.Call("r_get_time_machine_types", PACKAGE="diversitree")]
}

## TODO: Check that we span the limits properly.
check.spline.data <- function(spline.data, t.range) {
  names.spline.data <- c("t", "y", "deriv")  
  if ( is.null(names(spline.data)) && length(spline.data) == 3 )
    names(spline.data) <- names.spline.data
  if ( !is.list(spline.data) || length(spline.data) != 3 )
    stop("spline.data must be a list of length 3")
  if ( !identical(names(spline.data), names.spline.data) )
    stop("names(spline.data) must be ",
         paste(dQuote(names.spline.data), collapse=", "))
  t <- spline.data$t
  y <- spline.data$y
  if ( length(t) < 3 )
    stop("Must have at least three points") # ? or 2?
  if ( length(t) != length(y) )
    stop("Lengths of t and y in spline.data must be equal")
  if ( any(is.na(t)) || any(is.na(y)) )
    stop("Neither t nor y may contain NA values")

  deriv <- spline.data$deriv
  if ( length(deriv) != 1 )
    stop("spline.data$deriv must be an scalar")
  deriv <- check.integer(spline.data$deriv)

  ## Renormalise y data onto [0,1], after projecting what the actual
  ## minimum and maximum values are.
  ## 
  ## TODO: This could be prone to error, especially as I've just
  ## hardcoded a "large" number of x spacing in here.  Ideally we
  ## would be able to provide this with the spline.data?
  r <- range(spline(t, y,
                    xout=seq(t.range[1], t.range[2], length=10001))$y)
  y <- (y - r[1]) / (r[2] - r[1])

  if ( !is.null(t.range) )
    if ( min(t) > min(t.range) || max(t) < max(t.range) ) {
      msg <-
        sprintf("Spline data must span time range: [%s, %s]",
                min(t.range), max(t.range))
      stop(msg)
    }
  
  list(t=as.numeric(t), y=as.numeric(y), deriv=deriv)
}
