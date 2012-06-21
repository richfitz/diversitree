## TODO:
##   * accept functions, and roll back to compiled versions where
##     available.
##   * think about split functions (global variable screws with this,
##     but shifting to reassigning the pointer at initmod_*_t would
##     get around this.

## I'm a little unhappy as to the amount of tedious book-keeping
## here.  The basic idea is very simple!  However, this is glue to
## some finicky C code, so perhaps it is unavoidable.

## Still pending is a spline+linear and spline+exponential function,
## that will be useful for the climate analysis.  However, that
## minimally affects the R code, aside from the time.machine.types()
## function.

## Another possible inclusion is a generic R function replacement for
## the original time functions (which will be replaced).  Then
## functions can be given as targets.

make.time.machine <- function(functions, t.range, nonnegative=TRUE,
                              spline.data=NULL, k=NULL) {
  types <- check.time.machine.functions(functions)
  
  ## Compute the starting position for each function within the input
  ## parameter vector.
  info.t <- time.machine.types()  
  n <- sapply(info.t[types], length)
  start <- cumsum(n) - n + 1

  ## Argument names.
  argnames <- mapply(paste, names(functions), info.t[types], sep=".",
                     SIMPLIFY=FALSE)
  is.constant <- functions == "constant.t"
  argnames[is.constant] <- names(functions)[is.constant]
  argnames <- unlist(argnames)
  names(argnames) <- NULL
  if ( any(duplicated(argnames)) )
    stop("Duplicate argument names: consider different prefixes?")

  ## Check and expand non-negative checks.
  if ( is.null(nonnegative) )
    nonnegative <- rep(TRUE, length(functions))
  else if ( !is.logical(nonnegative) )
    stop("nonnegative must be logical")
  else if ( length(nonnegative) == 1 )
    nonnegative <- rep(nonnegative, length(functions))
  else if ( length(nonnegative) != length(functions) )
    stop("Invalid length for nonnegative")

  ret <- list(np.in=as.integer(length(argnames)),
              np.out=as.integer(length(functions)),
              types=types,
              start=start,
              argnames=argnames,
              functions=functions,
              funnames=names(functions),
              nonnegative=nonnegative,
              t.range=t.range)

  if ( any(functions %in% c("spline.t", "spline.linear.t")) )
    ret$spline.data <- check.spline.data(ret, spline.data)
  else if ( !is.null(spline.data) )
    warning("Ignoring spline.data -- no spline function specified")

  ret <- check.q.info(ret, k)

  ret$ptr <- ptr <- make.time.machine.ptr(ret)

  ## This interface is still a bit hairy.  Not sure which will stay.
  n.args <- length(argnames)
  ret$set <- function(pars) {
    if ( length(pars) != n.args )
      stop(sprintf("Invalid argument length: expected %d", n.args))
    .Call("r_init_time_machine", ptr, pars, PACKAGE="diversitree")
    invisible(TRUE)
  }

  ## Note that $get() does not set names, as this needs to be fast
  ## (called by initial.conditions.t).
  ret$get <- function(t)
    .Call("r_run_time_machine", ptr, t, PACKAGE="diversitree")

  q.info <- ret$q.info
  ret$getv <- function(t, pars=NULL) {
    if ( !is.null(pars) )
      ret$set(pars)
    ret <- t(sapply(t, ret$get))
    if ( !is.null(q.info) )
      ret <- ret[,-q.info$drop,drop=FALSE]
    colnames(ret) <- names(functions)
    ret
  }

  ret$get1 <- function(p, t, i) {
    ret$set(p)
    ret$get(t)[i]
  }

  ret
}

check.q.info <- function(obj, k) {
  if ( is.null(k) ) {
    obj$target <- seq_len(obj$np.out)
  } else {
    k <- check.integer(k)

    ## Increase the number of output parameters by k:
    obj$np.out <- obj$np.out + k

    idx.qmat <- cbind(rep(1:k, each=k - 1),
                      unlist(lapply(1:k, function(i) (1:k)[-i])))
    obj$target <- c(seq_len(2*k),
                    matrix(seq_len(k*k) + 2*k, k, k)[idx.qmat])

    ## Determine if any 'from' states are constant:
    functions <- obj$functions
    i.q <- (length(functions) - k*(k-1) + 1):length(functions)
    const.q <- apply(matrix(functions[i.q], k-1) == "constant.t", 2, all)

    drop <- as.integer(seq(obj$np.out-k^2 + 1, obj$np.out, by=k+1))
    obj$q.info <- list(k=k, const.q=const.q, drop=drop)
  }
  obj
}

## TODO: probably I should just move to rewriting the .Call() function
## not to take a list.
make.time.machine.ptr <- function(obj) {
  obj.C <- obj
  obj.C$types  <- toC.int(obj$types)
  obj.C$start  <- toC.int(obj$start)
  obj.C$target <- toC.int(obj$target)

  ptr <- obj$ptr <-
    .Call("r_make_time_machine", obj.C, PACKAGE="diversitree")
}

check.spline.data <- function(obj, spline.data) {
  if ( !all(c("t", "y") %in% names(spline.data)) )
    stop("spline.data must contain 't' and 'y'")
  t <- spline.data$t
  y <- spline.data$y
  if ( length(t) < 3 )
    stop("Must have at least three points") # ? or 2?
  if ( length(t) != length(y) )
    stop("Lengths of t and y in spline.data must be equal")
  if ( any(is.na(t)) || any(is.na(y)) )
    stop("Neither t nor y may contain NA values")
  if ( any(duplicated(t)) )
    stop("t cannot contain duplicated values")
  t.range <- obj$t.range
  if ( min(t) > t.range[1] || max(t) < t.range[2] )
    stop(sprintf("Spline data must span time range: [%s, %s]",
                 min(t.range), max(t.range)))
  
  deriv <- if (is.null(spline.data$deriv))
    0L else check.integer(spline.data$deriv)
  if ( length(deriv) != 1 )
    stop("spline.data$deriv must be an scalar")

  ## Renormalise y data onto [0,1], after projecting what the actual
  ## minimum and maximum values are.
  ## 
  ## TODO: This could be prone to error, especially as I've just
  ## hardcoded a "large" number of x spacing in here.  Ideally we
  ## would be able to provide this with the spline.data?  But the
  ## spline range may slightly exceed the values that are predicted
  ## analytically.
  tt <- seq(t.range[1], t.range[2], length.out=10001)
  r <- range(spline(t, y, xout=tt)$y)
  y <- (y - r[1]) / (r[2] - r[1])
  
  list(t=as.numeric(t), y=as.numeric(y), deriv=deriv)
}

## Information about possible types.
time.machine.types <- function() {
  info.t <- list(constant.t="c",
                 linear.t=c("c", "m"),
                 stepf.t=c("y0", "y1", "tc"),
                 sigmoid.t=c("y0", "y1", "tmid", "r"),
                 spline.t=c("y0", "y1"),
                 spline.linear.t=c("y0", "y1", "m"))
  info.t[.Call("r_get_time_machine_types", PACKAGE="diversitree")]
}

check.time.machine.functions <- function(functions) {
  info.t <- time.machine.types()
  if ( !is.character(functions) )
    stop("'functions' must be a character vector")
  if ( is.null(names(functions)) )
    stop("Functions must be named")
  info.t <- time.machine.types()
  types <- match(functions, names(info.t))
  if ( any(is.na(types)) )
    stop("Unknown time-varying function")
  types
}
