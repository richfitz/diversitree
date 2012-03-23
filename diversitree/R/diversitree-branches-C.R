make.all.branches.C <- function(cache, model, dll, neq, np, comp.idx,
                                control) {
  check.cvodes(error=TRUE)

  tol <- control$tol
  eps <- control$eps

  cache <- toC.cache(cache, comp.idx)
  neq <- as.integer(neq)
  np  <- as.integer(np)

  rhs.name <- sprintf("derivs_%s_cvode", model)
  ic.name  <- sprintf("initial_conditions_%s", model)

  rhs <- getNativeSymbolInfo(rhs.name, PACKAGE=dll)$address
  ic  <- getNativeSymbolInfo(ic.name,  PACKAGE=dll)$address

  rtol <- as.numeric(tol)
  atol <- rep(as.numeric(tol), neq)

  ptr <- .Call("r_make_dt_obj", cache, neq, np, rhs, ic, rtol, atol,
               eps, PACKAGE="diversitree")

  function(pars, intermediates=FALSE, preset=NULL) {
    if ( !is.null(preset) )
      stop("Don't know how to deal with preset values yet")
    res <- .Call("r_all_branches", ptr, pars, PACKAGE="diversitree")
    names(res) <- c("lq", "vals")
    if ( intermediates ) {
      vals <- res$vals
      res <- .Call("r_get_vals", ptr, PACKAGE="diversitree")
      names(res) <- c("init", "base", "lq")
      res$vals <- vals
    }
    res
  }
}

## For all the xxSSE models, this converts and tidies the cache ready
## to be used in the C code.
toC.cache <- function(cache, comp.idx) {
  ## Translate tips...
  if ( is.null(cache$y) || length(cache$y) == 0 )
    stop("Can't do tipless yet...")
  f <- function(x)
    list(tip.y     = as.numeric(x$y),
         tip.len   = x$t,
         tip.target= toC.int(x$target))

  cache$y <- lapply(cache$y, f)

  cache$children <- toC.int(t(cache$children))
  cache$parent   <- toC.int(cache$parent)
  cache$order    <- toC.int(cache$order)
  cache$root     <- toC.int(cache$root)
  cache$n.tip    <- as.integer(cache$n.tip)
  cache$tips     <- toC.int(cache$tips)
  cache$comp.idx <- toC.int(comp.idx)

  cache
}

