## make.all.branches.C <- function(cache, model, dll, neq, np, comp.idx,
##                                 control) {
make.all.branches.C <- function(cache, control) {
  check.cvodes(error=TRUE)
  info     <- check.info.ode(cache$info, control)
  model    <- info$name.ode
  neq      <- as.integer(info$ny)
  np       <- as.integer(info$np)
  comp.idx <- as.integer(info$idx.d)
  dll      <- info$dll
  rtol <- as.numeric(control$tol)
  atol <- rep(rtol, neq)
  eps  <- as.numeric(control$eps)

  ## This is depended on by make.do.asr.marginal, which reaches in
  ## here to access the cache.
  cache.C <- toC.cache(cache, comp.idx)

  rhs.name <- sprintf("derivs_%s_cvode", model)
  ic.name  <- sprintf("initial_conditions_%s", model)

  rhs <- getNativeSymbolInfo(rhs.name, PACKAGE=dll)$address
  ic  <- getNativeSymbolInfo(ic.name,  PACKAGE=dll)$address

  ptr <- .Call("r_make_dt_obj", cache.C, neq, np, rhs, ic,
               rtol, atol, eps, PACKAGE="diversitree")

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
