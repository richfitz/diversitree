## For all the xxSSE models, this converts and tidies the cache ready
## to be used in the C code.  The tip treatment is hardest, but I
## think that it is correct.
toC.int <- function(x) {
  x <- x - 1
  storage.mode(x) <- "integer"
  x
}

toC.cache <- function(cache, comp.idx) {
  ## Translate tips...
  if ( is.null(cache$y) || length(cache$y) == 0 )
    stop("Can't do tipless yet...")
  f <- function(x) {
    len <- x$t.uniq[x$unpack]
    list(tip.y = as.numeric(x$y),
         tip.len   = sort(len),
         tip.target= toC.int(x$target[order(len)]))
  }
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

  function(pars, intermediates=FALSE) {
    vals <- .Call("r_all_branches", ptr, pars, PACKAGE="diversitree")
    names(vals) <- c("lq", "vals")
    if ( intermediates ) {
      int <- .Call("r_get_vals", ptr, PACKAGE="diversitree")
      names(int) <- c("init", "base", "lq")
      vals$intermediates <- int
    }
    vals
  }
}

ll.xxsse.C <- function(pars, all.branches,
                       condition.surv=TRUE, root=ROOT.OBS,
                       root.p=NULL, intermediates=FALSE) {
  ans <- all.branches(pars, intermediates)
  vals <- ans[[2]]
  root.p <- root.p.xxsse(vals, pars, root, root.p)
  loglik <- root.xxsse(vals, pars, ans[[1]], condition.surv, root.p)

  if ( intermediates ) {
    ans$intermediates$root.p <- root.p
    attr(loglik, "intermediates") <- ans$intermediates
    attr(loglik, "vals") <- vals
  }

  loglik
}
