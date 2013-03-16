check.control.continuous <- function(control) {
  defaults <- list(method="vcv", backend="R")
  control <- modifyList(defaults, control)

  if ( length(control$method) != 1 )
    stop("control$method must be a scalar")
  if ( identical(control$method, "direct") ) {
    warning('method="direct" is deprecated, please use method="pruning"')
    control$method <- "pruning"
  }
    
  methods <- c("vcv", "pruning")
  if ( !(control$method %in% methods) )
    stop(sprintf("control$method must be in %s",
                 paste(methods, collapse=", ")))

  if ( length(control$backend) != 1 )
    stop("control$backend must be a scalar")
  backends <- c("C", "R")
  if ( !(control$backend %in% backends) )
    stop(sprintf("control$backend must be in %s",
                 paste(backends, collapse=", ")))
  
  control
}

make.all.branches.continuous <- function(cache, control) {
  cache.C <- toC.cache.continuous(cache)

  br.name <- sprintf("branches_%s", cache$info$name)
  br <- getNativeSymbolInfo(br.name)$address
  ## Hard coded now, as bm and ou share.
  ic <- getNativeSymbolInfo("initial_conditions_bm")$address
  
  ptr <- .Call("r_make_dt_obj_cont", cache.C, ic, br,
               PACKAGE="diversitree")
  
  function(pars, intermediates=FALSE, preset=NULL) {
    if ( !is.null(preset) )
      stop("Don't know how to deal with preset values yet")
    res <- .Call("r_all_branches_cont", ptr, pars,
                 PACKAGE="diversitree")
    names(res) <- c("lq", "vals")
    if ( intermediates ) {
      vals <- res$vals
      res <- .Call("r_get_vals_cont", ptr, PACKAGE="diversitree")
      names(res) <- c("init", "base", "lq")
      res$vals <- vals
    }
    res    
  }
}

toC.cache.continuous <- function(cache) {
  ## This is super awkward for now, but may change in future.
  ## Basically, for this version, I don't want the indices ordered at
  ## all, and assert that the order should be 1..ntips.
  if ( length(cache$y$target) != cache$n.tip )
    stop("Missing tips")
  i <- order(cache$y$target)
  cache$y$target <- toC.int(cache$y$target[i]) # now seq_len(n.tip)-1L
  cache$y$y      <- cache$y$y[,i,drop=FALSE]
  cache$y$t      <- cache$y$t[i]
  cache$y$type   <- NULL

  cache$children <- toC.int(t(cache$children))
  cache$parent   <- toC.int(cache$parent)
  cache$order    <- toC.int(cache$order)
  cache$root     <- toC.int(cache$root)
  cache$n.tip    <- as.integer(cache$n.tip)
  cache$tips     <- toC.int(cache$tips)

  cache
}
