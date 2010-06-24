protect <- function(f, fail.value.default=NULL) {
  function(..., fail.value=fail.value.default, finite=TRUE) {
    if ( is.null(fail.value) )
      f(...)
    else {
      ret <- try(f(...), silent = TRUE)
      failed <- (inherits(ret, "try-error") ||
                 (finite && (is.na(ret) || !is.finite(ret))))
      if ( failed )
        fail.value
      else
        ret
    }
  }
}

invert <- function(f) function(...) -f(...)

## Box constraints
boxconstrain <- function(f, lower, upper, fail.value=-Inf) {
  function(x, ...) {
    if ( any(x < lower | x > upper) )
      fail.value
    else
      f(x, ...)
  }
}

big.brother <- function(f, interval=1) {
  .x.eval <- list()
  .y.eval <- list()
  function(x, ...) {
    i <- length(.x.eval) + 1
    if ( i %% interval == 0 )
      cat(sprintf("[%s]", paste(formatC(x, 5), collapse=", ")))
    else
      cat(".")
    .x.eval[[i]] <<- x
    .y.eval[[i]] <<- ans <- f(x, ...)
    if ( i %% interval == 0 )
      cat(sprintf("\t -> %2.5f\n", ans))
    ans
  }
}

## Known to work on 1.3-1.6
## TODO: fall back on "safe" with a warning if the an unknown version.
make.ode <- function(func, dllname, initfunc, ny, safe=FALSE) {
  if (!is.character(func)) 
    stop("`func' must be a character vector")
  if (!is.character(dllname))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")

  if ( safe ) {
    function(y, times, parms, rtol, atol)
      t(lsoda(y, times, parms, rtol=rtol, atol=atol,
              initfunc=initfunc, func=func, dllname=dllname))
  } else {
    initfunc <- getNativeSymbolInfo(initfunc, PACKAGE=dllname)$address
    derivfunc <- getNativeSymbolInfo(func, PACKAGE=dllname)$address
    jacfunc <- NULL

    maxordn <- 12
    maxords <- 5
    lrw <- as.integer(max(20 + ny * (maxordn + 1) + 3 * ny,
                          20 + ny * (maxords + 1) + 3 * ny + ny^2 + 2))
    liw <- as.integer(20 + ny)

    iwork <- vector("integer", 20)
    iwork[1] <- as.integer(1)    # banddown
    iwork[2] <- as.integer(1)    # bandup
    iwork[6] <- as.integer(5000) # maxsteps

    rwork <- vector("double", 20) # TODO: 5, 6 do nothing.
    rwork[5] <- 0                     # hini
    rwork[6] <- 10                    # hmax (consider 0)
    rwork[7] <- 0                     # hmin

    INTZERO <- as.integer(0)
    INTONE <- as.integer(1)
    INTTWO <- as.integer(2)

    ## "Forcings" are new to deSolve version 1.5; this is the required
    ## argument where they are not used (I am not using them)
    flist <- list(fmat=0, tmat=0, imat=0, ModelForc=NULL)

    elist <- list()
    eventfunc <- NULL
    elag <- list(islag=0L)

    vers <- packageDescription("deSolve")$Version
    vers <- strsplit(vers, "-")[[c(1,1)]]

    f.1.5 <- function(y, times, parms, rtol, atol) {
      if (!is.numeric(y)) 
        stop("`y' must be numeric")
      if (!is.numeric(times)) 
        stop("`times' must be numeric")
      storage.mode(y) <- storage.mode(times) <- "double"

      ret <- 
        .Call("call_lsoda", y, times, derivfunc, parms, rtol, atol,
              NULL,      # rho: environment
              NULL,      # tcrit: critical times
              jacfunc, 
              initfunc,
              INTZERO,   # verbose (false)
              INTONE,    # itask
              rwork,
              iwork,
              INTTWO,    # jT: Jacobian type (fullint)
              INTZERO,   # nOut (no global variables)
              lrw,       # size of workspace (real)
              liw,       # size of workspace (int)
              INTONE,    # Solver
              NULL,      # rootfunc
              INTZERO,   # nRoot
              0,         # rpar: no extra real parameters
              INTZERO,   # ipar: no extra integer parameters
              INTZERO,   # Type
              flist,     # [New in 1.5]
              PACKAGE="deSolve")
      if ( max(abs(ret[1,] - times)) > 1e-6 )
        stop("Integration error: integration stopped prematurely")
      ret
    }
    f.1.6 <- function(y, times, parms, rtol, atol) {
      if (!is.numeric(y)) 
        stop("`y' must be numeric")
      if (!is.numeric(times)) 
        stop("`times' must be numeric")
      storage.mode(y) <- storage.mode(times) <- "double"

      ret <- 
        .Call("call_lsoda", y, times, derivfunc, parms, rtol, atol,
              NULL,      # rho: environment
              NULL,      # tcrit: critical times
              jacfunc, 
              initfunc,
              eventfunc, # [New in 1.6]
              INTZERO,   # verbose (false)
              INTONE,    # itask
              rwork,
              iwork,
              INTTWO,    # jT: Jacobian type (fullint)
              INTZERO,   # nOut (no global variables)
              lrw,       # size of workspace (real)
              liw,       # size of workspace (int)
              INTONE,    # Solver
              NULL,      # rootfunc
              INTZERO,   # nRoot
              0,         # rpar: no extra real parameters
              INTZERO,   # ipar: no extra integer parameters
              INTZERO,   # Type
              flist,     # [New in 1.5]
              elist,     # [New in 1.6]
              PACKAGE="deSolve")
      if ( max(abs(ret[1,] - times)) > 1e-6 )
        stop("Integration error: integration stopped prematurely")
      ret
    }

    f.1.7 <- function(y, times, parms, rtol, atol) {
      if (!is.numeric(y)) 
        stop("`y' must be numeric")
      if (!is.numeric(times)) 
        stop("`times' must be numeric")
      storage.mode(y) <- storage.mode(times) <- "double"

      ret <- 
        .Call("call_lsoda", y, times, derivfunc, parms, rtol, atol,
              NULL,      # rho: environment
              NULL,      # tcrit: critical times
              jacfunc, 
              initfunc,
              eventfunc, # [New in 1.6]
              INTZERO,   # verbose (false)
              INTONE,    # itask
              rwork,
              iwork,
              INTTWO,    # jT: Jacobian type (fullint)
              INTZERO,   # nOut (no global variables)
              lrw,       # size of workspace (real)
              liw,       # size of workspace (int)
              INTONE,    # Solver
              NULL,      # rootfunc
              INTZERO,   # nRoot
              0,         # rpar: no extra real parameters
              INTZERO,   # ipar: no extra integer parameters
              INTZERO,   # Type
              flist,     # [New in 1.5]
              elist,     # [New in 1.6]
              elag,      # [New in 1.7]
              PACKAGE="deSolve")
      if ( max(abs(ret[1,] - times)) > 1e-6 )
        stop("Integration error: integration stopped prematurely")
      ret
    }

    switch(vers,
           "1.5"=f.1.5, "1.5-1"=f.1.5, "1.6"=f.1.6, "1.7"=f.1.7,
           stop("Cannot use diversitree with deSolve version ", vers))
  }
}



quadratic.roots <- function(a, b, c)
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)


discretize <- function(x, n, r=range(x)) {
  at <- seq(r[1], r[2], length=n+1)
  as.integer(cut(x, at, include.lowest=TRUE, labels=FALSE))
}

make.prior.exponential <- function(r) {
  function(pars)
    sum(log(r) - pars * r)
}


## This is the same as the function in quasse2
descendants <- function(node, edge) {
  ans <- node
  repeat {
    node <- edge[edge[,1] %in% node,2]
    if ( length(node) > 0 )
      ans <- c(ans, node)
    else
      break
  }

  unlist(ans)
}

ancestors <- function(phy, i=seq_along(phy$tip.label)) {
  anc <- i
  edge <- phy$edge
  while ( any(!is.na(i)) ) {
    i <- edge[match(i, edge[,2]),1]
    anc <- cbind(anc, i, deparse.level=0)
  }

  apply(anc, 1, function(x)
        c(rev(x[!is.na(x)]), rep(NA, sum(is.na(x)))))
}

## Compute the MRCA of tips with indices in 'tips'
mrca.tipset <- function(phy, tips) {
  if ( is.character(tips) )
    tips <- match(tips, phy$tip.label)
  
  if ( length(tips) == 1 )
    tips
  else {
    anc <- ancestors(phy, tips)
    j <- which(apply(anc, 1, function(x) length(unique(x))) > 1)[1]
    anc[j-1,1]
  }
}

## Similar to ape's branching.times(), but returns the height above
## the root node, even for non-ultrametric trees.  Includes tip times.
branching.heights <- function(phy) {
  if (!inherits(phy, "phylo"))
    stop('object "phy" is not of class "phylo"')

  edge <- phy$edge
  n.node <- phy$Nnode
  n.tip <- length(phy$tip.label)

  ht <- numeric(n.node + n.tip) # zero'd
  for (i in seq_len(nrow(edge)) )
    ht[edge[i, 2]] <- ht[edge[i, 1]] + phy$edge.length[i]

  ## Ugly, but fairly compatible with branching.times()
  names.node <- phy$node.label
  if (is.null(names.node))
    names.node <- (n.tip + 1):(n.tip + n.node)
  names(ht) <- c(phy$tip.label, names.node)

  ht
}

## Convert a matrix to a list by row.
matrix.to.list <- function(m) {
  n <- nrow(m)
  out <- vector("list", n)
  for ( i in seq_len(n) )
    out[[i]] <- m[i,]
  out
}
