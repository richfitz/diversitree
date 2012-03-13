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
  f <- f # force argument to prevent recursion (pass by value)
  .x.eval <- list()
  .y.eval <- list()
  if ( interval < 0 )
    stop("Interval must be >= 0")
  function(x, ...) {
    i <- length(.x.eval) + 1
    if ( interval > 0 && i %% interval == 0 )
      cat(sprintf("[%s]", paste(formatC(x, 5), collapse=", ")))
    else if (interval > 0 )
      cat(".")
    .x.eval[[i]] <<- x
    .y.eval[[i]] <<- ans <- f(x, ...)
    if ( interval > 0 && i %% interval == 0 )
      cat(sprintf("\t -> %2.5f\n", ans))
    ans
  }
}

count.eval <- function(f) {
  n <- 0
  function(...) {
    n <<- n + 1
    f(...)
  }
}

## Known to work on 1.3-1.6
make.ode <- function(func, dllname, initfunc, ny, safe=FALSE) {
  if (!is.character(func)) 
    stop("`func' must be a character vector")
  if (!is.character(dllname))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")

  if ( safe ) {
    function(y, times, parms, rtol, atol) {
      ret <- t(lsoda(y, times, parms, rtol=rtol, atol=atol,
                     initfunc=initfunc, func=func, dllname=dllname))
      dimnames(ret) <- NULL
      ret
    }
  } else {
    initfunc.addr <- getNativeSymbolInfo(initfunc, PACKAGE=dllname)$address
    derivfunc.addr <- getNativeSymbolInfo(func, PACKAGE=dllname)$address
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

    rwork <- vector("double", 20)
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

    sol <- function(y, times, parms, rtol, atol) {
      if (!is.numeric(y)) 
        stop("`y' must be numeric")
      if (!is.numeric(times)) 
        stop("`times' must be numeric")
      storage.mode(y) <- storage.mode(times) <- "double"

      ret <- 
        .Call("call_lsoda", y, times, derivfunc.addr, parms,
              rtol, atol,
              NULL,      # rho: environment
              NULL,      # tcrit: critical times
              jacfunc, 
              initfunc.addr,
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


    ## Temporary fix so that I can work on the cluster.  This will be
    ## removed and DESCRIPTION updated to require R 2.12.0 or
    ## greater.
    if ( getRversion() >= "2.12.0" ) {
      vers <-  packageVersion("deSolve")
    } else {
      vers <- package_version(packageDescription("deSolve",
                                                 fields="Version"))
    }
    ## Update here when deSolve is updated.
    max.deSolve <- package_version("1.10-3")
    if ( vers > max.deSolve ) {
      str <- paste("diversitree is not known to work with deSolve > ",
                   max.deSolve, "\n\tfalling back on slow version")
      warning(str)
      make.ode(func, dllname, initfunc, ny, safe=TRUE)
    } else {
      sol
    }
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

## This is still experimental, and will not work nicely unless
## everything is nicely paired (it will not work well with constrained
## models, for example).
make.prior.ExpBeta <- function(r, beta) {
  to.pars2 <- function(pars) {
    m <- matrix(pars, 2)
    pars.mean <- colMeans(m)
    d <- 1 - (m[1,] / (2*pars.mean))
    rbind(pars.mean, d)
  }
  function(pars) {
    pars2 <- to.pars2(pars)
    sum(dexp(pars2[1,], r, log=TRUE)) +
      sum(dbeta(pars2[2,], beta, beta, log=TRUE))
  }
}

## TODO: Allow vector lower and upper here...
make.prior.uniform <- function(lower, upper) {
  if ( length(lower) == 2 && missing(upper) ) {
    upper <- lower[2]
    lower <- lower[1]
  }
  n <- length(lower)
  if ( length(upper) != n )
    stop("'lower' and 'upper' both be the same length")
  p <- log(1/(upper - lower))
  function(x) {
    ret <- rep(p, length=length(x))
    ret[x < lower | x > upper] <- 0
    ret
  }
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

## This only works for ultrametric trees:
branching.depth <- function(len, children, order, tips) {
  depth <- numeric(nrow(children))
  depth[tips] <- 0
  for ( i in order )
    depth[i] <- depth[children[i,1]] + len[children[i,1]]
  depth
}

## Convert a matrix to a list by row.
matrix.to.list <- function(m) {
  n <- nrow(m)
  out <- vector("list", n)
  for ( i in seq_len(n) )
    out[[i]] <- m[i,]
  out
}

matrix.to.list <- function(m) {
  storage.mode(m) <- "double"
  .Call("matrix_to_list", m)
}

argnames.twopart <- function(x, base, n.level) {
  obj <- attr(x, "argnames")
  if ( is.null(obj) )
    obj <- list(base=base, levels=seq_len(n.level))

  paste(obj$base, rep(obj$levels, each=length(obj$base)), sep=".")
}

argnames.twopart.set <- function(x, value, n.base, n.level) {
  if ( !is.list(value) || length(value) != 2 )
    stop("'value' must be a list of length 2")
  if ( length(value[[1]]) != n.base || length(value[[2]] != n.level) )
    stop(sprintf("value's elements must be of length %d, %d",
                 n.base, n.level))

  names(value) <- c("base", "levels")
  attr(x, "argnames") <- value
  x
}

set.defaults <- function(f, defaults) {
  if ( is.null(defaults) )
    return(f)
  if ( !all(names(defaults) %in% names(formals(f))) )
    stop("Unknown defaults")
  formals(f)[names(defaults)] <- defaults
  f
}
