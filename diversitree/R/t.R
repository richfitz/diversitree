## A special function will still be needed for cvodes and CVODES, as
## the normal branches function possibly won't work?
make.all.branches.t.dtlik <- function(cache, control,
                                      initial.conditions.base) {
  control <- check.control.ode(control)

  if ( control$backend == "CVODES" ) {
    make.all.branches.C(cache, control)
  } else {
    branches <- make.branches.dtlik(cache$info, control)
    initial.conditions <-
      make.initial.conditions.t(cache, initial.conditions.base)
    function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache, initial.conditions,
                          branches, preset)
  }
}

update.cache.t <- function(cache, functions, spline.data, with.q=FALSE) {
  info <- cache$info
  n.args <- length(info$argnames)

  ## This could be check.functions.t():
  if ( !is.character(functions) )
    stop("'functions' must be characters [will relax soon]")
  if ( length(functions) == 1L )
    functions <- rep(functions, n.args)
  else if ( length(functions) != n.args )
    stop(sprintf("Expected %d functions", n.args))
  if ( is.null(names(functions)) )
    names(functions) <- info$argnames

  t.range <- range(0, cache$depth[cache$root])
  nonnegative <- cache$nonnegative
  k <- if ( with.q ) cache$info$k else NULL
  cache$time.machine <-
    make.time.machine(functions, t.range, nonnegative, spline.data, k)

  info$time.varying <- TRUE
  info$argnames <- cache$time.machine$argnames
  info$name.ode <- sprintf("%s_t", cache$info$name)
  info$name.pretty <- sprintf("%s (time-varying[v2])", info$name.pretty)
  info$name <- sprintf("%s.t", cache$info$name)

  cache$info <- info
  cache
}

make.initial.conditions.t <- function(cache, initial.conditions) {
  pars.t <- cache$time.machine$get
  function(init, pars, t, idx)
    initial.conditions(init, pars.t(t), t, idx)
}

## TODO/NEW This is somewhat tedious, as we really should check for
## 'root' not being ROOT.EQUI, as that can't be done in time-dependent
## models.  However, because this is used by different functions with
## different argument lists, that is hard to do.  But then, that check
## is duplicated in too many functions.  For now, I'm skipping this.
make.rootfunc.t <- function(cache, rootfunc) {
  pars.t <- cache$time.machine$get
  t.root <- cache$depth[cache$root]
  function(ans, pars, ...)             # pars here is ignored...
    rootfunc(ans, pars.t(t.root), ...) # ...because tm version used.
}

make.prep.all.branches.t <- function(cache, backend) {
  time.machine <- cache$time.machine
  tm.ptr <- time.machine$ptr
  setfunc <- sprintf("r_set_tm_%s", cache$info$name.ode)
  dummy <- rep(0.0, cache$info$np)
  attr(dummy, "dummy") <- TRUE # for identical() to work.

  function(pars) {
    if ( !identical(pars, dummy) )
      time.machine$set(pars)
    else
      message("Not setting parameters - must already be done?")
    .Call(setfunc, tm.ptr)
    dummy
  }
}

## Making the output useful.
predict.dtlik.t <- function(object, p, t, nt=101, v=NULL, thin=10,
                            alpha=1/20, everything=FALSE, ...) {
  tm <- get.cache(object)$time.machine
  if ( inherits(p, "fit.mle") || inherits(p, "mcmcsamples") )
    ## TODO: Improve the coef.mcmcsamples to allow full, here, then
    ## use full.  Possibly add a 'lik' function in that can do the
    ## resolution below?
    p <- coef(p)
  if ( missing(t) )
    t <- seq(tm$t.range[1], tm$t.range[2], length=nt)
  if ( is.null(v) )
    v <- tm$funnames
  is.matrix <- !is.null(dim(p)) && nrow(p) > 1
  ## Thin the chain to speed things up?
  if ( is.matrix && thin > 1 )
    p <- p[seq(1, nrow(p), by=thin),,drop=FALSE]
  if ( inherits(object, "constrained") ) {
    if ( is.matrix && ncol(p) == length(argnames(object)) )
      p <- t(apply(p, 1, object, pars.only=TRUE))
    else if ( !is.matrix && length(p) == length(argnames(object)) )
      p <- object(p, pars.only=TRUE)
  }

  if ( is.matrix ) {
    if ( everything ) {
      np <- if ( is.matrix ) nrow(p) else 1
      ret <- lapply(t, function(ti)
                    t(apply(p, 1, tm$get1, ti, i=match(v, tm$funnames))))
      ret <- array(unlist(ret), c(np, length(v), length(t)))
      tmp <- aperm(ret, c(2, 3, 1))
      dimnames(tmp) <- list(v, NULL, NULL)
      ret <- vector("list", length(v))
      for ( i in v )
        ret[[i]] <- tmp[i,,]
    } else { 
      ret <- lapply(match(v, tm$funnames), function(i)
                    average.over.mcmc(p, tm$get1, t, i=i))
      names(ret) <- v
    }
  } else {
    ret <- tm$getv(t, p)[,v,drop=FALSE]
  }
  list(t=t, y=ret)
}

plot.dtlik.t <- function(x, p, xlab="Time", ylab="Parameter",
                         lty=1, lwd=1, v=NULL, col=NULL,
                         nt=101, legend.pos=NULL, thin=10, ...) {
  xy <- predict(x, p, v=v, nt=nt, thin=thin)
  is.matrix <- !is.null(dim(p)) && nrow(p) > 1
  xlim <- rev(range(xy$t))  
  if ( is.matrix ) {
    v <- names(xy$y)
    if ( is.null(col) )
      col <- seq_along(v)
    fill <- add.alpha(col, .5)
    names(col) <- names(fill) <- v
    ylim <- range(lapply(xy$y, range))
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, las=1)
    tt <- c(xy$t, rev(xy$t))
    for ( i in rev(v) )
      polygon(tt, c(xy$y[[i]][,"lower"], rev(xy$y[[i]][,"upper"])),
              col=fill[i], border=NA)
    for ( i in rev(v) )
      lines(xy$t, xy$y[[i]][,"mean"], col=col[i])
  } else {
    v <- colnames(xy$y)
    if ( is.null(col) )
      col <- seq_along(v)
    matplot(xy$t, xy$y, type="l", xlim=xlim, las=1, lty=lty,
            col=col, xlab=xlab, ylab=ylab, lwd=lwd, ...)
  }
  if ( !is.null(legend.pos) )
    legend(legend.pos, v, col=col, lty=lty, lwd=lwd, bty="n")
  invisible(v)
}

## Given a function that takes arguments
##   f(p, x, ...)
## where p is the parameter vector and 'x' is the domain position.
## Additional arguments are passed in.
average.over.mcmc <- function(p, f, xx, ..., alpha=1/20) {
  g <- function(xi) {
    y <- apply(p, 1, f, xi, ...)
    c(mean(y), hdr.new(y, alpha))
  }
  if ( is.data.frame(p) )
    p <- as.matrix(p)
  ret <- t(sapply(xx, g))
  colnames(ret) <- c("mean", "lower", "upper")
  ret
}
