## A special function will still be needed for cvodes and CVODES, as
## the normal branches function possibly won't work?
make.all.branches.t2.dtlik <- function(cache, control,
                                       initial.conditions) {
  control <- check.control.ode(control)
  time.machine <- cache$time.machine
  tm.ptr <- time.machine$ptr

  if ( control$backend == "CVODES" ) {
    setfunc <- sprintf("r_set_tm_%s", cache$info$name.ode)
    dummy <- rep(0.0, cache$info$np)
    all.branches <- make.all.branches.C(cache, control)
    function(pars, intermediates, preset=NULL) {
      time.machine$set(pars)
      .Call(setfunc, tm.ptr)
      all.branches(dummy, intermediates, preset)
    }
  } else {
    branches <- make.branches.dtlik(cache$info, control)
    initial.conditions.t <-
      make.initial.conditions.t2(cache, initial.conditions)

    if ( control$backend == "deSolve" ) {
      function(pars, intermediates, preset=NULL) {
        time.machine$set(pars)
        all.branches.matrix(tm.ptr, cache, initial.conditions.t,
                            branches, preset)
      }
    } else { # cvodes
      setfunc <- sprintf("r_set_tm_%s", cache$info$name.ode)
      dummy <- rep(0.0, cache$info$np)
      function(pars, intermediates, preset=NULL) {
        time.machine$set(pars)
        .Call(setfunc, tm.ptr)
        all.branches.matrix(dummy, cache, initial.conditions.t,
                            branches, preset)
      }
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

## Hard coded for bd right now, but will generalise really soon.
predict.bd.t2 <- function(object, p, t, nt=101, v=NULL,
                          alpha=1/20, everything=FALSE, ...) {
  tm <- get.cache(object)$time.machine
  if ( inherits(p, "fit.mle") || inherits(p, "mcmcsamples") )
    p <- coef(p)
  if ( missing(t) )
    t <- seq(tm$t.range[1], tm$t.range[2], length=nt)
  if ( is.null(v) )
    v <- tm$funnames
  is.matrix <- !is.null(dim(p)) && nrow(p) > 1

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
    tm$set(p)
    ret <- tm$getv(t)[,v,drop=FALSE]
  }
  ret <- list(t=t, y=ret)
}

plot.bd.t2 <- function(x, p, xlab="Time", ylab="Parameter",
                       lty=1, lwd=1, v=NULL, col=NULL,
                       nt=101, legend.pos=NULL, ...) {
  xy <- predict(x, p, v=v, nt=nt)
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
    legend(legend.pos, v, col=col, lty=lty, lwd=lwd)
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
