## There is enough weird stuff in the TGP optimiser that I have pulled
## it into its own file so that it does not pollute the relative order
## in mle.R

## This one behaves a little differently to the others.
##
## The tgp package provides tools for statistical optimisation; it
## will require substantially more tweaking than normal to make this
## work well.  An iterative strategy may be required, and it may not
## be possible to run this totally unsupervised.
do.mle.search.tgp <- function(func, x.init, control, lower, upper) {
  if ( !require(tgp) )
    stop("This method requires the tgp package")
  make.f.tgp <- function(func, parallel=FALSE) {
    if ( parallel && require(multicore) )
      lapply <- mclapply
    ## TODO: This should be OK when x is not a matrix, but currently
    ## fails.
    function(x)
      -as.numeric(unlist(lapply(matrix.to.list(x),
                                function(x) func(x))))
  }

  prev <- control$prev
  if ( is.null(prev) ) {
    ## Brand new
    if ( is.vector(x.init) )
      x.init <- matrix(x.init, 1, length(x.init),
                       dimnames=list(NULL, names(x.init)))
    z.init <- NULL
  } else {
    ## Pick up where we left off.
    if ( prev$method != "tgp" )
      stop("Only tgp fits can (currently) be continued")

    if ( missing(lower) ) lower <- prev$aux$rect[,1]
    if ( missing(upper) ) upper <- prev$aux$rect[,2]

    if ( missing(x.init) )
      x.init <- prev$aux$x
    else
      x.init <- rbind(prev$aux$x, x.init)
    z.init <- prev$aux$z
    npar <- ncol(x.init)
  }

  npar <- ncol(x.init)

  control <- modifyList(list(tgp.method=btgp,
                             improv.scope=2,
                             improv.n=20,
                             BTE=c(1, 1000, 1),
                             verbose=FALSE,
                             iterations=5,
                             n.init=10,
                             np=40 * npar,
                             parallel=FALSE,
                             run.optim=TRUE),
                        control)

  func2 <- make.f.tgp(func, control$parallel)
  lower <- check.par.length(lower, npar)
  upper <- check.par.length(upper, npar)
  rect <- cbind(lower, upper)
  if ( any(!is.finite(rect)) )
    stop("Finite bounds must be given")

  if ( nrow(x.init) < control$n.init )
    x <- rbind(x.init, lhs(control$n.init - nrow(x.init), rect))
  else
    x <- x.init

  ## Evaluate initial starting positions
  if ( is.null(z.init) )
    z <- func2(x)
  else if ( nrow(x.init) > length(z.init) )
    z <- c(z.init, func2(x.init[-seq_along(z.init),,drop=FALSE]))
  else
    z <- z.init

  improv <- c(control$improv.scope, control$improv.n)
  verbose <- control$verbose

  progress <- matrix(NA, control$iterations, npar + 1)
  colnames(progress) <- c(colnames(x.init), "p")

  for ( i in seq_len(control$iterations) ) {
    cat(sprintf("Iteration %d ... ", i))
    ok <- is.finite(z)
    x <- x[ok,,drop=FALSE]
    z <- z[ok]
    
    xx <- lhs(control$np, rect)
    fit.tgp <- control$tgp.method(X=x, Z=z, XX=xx, improv=improv,
                                  BTE=control$BTE, verb=verbose)

    x.new <- xx[fit.tgp$improv[,2] <= improv[2],]

    if ( control$run.optim ) {
      x.all <- rbind(x, xx)
      x.best <- x.all[which.min(c(fit.tgp$Zp.mean, fit.tgp$ZZ.mean)),]
      
      opt <- optim.tgp(fit.tgp, x.best, rect)

      ## Here are the new suggested x values; those from the
      ## improvement, plus the one run against the predicted surface.
      x.new <- rbind(opt$par, x.new)
    }

    x <- rbind(x, x.new)
    z <- c(z, func2(x.new))

    j <- which.min(z)
    progress[i,] <- c(x[j,], -z[j])
    cat(sprintf("best point: %2.5f\n", -z[j]))
  }

  if ( !is.null(prev) )
    progress <- rbind(prev$aux$progress, progress)

  list(par=progress[i,seq_len(npar)],
       lnLik=progress[i,npar+1],
       counts=length(z),
       iterations=control$iterations,
       convergence=NA,
       message="",
       hessian=NULL,
       aux=list(x=x, z=z, progress=progress, fit=fit.tgp, rect=rect))
}

## Rewritten version of the function to do optimisation on a tgp
## object.
optim.tgp <- function(obj, start, rect, dx=1e-5, ...) {
  npar <- nrow(rect)
  f <- function(x) {
    x <- matrix(x, ncol=npar)
    as.vector(predict(obj, XX=x, pred.n=FALSE)$ZZ.km)
  }
  dx <- check.par.length(dx, npar)
  gr <- make.grad.fd(f, dx, rect[,1], rect[,2])

  optim(start, f, gr, method="L-BFGS-B",
        lower=rect[,1], upper=rect[,2], ...)  
}

## Finite difference gradient calculation, using the same algorithm as
## optim() uses internally.
##
## Why is this useful?  Because the list of points can be generated
## ahead of time, if the cost of evaluating many points at once is not
## much greater than the cost of evaluating a single point.  I have
## done this just for the bounded version, but removing these should
## be easy.
make.grad.fd <- function(f, dx, lower, upper, ...) {
  func <- function(x) f(x, ...)
  npar <- length(dx)
  function(x) {
    x0 <- x - dx
    x1 <- x + dx

    shrink.0 <- x0 < lower
    shrink.1 <- x1 > upper

    x0[shrink.0] <- lower[shrink.0]
    x1[shrink.1] <- upper[shrink.1]

    eps <- (dx * !shrink.0) + (dx * !shrink.1)

    ## build the list of coordinates:
    g <- function(x, i, v) {x[i] <- v; x}
    xx <- c(lapply(seq_len(npar), function(i) g(x, i, x0[i])),
            lapply(seq_len(npar), function(i) g(x, i, x1[i])))
    xx <- matrix(unlist(xx), 2*npar, npar, TRUE)

    ## Evaluate and compute the FD approximation to the gradient.
    tmp <- matrix(f(xx), npar, 2)
    (tmp[,2] - tmp[,1]) / eps
  }
}


