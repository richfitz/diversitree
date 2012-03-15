## Functions that provide useful lambda and mu functions.
sigmoid.x <- function(x, y0, y1, xmid, r)
  y0 + (y1 - y0)/(1 + exp(r * (xmid - x)))
sigmoid2.x <- function(x, y0, y1, xmid, r) 
  y0 + (y1 - y0)/(1 + exp(4 * r * (xmid - x) / (y1 - y0)))
constant.x <- function(x, c) rep(c, length(x))
noroptimal.x <- function(x, y0, y1, xmid, s2)
  y0 + (y1-y0)*exp(-(x - xmid)^2/(2 * s2))
make.linear.x <- function(x0, x1) {
  if ( is.null(x1) ) {
    function(x, c, m) {
      x1 <- length(x) - x0 + 1
      x[seq_len(x0)]  <- x[x0]
      x[x1:length(x)] <- x[x1]
      ans <- m * x + c
      ans[ans < 0] <- 0
      ans
    }
  } else {
    function(x, c, m) {
      x[x < x0] <- x0
      x[x > x1] <- x1
      ans <- m * x + c
      ans[ans < 0] <- 0
      ans
    }
  }
}
stepf.x <- function(x, y0, y1, xmid)
  ifelse(x < xmid, y0, y1)

normalise <- function(x) x / sum(x)

starting.point.quasse <- function(tree, states, states.sd=NULL)
  c(starting.point.bd(tree),
    diffusion=as.numeric(coef(find.mle(make.bm(tree, states,
      states.sd), .1))))

load.wisdom <- function(file="wisdom") {
  w <- paste(readLines(file), collapse="\n")
  .Call("r_set_wisdom", w, PACKAGE="diversitree")
}

save.wisdom <- function(file="wisdom") {
  w <- .Call("r_get_wisdom", PACKAGE="diversitree")
  write(w, file)
}

## Checking and sanitisation code:
check.f.quasse <- function(f) {
  args <- names(formals(f))
  if ( args[1] != "x" )
    stop("First argument of speciation/extinction function must be x")
  length(args) - 1
}

check.states.quasse <- function(tree, states, states.sd) {
  states <- check.states(tree, states)

  if ( length(states.sd) == 1 )
    states.sd <- structure(rep(states.sd, length(states)),
                           names=names(states))
  else
    states.sd <- check.states(tree, states.sd)
  
  list(states=states, states.sd=states.sd)
}  

check.control.quasse <- function(control, tree, states) {
  tree.length <- max(branching.times(tree))
  xr <- diff(range(states))
  xr.mult <- if ( "xr.mult" %in% names(control) )
    control$xr.mult else 5
  defaults <- list(tc=tree.length/10,
                   dt.max=tree.length/1000,
                   nx=1024,
                   dx=xr * xr.mult / 1024,
                   r=4,
                   xmid=mean(range(states)),
                   w=5,
                   method="fftC",
                   tips.combined=FALSE,
                   flags=FFTW.MEASURE, # fftC only
                   atol=1e-6, # mol only
                   rtol=1e-6, # mol only
                   eps=1e-6,  # perhaps scale with dx?
                   verbose=FALSE)

  nx.changed <- "nx" %in% names(control)
  dx.changed <- "dx" %in% names(control)
  control <- if ( is.null(control) )
    defaults else modifyList(defaults, control)
  if ( dx.changed && !nx.changed )
    control$nx <- 2^ceiling(log2(xr * xr.mult / control$dx))
  else if ( nx.changed && !dx.changed )
    control$dx <- xr * xr.mult / control$nx

  ## Eventually, this will contain "mol"
  method <- match.arg(control$method, c("fftC", "fftR", "mol"))

  if ( control$tips.combined && method != "fftC" )
    stop("'tips.combined' can only be used with method 'fftC'")

  if ( control$tc <= 0 || control$tc >= tree.length )
    stop(sprintf("tc must lie in (0, %2.2f)", tree.length))
  if ( log2(control$nx) %% 1 != 0 )
    stop("nx must be a power of two")
  if ( log2(control$r) %% 1 != 0 )
    stop("r must be a power of two")

  rr <- with(control, xmid + c(-1,1) * dx * nx / 2)
  rmin <- min(c(1, -1) * (mean(range(states)) - rr) / (xr / 2))
  if ( rmin - xr.mult < -1e-5 )
    warning("Range does not look wide enough - be careful!")
  else if ( rmin < 2 )
    stop("Range is not wide enough")

  ## These will be passed through to some C code, so type safety is
  ## important.
  ctrl.int <- c("nx", "flags", "verbose")
  ctrl.num <- c("tc", "dt.max", "r", "xmid", "w", "atol", "rtol")
  control[ctrl.int] <- sapply(control[ctrl.int], as.integer)
  control[ctrl.num] <- sapply(control[ctrl.num], as.numeric)

  control
}

## I use a number of elements of pars.
## pars[[i]]${lambda,mu,drift,diffusion,padding}
## pars$tr
expand.pars.quasse <- function(lambda, mu, args, ext, pars) {
  pars.use <- vector("list", 2)
  for ( i in c(1,2) ) {
    x <- list()
    pars.use[[i]] <-
      list(x=ext$x[i],
           lambda=do.call(lambda, c(ext$x[i], pars[args$lambda])),
           mu=do.call(mu, c(ext$x[i], pars[args$mu])),
           drift=pars[args$drift],
           diffusion=pars[args$diffusion],
           padding=ext$padding[i,])
  }
  names(pars.use) <- c("hi", "lo")
  pars.use$tr <- ext$tr
  pars.use
}

quasse.extent <- function(control, drift, diffusion) {
  nx <- control$nx
  dx <- control$dx
  dt <- control$dt.max
  xmid <- control$xmid
  r  <- control$r
  w <- control$w

  if ( control$method == "mol" ) {
    ndat <- nx*c(r, 1)
    padding <- NULL
  } else {
    mean <- drift * dt
    sd   <- sqrt(diffusion * dt)

    ## Another option here is to compute all the possible x values and
    ## then just drop the ones that are uninteresting?
    nkl <- max(ceiling(-(mean - w * sd)/dx)) * c(r, 1)
    nkr <- max(ceiling( (mean + w * sd)/dx)) * c(r, 1)
    ndat <- nx*c(r, 1) - (nkl + 1 + nkr)

    padding <- cbind(nkl, nkr)
    storage.mode(padding) <- "integer"
  }

  x0.2 <- xmid - dx*ceiling((ndat[2] - 1)/2)
  x0.1 <- x0.2 - dx*(1 - 1/r)

  ## Concatenate the x values, so that the lambda(x), mu(x)
  ## calculations work for both spaces simultaneously.
  x <- list(seq(x0.1, length.out=ndat[1], by=dx/r),
            seq(x0.2, length.out=ndat[2], by=dx))

  tr <- seq(r, length.out=ndat[2], by=r)

  list(x=x, padding=padding, ndat=ndat, tr=tr, nx=c(nx*r, nx))
}

