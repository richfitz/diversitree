## Slice sampling

## These functions take a function 'f' that evaluates an argument 'x'
## which is a vector of parameters.  Across a single iteration,
## take the input x and return a vector of parameters 'y'
## corresponding to a new position.
TYPE.RAW <- 1
TYPE.LOG <- 2
TYPE.NLOG <- 3

mcmc <- function(f, x0, y0, nsteps, w, lower, upper, type=TYPE.LOG,
                 ...) {
  fn <- function(x) f(x, ...)

  if ( missing(y0) )
    y0 <- fn(x0)

  hist <- vector("list", nsteps)
  for ( i in seq_len(nsteps) ) {
    hist[[i]] <- tmp <- slice.nd(fn, x0, y0, w, lower, upper, type)
    x0 <- tmp[[1]]
    y0 <- tmp[[2]]
    cat(sprintf("%d: {%s} -> %2.5f\n", i,
                paste(sprintf("%2.4f", tmp[[1]]), collapse=", "),
                tmp[[2]]))
  }

  hist <- cbind(i=seq_along(hist),
                as.data.frame(t(sapply(hist, unlist))))

  names(hist)[ncol(hist)] <- "p"
  hist
}

## Here, w, lower and upper are vectors
slice.nd <- function(f, x0, y0, w, lower, upper, type=TYPE.LOG) {
  if ( is.na(y0) )
    y0 <- f(x0)
  for ( i in seq_along(x0) ) {
    xy <- slice.1d(make.g(f, x0, i), x0[i], y0, w[i], lower[i],
                   upper[i], type)
    x0[i] <- xy[1]
    y0 <- xy[2]
  }
  list(x0, y0)
}

## Here, w, lower and upper are scalars
slice.1d <- function(g, x0, y0, w, lower=0, upper=Inf,
                     type=TYPE.LOG) {
  if ( type %in% c(TYPE.LOG, TYPE.NLOG) )
    z <- y0 - rexp(1)
  else
    z <- runif(1) * y0
  r <- isolate.step(g, x0, y0, z, w, lower, upper)
  take.sample(g, x0, z, r)
}

take.sample <- function(g, x0, z, r) {
  r0 <- r[1]
  r1 <- r[2]

  repeat {
    xs <- runif(1, r0, r1)
    ys <- g(xs)
    if ( ys > z )
      break
    if ( xs < x0 )
      r0 <- xs
    else
      r1 <- xs
  }
  c(xs, ys)
}

isolate.step <- function(g, x0, y0, z, w, lower, upper) {
  u <- runif(1) * w
  L <- x0 - u
  R <- x0 + (w-u)

  while ( L > lower && g(L) > z )
    L <- L - w
  while ( R < upper && g(R) > z )
    R <- R + w

  c(max(L, lower), min(R, upper))
}

make.g <- function(f, x, i) {
  function(z) {
    x[i] <- z
    f(x)
  }
}
