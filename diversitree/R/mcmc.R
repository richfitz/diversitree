## Slice sampling

## These functions take a function 'f' that evaluates an argument 'x'
## which is a vector of parameters.  Across a single iteration,
## take the input x and return a vector of parameters 'y'
## corresponding to a new position.
mcmc <- function(f, x0, nsteps, w, lower=-Inf, upper=Inf,
                 fail.value=-Inf, print.every=1, ...) {
  fn <- protect(function(x) f(x, ...), fail.value)

  y0 <- try(f(x0, ...))
  if ( inherits(y0, "try-error") )
    stop("The above error occured when testing the starting position")
  if ( !is.finite(y0) )
    stop("Starting point must have finite probability")

  if ( length(lower) == 1 ) lower <- rep(lower, length(x0))
  if ( length(upper) == 1 ) upper <- rep(upper, length(x0))

  if ( is.null(names(x0)) )
    try(names(x0) <- argnames(f), silent=TRUE)

  hist <- vector("list", nsteps)
  for ( i in seq_len(nsteps) ) {
    hist[[i]] <- tmp <- slice.nd(fn, x0, y0, w, lower, upper)
    x0 <- tmp[[1]]
    y0 <- tmp[[2]]
    if ( print.every > 0 && i %% print.every == 0 )
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
slice.nd <- function(f, x0, y0, w, lower, upper) {
  if ( is.na(y0) )
    y0 <- f(x0)
  for ( i in seq_along(x0) ) {
    xy <- slice.1d(make.g(f, x0, i), x0[i], y0, w[i], lower[i],
                   upper[i])
    x0[i] <- xy[1]
    y0 <- xy[2]
  }
  list(x0, y0)
}

## Here, w, lower and upper are scalars
slice.1d <- function(g, x0, y0, w, lower, upper) {
  z <- y0 - rexp(1)
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
