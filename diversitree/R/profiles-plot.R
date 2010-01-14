## Some simple MCMC plotting assistance.
## TODO: This needs to be fixed to allow for partial bars to be shaded.
add.profile.shading <- function(h, ci, col) {
  dx <- diff(h$mids[1:2])
  i <- which(with(h, mids > ci[1] & mids < ci[2]))
  if ( length(i) )
    with(h, polygon(rep(c(mids[c(i, i[length(i)]+1)] - dx/2), each=2),
                    c(0, rep(density[i], each=2), 0),
                    col=col, border=NA))
  else {
    ## The density is really constrained - just colour in the peak?
    j <- which.max(h$density)
    with(h, polygon(rep(c(mids[c(j, j+1)] - dx/2), each=2),
                    c(0, rep(density[j], each=2), 0),
                    col=col, border=NA))
  }
}
add.profile.outline <- function(h, col, vertical=FALSE) {
  dx <- diff(h$mids[1:2])
  if ( vertical )
    lines(h, freq=FALSE, col=col)
  else
    with(h, lines(mids - dx/2, density, type="s", col=col))
}
hdr.uniroot <- function(z, p=0.95) {
  xx <- c(0, seq(min(z), max(z), length=1024))
  ez <- ecdf(z)
  f <- suppressWarnings(approxfun(ez(xx), xx))
  ci <- optimize(function(x) f(x + p) - f(x), c(0, 1-p))$min
  f(c(ci, ci+p))
}
profiles.plot <- function(y, col.line, col.fill, xlim=NULL, n.br=50,
                          ...) {
  if ( missing(col.fill) )
    col.fill <- add.alpha(col.line, .5)
  r <- range(unlist(y))
  br <- seq(r[1], r[2], length=n.br)
  hh <- lapply(y, hist, br, plot=FALSE)
  ci <- lapply(y, hdr.uniroot)

  if ( is.null(xlim) ) xlim <- r
  ylim <- c(-0.075, 1.05) * max(sapply(hh, function(x) max(x$density)))

  plot(NA, xlim=xlim, ylim=ylim, type="n", yaxs="i", ...)
  for ( i in seq_along(y) )
    add.profile.shading(hh[[i]], ci[[i]], col.fill[i])
  for ( i in seq_along(y) )
    add.profile.outline(hh[[i]], col.line[i])

  z <- seq(0, 1, length=length(y) + 2)[-1] * par("usr")[3]
  for ( i in seq_along(y) )
    arrows(ci[[i]][1], z[i], ci[[i]][2], z[i], code=3, angle=90,
           length=0.02, col=col.line[i])
}

add.alpha <- function(col, alpha=.5) {
  tmp <- col2rgb(col)/255
  rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
}
