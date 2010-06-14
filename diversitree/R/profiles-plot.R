## TODO: I still do not deal with the case where there is only a
## single unique y value very nicely.  I should probably draw a delta
## function?
profiles.plot <- function(y, col.line, col.fill, xlim=NULL, ymax=NULL,
                           n.br=50, ...) {
  if ( missing(col.fill) )
    col.fill <- add.alpha(col.line, .5)
  if ( is.null(xlim) )
    r <- range(unlist(y))
  else
    r <- xlim
  br <- seq(r[1], r[2], length=n.br)
  hh <- lapply(y, hist, br, plot=FALSE)
  ci <- lapply(y, hdr)

  if ( is.null(xlim) ) xlim <- r
  if ( is.null(ymax) )
    ymax <- max(sapply(hh, function(x) max(x$density)))
  ylim <- c(-0.075, 1.05) * ymax
  
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

add.profile.shading <- function(h, ci, col) {
  dx <- diff(h$mids[1:2])
  xx <- c(h$mids - dx / 2, h$mids[length(h$mids)] + dx / 2)
  i <- which(xx > ci[1] & xx < ci[2])
  xs <- rep(c(ci[1], xx[i], ci[2]), each=2)
  j <- if ( length(i) > 1 ) min(i) - 1 else 1
  ys <- c(0, rep(h$density[c(j, i)], each=2), 0)
  polygon(xs, ys, col=col, border=NA)
}
add.profile.outline <- function(h, col, vertical=FALSE) {
  dx <- diff(h$mids[1:2])
  if ( vertical )
    lines(h, freq=FALSE, col=col)
  else {
    xx <- rep(with(h, c(mids-dx/2, mids[length(mids)]+dx/2)), each=2)
    yy <- c(0, rep(h$density, each=2), 0)
    lines(xx, yy, col=col)
  }
}

hdr.uniroot <- function(z, p=0.95, lim=NULL) {
  xx <- sort(c(lim, seq(min(z), max(z), length=1024)))
  ez <- ecdf(z)
  f <- suppressWarnings(approxfun(ez(xx), xx))
  fit <- suppressWarnings(optimize(function(x)
                                   f(x + p) - f(x), c(0, 1-p)))
  if ( inherits(fit, "try-error") || is.na(fit$objective) )
    stop("HDR interval failure")
  ci <- fit$min
  f(c(ci, ci+p))
}

hdr <- function(z, p=0.95, lim=NULL) {
  ci <- try(hdr.uniroot(z, p, lim), silent=TRUE)
  if ( inherits(ci, "try-error") ) {
    warning("HDR falling back on quantile-based intervals")
    ci <- as.numeric(quantile(z, c((1-p)/2, 1/2 + p/2)))
  }
  ci
}

## hdr.uniroot <- function(z, p=0.95, lim=NULL) {
##   xx <- sort(c(lim, seq(min(z), max(z), length=1024)))
##   ez <- ecdf(z)
##   f <- suppressWarnings(approxfun(ez(xx), xx))
##   fit <- try(suppressWarnings(optimize(function(x) f(x + p) - f(x), c(0, 1-p))))
##   if ( inherits(fit, "try-error") || is.na(fit$objective) ) {
##     warning("HDR falling back on quantile-based intervals")
##     ci <- as.numeric(quantile(z, c((1-p)/2, 1/2 + p/2)))
##   } else {
##     ci <- fit$min
##     f(c(ci, ci+p))
##   }
## }

add.alpha <- function(col, alpha=.5) {
  tmp <- col2rgb(col)/255
  rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
}
