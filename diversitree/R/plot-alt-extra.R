## TODO: phylogram version.

## obj: result of plot2.phylo
## lab: vector along the tip labels with a "group" level.  There will
## be repetition.
## col.bar and col.lab: vectors along sort(unique(lab)) with colours
## for the bar and label.
group.label.tip.rad <- function(obj, lab, col.bar, col.lab, lwd=1,
                                offset.bar=0, offset.lab=0, cex=1,
                                font=1, ...) {
  n.taxa <- obj$n.taxa
  n <- obj$n.spp
  if ( is.null(n.taxa) )
    dt <- 1/6 / n * 2 * pi
  else
    dt <- (n.taxa/2 - .5 + 1/6) / n * 2 * pi

  theta <- obj$xy$theta[seq_len(obj$Ntip)]

  t0 <- tapply(theta - dt, lab, min)
  t1 <- tapply(theta + dt, lab, max)
  tm <- (t0 + t1) / 2
  r.bar <- rep(max(obj$xx) + offset.bar, length(t0))
  r.lab <- rep(max(obj$xx) + offset.lab, length(t0))

  arcs(t0, t1, r.bar, col=col.bar, lwd=2)

  if ( any(!is.na(col.lab)) )
    radial.text(r.lab, tm, sort(unique(lab)),
                col=col.lab, font=font, cex=cex, ...)
}
