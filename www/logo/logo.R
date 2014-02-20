#!/usr/bin/env Rscript
logo <- function(cleanup=FALSE) {
  require(diversitree)
  require(gridBase)
  require(grImport)
  cols <- rev(c("#7f1d80",
                "#005083",
                "#007eb9",
                "#93c13a",
                "#f8d50a",
                "#ed7624",
                "#ca171d"))
  
  if ( !file.exists("beetles.eps") )
    download.file("http://www.zoology.ubc.ca/prog/diversitree/beetles.eps",
                  "beetles.eps")
  if ( cleanup )
    on.exit(file.remove("beetles.eps"))
  
  PostScriptTrace("beetles.eps")
  file.remove("capturebeetles.eps")
  beetles <- readPicture("beetles.eps.xml")
  if ( cleanup )
    on.exit(file.remove(c("beetles.eps", "beetles.eps.xml")))
  
  set.seed(23)
  phy <- ladderize(tree.bd(c(.1, 0), max.taxa=17), FALSE)
  phy$edge.length <- phy$edge.length / max(branching.times(phy)) * .8

  xy <- diversitree:::pp.node.coords(phy)
  xy$xx <- xy$xx + .2
  xy.seg <- diversitree:::pp.coords.phylogram(phy, xy)
  n <- length(phy$tip.label)
  xy.seg$theta0 <- (xy.seg$y0-1)/(n-1) * pi
  xy.seg$theta1 <- (xy.seg$y1-1)/(n-1) * pi
  xy.seg$r0 <- xy.seg$x0
  xy.seg$r1 <- xy.seg$x1

  par(mar=rep(0, 4))
  plot.new()
  plot.window(c(-1.3, 1.3), c(-.15, 1.2), asp=1)
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)

  n <- beetles@summary@numPaths
  xx <- sapply(seq_len(n), function(i) beetles[[i]]@summary@xscale)
  beetles <- beetles[order(-xx[1,])]
  w <- (xx[2,] - xx[1,])[order(-xx[1,])]
  w <- w / max(w) * unit(.2, "native")
  r <- 1.1
  
  for ( i in seq_len(n) ) {
    theta  <- pi / n * (i - 1/2)
    vp <- viewport(x=unit(r * cos(theta), "native"),
                   y=unit(r * sin(theta), "native"),
                   angle=(theta - pi/2) / (2 * pi) * 360)
    pushViewport(vp)
    grid.picture(colour.picture(beetles[[i]], cols[i]),
                 x=unit(.5, "npc"), y=unit(.5, "npc"),
                 just=c("centre", "centre"), width=w[i])
    popViewport()
  }
  popViewport(3)
  diversitree:::pp.segments.fan(phy, xy.seg, "black", 2, 1)
  with(xy.seg[1,], diversitree:::arcs(theta0, theta1, r1, lwd=2, n=200))
  text(0, .03, "diversitree", col="#291441",
       cex=1.8/strwidth("diversitree"), adj=c(.525, 1))
}

logo.hardcopy <- function(filename, type="pdf") {
  if ( missing(filename) )
    filename <- sprintf("logo.%s", type)
  require(Cairo)
  r <- 140 / 200
  w   <- if (type == "png") 250 else 5
  dev <- if (type == "png") CairoPNG else CairoPDF

  CairoFonts(regular="Oswald")
  dev(filename, width=w, height=w*r)
  on.exit(dev.off())
  logo()
}

colour.picture <- function(picture, col) {
  for ( j in seq_along(picture@paths) )
    picture@paths[[j]]@rgb <- col
  picture
}

if ( !interactive() ) {
  logo.hardcopy(type="pdf")
  logo.hardcopy(type="png")
  invisible(file.remove("beetles.eps.xml"))
}
