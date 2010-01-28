## This is code for plotting phylogenies separate from ape.  However,
## I am aiming to stay fairly compatible with ape plots.  We shall see
## how this goes.
plot2.phylo <- function(x, xlim=NULL, ylim=NULL,
                        show.tip.label=TRUE, show.node.label=FALSE,
                        cex=1, font=3, srt=0, adj=0, label.offset=0,
                        tip.color="black", edge.color="black", ...) {
  if ( length(x$tip.label) < 2 )
    stop("Cannot plot tree with < 2 tips")
  if ( any(tabulate(x$edge[, 1]) == 1) )
    stop("there are single (non-splitting) nodes in your tree;\n",
         "\tYou may need to use collapse.singles()")

  xy <- plot.phylo.node.coords(x)
  xy.seg <- plot.phylo.coords(x, xy)

  plot.phylo.prepare(x, xy, xlim, ylim, cex, show.tip.label,
                     label.offset, ...)
  with(xy.seg, segments(x0, y0, x1, y1))
  if ( show.tip.label )
    plot.phylo.tiplabel(x, xy, label.offset,
                        cex=cex, font=font, srt=srt, adj=adj)
  if ( show.node.label )
    plot.phylo.nodelabel(x, xy, label.offset,
                         cex=cex, font=font, srt=srt, adj=adj)
  plot.phylo.cleanup(x, xy, show.tip.label, show.node.label,
                     font, cex, adj, srt, label.offset, xlim, ylim)
}

## Returns x/y coordinates of the nodes in a phylogeny.  The tree must
## be in pruningwise order.
plot.phylo.node.coords <- function(x) {
  x.p <- reorder(x, "pruningwise")
  if ( is.null(x.p$edge.len) )
    xx <- node.depth(x.p)
  else
    xx <- node.depth.edgelength(x.p)
  yy <- node.height(x.p)
  list(xx=xx, yy=yy)
}

## This returns a matrix with x/y coordinates for the different
## segments of a tree.
plot.phylo.coords <- function(phy, xy) {
  phy.p <- reorder(phy, "pruningwise")
  phy.c <- reorder(phy, "cladewise")
  
  edge <- phy.c$edge
  n.node <- phy$Nnode
  n.tip <- length(phy$tip.label)
  nodes <- (n.tip + 1):(n.tip + n.node)

  ## First, grab the easy ones:
  x0v <- xy$xx[nodes]
  y0v <- y1v <- numeric(n.node)
  x0h <- xy$xx[edge[,1]]
  x1h <- xy$xx[edge[,2]]
  y0h <- xy$yy[edge[,2]]

  for ( i in nodes ) {
    tmp <- range(xy$yy[edge[edge[,1] == i, 2]])
    y0v[i - n.tip] <- tmp[1]
    y1v[i - n.tip] <- tmp[2]
  }

  h <- data.frame(x0=x0h, y0=y0h, x1=x1h, y1=y0h, horiz=TRUE,
                  idx=seq_len(nrow(edge)) + n.tip + n.node,
                  parent=edge[,1])
  v <- data.frame(x0=x0v, y0=y0v, x1=x0v, y1=y1v, horiz=FALSE,
                  idx=seq_len(n.node) + n.tip,
                  parent=h$idx[match(nodes, edge[,2])])
  rbind(v, h)
}

## Basic preparing.
plot.phylo.prepare <- function(x, xy, xlim, ylim, cex,
                               show.tip.label, label.offset, ...) {
  n.tip <- length(x$tip.label)
  if ( is.null(xlim) ) {
    xlim <- c(0, NA)
    ## TODO: Add label offset
    xlab <- if ( show.tip.label ) 
      nchar(x$tip.label) * 0.018 * max(xy$xx) * cex else 0
    xlim[2] <- max(xy$xx[1:n.tip] + xlab)
    xlim[2] <- xlim[2] + label.offset
  }
  if ( is.null(ylim) )
    ylim <- c(1, n.tip)
  
  plot(NA, type="n", xlim=xlim, ylim=ylim, xlab="",
       ylab="", xaxt="n", yaxt="n", bty="n", ...)
}


plot.phylo.tiplabel <- function(x, xy, label.offset, cex, adj, ...) {
  n.tip <- length(x$tip.label)
  wmax <- max(strwidth(x$tip.label, cex = cex))
  text(xy$xx[1:n.tip] + label.offset + wmax * 1.05 * adj,
       xy$yy[1:n.tip], x$tip.label, cex=cex, adj=adj, ...)
}

plot.phylo.nodelabel <- function(x, xy, label.offset, ...) {
  root <- length(x$tip.label) + 1
  text(xy$xx[root:length(xy$xx)] + label.offset,
       xy$yy[root:length(xy$yy)], x$node.label, ...)
}

plot.phylo.nodepoints <- function(x, xy, pch=19, ...) {
  points(as.data.frame(xy)[-(1:length(x$tip.label)),],
         pch=pch, ...)
}

plot.phylo.tippoints <- function(x, xy, pch=19, ...) {
  with(as.data.frame(xy)[1:length(x$tip.label),],
       points(xx + 0.5, yy, pch=pch, ...))
}

plot.phylo.cleanup <- function(x, xy, show.tip.label, show.node.label,
                               font, cex, adj, srt, label.offset,
                               xlim, ylim) {
  ret <- list(use.edge.length=!is.null(x$edge.length), 
              show.tip.label=show.tip.label, 
              show.node.label=show.node.label,
              font=font, cex=cex, adj=adj, srt=srt,
              label.offset=label.offset, 
              x.lim=xlim, y.lim=ylim,
              Ntip=length(x$tip.label), Nnode=x$Nnode,
              xx=xy$xx, yy=xy$yy)
  assign("last_plot.phylo", c(ret, list(edge=x$edge, xx=xy$xx, yy=xy$yy)), 
         envir=.PlotPhyloEnv)
  invisible(ret)
}

## Utility functions.  These replace some C functions within ape.
node.depth <- function(x) {
  x <- reorder(x, "pruningwise")  
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode
  edge1 <- x$edge[,1]
  edge2 <- x$edge[,2]

  xx <- numeric(n.tip + n.node)
  xx[1:n.tip] <- 1
  for ( i in 1:n.edge)
    xx[edge1[i]] = xx[edge1[i]] + xx[edge2[i]]
  xx <- xx - 1
  max(xx) - xx
}

node.depth.edgelength <- function(x) {
  x <- reorder(x, "pruningwise")
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode
  edge1 <- x$edge[,1]
  edge2 <- x$edge[,2]
  edge.length <- x$edge.length

  xx <- numeric(n.tip + n.node)
  for ( i in n.edge:1 )
    xx[edge2[i]] <- xx[edge1[i]] + edge.length[i]
  xx
}

node.height <- function(x) {
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode

  if (!is.null(attr(x, "order")) && attr(x, "order") == "pruningwise")
    x <- reorder(x)

  xx <- reorder(x, order="pruningwise")
  edge1 <- xx$edge[,1]
  edge2 <- xx$edge[,2]

  ## These seem not to be in the right order.
  yy <- numeric(n.tip + n.node)
  tips <- x$edge[x$edge[, 2] <= n.tip, 2]
  yy[tips] <- 1:n.tip

  S <- n <- 0
  for ( i in 1:n.edge ) {
    S <- S + yy[edge2[i]]
    n <- n + 1
    if ( i == n.edge || edge1[i+1] != edge1[i] ) {
      yy[edge1[i]] <- S/n
      S <- n <- 0
    }
  }
  yy
}

## This is currently unused, but sets the vertical bar state/colour
## (which is the same as the node) depending on the state along an
## edge.  If the two edges subtended by a node disagree in
## state/colour, then this returns NA, otherwise they return the
## single unique state/colour.
plot.phylo.node.state <- function(phy.c, edge.state) {
  n.tip <- length(phy.c$tip.label)
  n.node <- phy.c$Nnode
  node.state <- rep(NA, n.node)
  for ( i in seq_len(n.node) ) {
    br <- which(phy.c$edge[, 1] == i + n.tip)
    state <- unique(edge.state[br])
    if ( length(state) == 1 ) 
      node.state[i] <- state
  }
  node.state
}
