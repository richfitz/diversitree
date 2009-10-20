## Plotting function for clade trees.  This duplicates much of the
## original ape function, and is only partially compatible.
plot.clade.tree <- function(x, use.edge.length=TRUE,
                            show.tip.label=TRUE, show.node.label=FALSE,
                            edge.color="black",  edge.width=1, font=3,
                            cex=par("cex"), adj=NULL, srt=0, 
                            root.edge=FALSE, label.offset=0,
                            underscore=FALSE, 
                            x.lim=NULL, y.lim=NULL,
                            tip.color = "black", clade.color="gray",
                            tip.label=x$tip.label,
                            f.clades=sqrt, ...) {
  n.tip <- length(x$tip.label)
  if (n.tip == 1) 
    stop("found only one tip in the tree!")
  n.edge <- nrow(x$edge)
  if (any(tabulate(x$edge[, 1]) == 1)) 
    stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles().")
  n.node <- x$Nnode
  root <- n.tip + 1
  if (is.null(x$edge.length)) 
    use.edge.length <- FALSE

  if (!is.null(attr(x, "order")) && attr(x, "order") == "pruningwise")
    x <- reorder(x)

  edge.color <- rep(edge.color, length.out=n.edge)
  edge.width <- rep(edge.width, length.out=n.edge)

  xe <- x$edge

  x <- reorder(x, order="pruningwise")
  ereorder <- match(x$edge[, 2], xe[, 2])
  edge.color <- edge.color[ereorder]
  edge.width <- edge.width[ereorder]

  if (!use.edge.length)
    xx <- node.depth(x)
  else
    xx <- node.depth.edgelength(x)

  if ( is.null(x$clades) )
    n.taxa <- NULL
  else {
    n.taxa <- f.clades(sapply(x$clades, length))
    names(n.taxa) <- names(x$clades)
  }
  
  yy <- node.height(x, n.taxa)

  if (is.null(x.lim)) {
    x.lim <- c(0, NA)
    tmp <- if (show.tip.label) 
      nchar(x$tip.label) * 0.018 * max(xx) * cex else 0
    x.lim[2] <- max(xx[1:n.tip] + tmp)
  }
  if (is.null(y.lim))
    y.lim <- c(1, sum(n.taxa - 1) + length(x$tip.label))

  plot(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", 
       ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       ...)

  new.phylogram.plot(x, n.tip, n.node, xx, yy,
                     edge.color, edge.width)
  if ( !is.null(x$clades) )
    add.clades(x, xx, yy, clade.color, n.taxa)
  
  if (show.tip.label) {
    if (!underscore) 
      x$tip.label <- gsub("_", " ", x$tip.label)
    if (is.null(adj))
      adj <- 0
    MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
    lox <- label.offset + MAXSTRING * 1.05 * adj
    loy <- 0
    text(xx[1:n.tip] + lox, yy[1:n.tip] + loy, tip.label, 
         adj = adj, font = font, srt = srt, cex = cex, 
         col = tip.color)
  }
  if (show.node.label) 
    text(xx[root:length(xx)] + label.offset, yy[root:length(yy)], 
         x$node.label, adj = adj, font = font, srt = srt, 
         cex = cex)

  L <- list(use.edge.length=use.edge.length, 
            show.tip.label=show.tip.label, 
            show.node.label=show.node.label, font=font, cex=cex, 
            adj=adj, srt=srt, label.offset=label.offset, 
            x.lim=x.lim, y.lim=y.lim,
            tip.color=tip.color, n.tip=n.tip, n.node=n.node,
            n.taxa=n.taxa, xx=xx, yy=yy)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  invisible(L)
}

## Pure R versions of node.height/node.depth etc.
node.height <- function(x, n.taxa=NULL) {
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

  if ( !is.null(n.taxa) ) {
    .yy <- yy
    .yy[tips] <- 1:n.tip

    idx1 <- match(names(n.taxa), x$tip.label)
    n.taxa.per.tip <- rep(1, n.tip)
    n.taxa.per.tip[.yy[idx1]] <- n.taxa

    yy[tips] <- cumsum(n.taxa.per.tip)
    yy[idx1] <- yy[idx1] - (n.taxa-1)/2
  } else {
    yy[tips] <- 1:n.tip
  }

  S <- n <- 0
  for ( i in 1:n.edge ) {
    S <- S + yy[edge2[i]]
    n <- n + 1
    if ( i == n.edge || edge1[i+1] != edge1[i] ) {
      yy[edge1[i]] = S/n
      S <- n <- 0
    }
  }
  yy
}

node.depth <- function(x) {
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
  n.tip <- length(x$tip.label)
  n.edge <- nrow(x$edge)
  n.node <- x$Nnode
  edge1 <- x$edge[,1]
  edge2 <- x$edge[,2]
  edge.length <- x$edge.length

  xx <- numeric(n.tip + n.node)
  for ( i in n.edge:1)
    xx[edge2[i]] = xx[edge1[i]] + edge.length[i];
  xx
}

## Tweaked phylogram.plot function
new.phylogram.plot <- function(x, Ntip, Nnode, xx, yy,
                               edge.color, edge.width) {
  edge <- x$edge
  nodes <- (Ntip + 1):(Ntip + Nnode)
  x0v <- xx[nodes]
  y0v <- y1v <- numeric(Nnode)
  for (i in nodes) {
    j <- edge[which(edge[, 1] == i), 2]
    y0v[i - Ntip] <- min(yy[j])
    y1v[i - Ntip] <- max(yy[j])
  }
  sq <- if (Nnode == 1) 
    1:Ntip
  else c(1:Ntip, nodes[-1])
  y0h <- yy[sq]
  x1h <- xx[sq]
  pos <- match(sq, edge[, 2])
  x0h <- xx[edge[pos, 1]]
  e.w <- unique(edge.width)
  if (length(e.w) == 1) 
    width.v <- rep(e.w, Nnode)
  else {
    width.v <- rep(1, Nnode)
    for (i in 1:Nnode) {
      br <- edge[which(edge[, 1] == i + Ntip), 2]
      width <- unique(edge.width[br])
      if (length(width) == 1) 
        width.v[i] <- width
    }
  }
  e.c <- unique(edge.color)
  if (length(e.c) == 1) 
    color.v <- rep(e.c, Nnode)
  else {
    color.v <- rep("black", Nnode)
    for (i in 1:Nnode) {
      br <- which(edge[, 1] == i + Ntip)
      color <- unique(edge.color[br])
      if (length(color) == 1) 
        color.v[i] <- color
    }
  }
  edge.width <- edge.width[pos]
  edge.color <- edge.color[pos]
  segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v)
  segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width)
}

## n.taxa is the 'scaled' number of taxa.
add.clades <- function(x, xx, yy, col, n.taxa) {
  to   <- match(names(n.taxa), x$tip.label)
  from <- x$edge[match(to, x$edge[,2]),1]

  ## This is hung up on the sister clades.
  ## n.taxa <- sqrt(sapply(x$clades, length))

  x0 <- xx[from]
  x1 <- xx[to]
  ym <- yy[to]
  y0 <- yy[to] - n.taxa/2 + .5
  y1 <- yy[to] + n.taxa/2 - .5

  polygon(rbind(x0, x1, x1, x0, NA), 
          rbind(ym, y0, y1, ym, NA),
          col=col)
}

## Tip-map function, similar to ape::tiplabels()
tipmap <- function(tree, states, plt, cols, offset=0, width=1,
                   offset.y=0) {
  ## To do this, what is the order in which the labels were plotted?
  i <- tree$tip.label[order(plt$yy[1:length(tree$tip.label)])]
  is.clade <- i %in% names(tree$clades)

  ## Categorise the binary data into three classes
  states.3 <- states + 1
  states.3[is.na(states)] <- 3

  yy <- sort(plt$yy[seq_along(tree$tip.label)])

  y0.tip <- yy[!is.clade] - .5
  y1.tip <- y0.tip + 1
  col.tip <- cols[states.3[i][!is.clade]]

  if ( !is.null(tree$clades) ) {
    clades <- tree$clades[match(names(tree$clades), i[is.clade])]
    tab <- t(sapply(clades, function(x)
                    tabulate(states.3[x],3)))
    loc <- t(apply(cbind(0, tab / rowSums(tab) * plt$n.taxa), 1, cumsum)) +
      yy[is.clade] - plt$n.taxa/2
    keep <- as.logical(tab > 0)
    y0.clade <- loc[,-4][keep]
    y1.clade <- loc[,-1][keep]
    col.clade <- cols[rep(1:3, each=length(clades))[keep]]
  } else {
    y0.clade <- y1.clade <- col.clade <- NULL
  }

  x0 <- max(plt$xx) + offset
  x1 <- x0 + width
  y0 <- c(y0.tip, y0.clade) + offset.y
  y1 <- c(y1.tip, y1.clade) - offset.y
  col <- c(col.tip, col.clade)
  rect(x0, y0, x1, y1, col=col, border=NA)
}
