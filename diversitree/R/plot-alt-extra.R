## TODO: phylogram version.

## obj: result of plot2.phylo
## lab: vector along the tip labels with a "group" level.  There will
## be repetition.
## col.bar and col.lab: vectors along sort(unique(lab)) with colours
## for the bar and label.
group.label.tip.rad <- function(obj, lab, col.bar, col.lab, lwd=1,
                                offset.bar=0, offset.lab=0, cex=1,
                                font=1, check=FALSE, quiet=FALSE,
                                ...) {
  n.taxa <- obj$n.taxa
  n <- obj$n.spp
  if ( is.null(n.taxa) )
    dt <- 1/6 / n * 2 * pi
  else
    dt <- (n.taxa/2 - .5 + 1/6) / n * 2 * pi

  theta <- obj$xy$theta[seq_len(obj$Ntip)]

  t0 <- tapply(theta - dt, lab, min)
  t1 <- tapply(theta + dt, lab, max)
  str <- names(t0)

  if ( check ) {
    ## Sort these in order around the circle (breaks vector colour/fonts
    ## though, which I think I used somewhere)
    i <- order(t0)
    t0 <- t0[i]
    t1 <- t1[i]
    str <- str[i]
    
    g <- integer(length(t0))
    g[1] <- j <- 1
    end <- t1[1]
    for ( i in seq_along(g)[-1] ) {
      if ( t0[i] > end ) {
        j <- j + 1
        end <- t1[i]
      } else {
        end <- max(end, t1[i])
      }
      g[i] <- j
    }

    tg <- table(g)
    if ( any(tg > 1) ) {
      if ( !quiet ) {
        err <- sapply(which(tg != 1), function(x)
                      paste(str[g == x], collapse=", "))
        warn <- c("Collapsing non-monophyletic groups:",
                  sprintf("\t%s", err))
        warning(paste(warn, collapse="\n"))
      }
      t0 <- tapply(t0, g, min)
      t1 <- tapply(t1, g, max)
      str <- as.character(tapply(str, g, collapse))
    }
  }

  tm <- (t0 + t1) / 2
  
  r.bar <- rep(max(obj$xx) + offset.bar, length(t0))
  r.lab <- rep(max(obj$xx) + offset.lab, length(t0))

  arcs(t0, t1, r.bar, col=col.bar, lwd=lwd)

  if ( any(!is.na(col.lab)) )
    radial.text(r.lab, tm, str,
                col=col.lab, font=font, cex=cex, ...)
}

trait.plot <- function(tree, dat, cols, lab=names(cols), str=0:1,
                       class=NULL, type="f", w=1/50,
                       legend=length(cols) > 1, cex.lab=.5,
                       font.lab=3, cex.legend=.75, margin=1/4,
                       check=TRUE, quiet=FALSE) {
  if ( type != "f" )
    stop("type != f not yet implemented")
  if ( is.null(names(cols)) )
    stop("'cols' must be named")
  if ( !all(names(cols) %in% names(dat)) )
    stop("Not all colours have data")
  if ( !is.null(class) && length(class) != length(tree$tip.label) )
    stop("'class' must be a vector along tree$tip.label")
  n <- length(cols)
  if ( n < 1 )
    stop("Need some colours")
  if ( !is.data.frame(dat) ) {
    if ( is.vector(dat) && n == 1 ) {
      dat <- matrix(dat)
      rownames(dat) <- names(dat)
    } else {
      stop("dat must be a matrix")
    }
  }
  if ( !all(tree$tip.label %in% rownames(dat)) )
    stop("All taxa must have entries in 'dat' (rownames)")
  dat <- dat[tree$tip.label,,drop=FALSE]

  par(mar=rep(0, 4))
  t <- max(branching.times(tree))
  w <- w * t
  if ( is.null(class) ) {
    plt <- plot2.phylo(tree, type="f", show.tip.label=TRUE,
                       label.offset=(n+2)*w, cex=cex.lab)
  } else {
    plt <- plot2.phylo(tree, type="f", show.tip.label=FALSE,
                       label.offset=t*margin)
    group.label.tip.rad(plt, class, "black", "black",
                        offset.bar=w*4, offset.lab=w*5, lwd=1.5,
                        cex=cex.lab, font=font.lab,
                        check=check, quiet=quiet)
  }

  xy <- plt$xy
  theta <- xy$theta[seq_along(tree$tip.label)]
  dt <- diff(sort(theta))[1]/2

  for ( i in seq_along(cols) ) {
    v <- names(cols)[i]
    filled.arcs(theta - dt, theta + dt, max(xy$x) + i * w, w,
                cols[[v]][dat[[v]]+1])
  }

  if ( legend ) {
    leg <- legend("topright", legend=rep(str, each=n), ncol=2, bty="n",
                  fill=c(do.call(rbind, cols)), cex=cex.legend)
    text(leg$rect$left, leg$text$y[1:n], sprintf("%s:", lab), adj=1,
         cex=cex.legend)
  }
}
environment(trait.plot) <- environment(make.bisse)

collapse <- function(x) {
  if ( length(x) == 1 )
    x
  else if ( length(x) == 2 )
    paste(x, collapse=" &\n")
  else
    sprintf("%s,\n& %s",
            paste(x[-length(x)], collapse=",\n"),
            x[length(x)])
}
