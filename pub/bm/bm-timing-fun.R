timeit.one <- function(phy, states, s2, method, n=1, t.target=1,
                       states.sd=NULL) {
  if ( method == "vcv" )
    control <- list(method="vcv")
  else if ( method == "pruning.R" )
    control <- list(method="pruning", backend="R")
  else if ( method == "pruning.C" )
    control <- list(method="pruning", backend="C")
  else
    stop()
  lik <- make.bm(phy, states, states.sd, control=control)

  repeat {
    t <- sum(system.time(replicate(n, lik(s2)))[1:2])
    if ( max(t) > t.target )
      return(c(n=n, t=t))
    else
      n <- n * 2
  }
}

timeit <- function(n.taxa, n.rep, method, seed=1, s2=.1, t.target=1,
                   states.sd=NULL) {
  cat(sprintf("%s: %d\n", method, n.taxa))
  set.seed(seed)
  trees <- trees(.1, "yule", n=n.rep, max.taxa=n.taxa)
  states <- lapply(trees, sim.character, s2)

  tmp <- timeit.one(trees[[1]], states[[1]], s2, method, 1, t.target)
  n <- tmp[["n"]]
  ans <- sapply(seq_len(n.rep)[-1], function(i)
                timeit.one(trees[[i]], states[[i]], s2, method, n,
                           t.target, states.sd))
  ans <- rbind(tmp, t(ans), deparse.level=0)
  cbind(ans, t.each=ans[,"t"] / ans[,"n"])
}

timings <- function(n.taxa, n.rep, seed=1, s2=.1, t.target=1,
                    states.sd=NULL) {
  res.vcv    <- lapply(n.taxa, timeit, n.rep, "vcv",
                       t.target=t.target, states.sd=states.sd)
  res.pruning.R <- lapply(n.taxa, timeit, n.rep, "pruning.R",
                         t.target=t.target, states.sd=states.sd)
  res.pruning.C <- lapply(n.taxa, timeit, n.rep, "pruning.C",
                         t.target=t.target, states.sd=states.sd)
  list(vcv=res.vcv, pruning.R=res.pruning.R, pruning.C=res.pruning.C)
}

load.local <- function(x) {
  v <- load(x)
  if ( length(v) != 1 )
    stop(".Rdata file must contain exactly one object")
  get(v)
}

## Simple caching; run 'expr' and save the output in 'filename'; if
## 'filename' already exists just load that.  If regenerate is TRUE,
## it always runs the expression.
run.cached <- function(filename, expr, regenerate=FALSE) {
  if ( file.exists(filename) && !regenerate ) {
    load.local(filename)
  } else {
    res <- eval.parent(substitute(expr))
    save(res, file=filename)
    res
  }
}

to.pdf <- function(filename, width, height, expr,
                   ..., pointsize=10, verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

## fig.timing <- function(n.taxa, obj) {
##   times <- sapply(obj, sapply, function(x) mean(x[,"t.each"]))
##   cols <- c("#ff7f00", "#1f78b4", "#33a02c")
##   pch <- c(19, 15, 18)
##   par(mar=c(4.1, 4.1, .5, .5))
##   matplot(n.taxa, times, log="xy", type="o", pch=pch, col=cols, lty=1,
##           xaxt="n", yaxt="n",
##           xlab="Number of taxa",
##           ylab="Elapsed time (s)")
##   axis(1, n.taxa, n.taxa, col=NA, col.ticks="black")
##   r <- par("usr")[3:4]
##   at <- pretty(r)
##   lab <- do.call(expression, lapply(at, function(i) bquote(10^.(i))))
##   axis(2, at=10^at, lab, col=NA, col.ticks="black", las=1)
##   legend("topleft", c("vcv", "pruning/R", "pruning/C"),
##          col=cols, pch=pch, lty=1, bty="n", inset=c(1/60, 0))
## }

fig.timing <- function(n.taxa, times, times.err) {
  ylim <- range(times, times.err)
  cols <- c("#ff7f00", "#1f78b4", "#33a02c")
  pch <- c(19, 15, 18)
  par(mar=c(4.1, 4.1, .5, .5))
  matplot(n.taxa, times, type="o", pch=pch, col=cols, lty=1,
          log="xy", ylim=ylim, xaxt="n", yaxt="n",
          xlab="Number of taxa",
          ylab="Elapsed time (s)")
  matlines(n.taxa, times.err, type="o", pch=pch, col=cols, lty=2)
  axis(1, n.taxa, n.taxa, col=NA, col.ticks="black")
  r <- par("usr")[3:4]
  at <- pretty(r)
  lab <- do.call(expression, lapply(at, function(i) bquote(10^.(i))))
  axis(2, at=10^at, lab, col=NA, col.ticks="black", las=1)
  legend("topleft", c("vcv", "pruning/R", "pruning/C"),
         col=cols, pch=pch, lty=1, bty="n", inset=c(1/60, 0))
}


