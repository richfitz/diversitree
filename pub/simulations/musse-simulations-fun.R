## Simulate some trees under MuSSE/multitrait, where state 1 in trait
## A doubles the speciation rates.  All traits are uncorrelated with
## transition rates of 0.01 (1/10 of the speciation rate).
sim.trees <- function(n.trait, n.taxa, n.trees, regenerate=FALSE,
                      seed=1) {
  filename <- sprintf("trees/trees-%d-%d.Rdata", n.trait, n.taxa)
  if ( file.exists(filename) && !regenerate )
    return(load.local(filename))
  
  tr <- musse.multitrait.translate(n.trait, 1)
  p <- rep(0, ncol(tr))
  names(p) <- colnames(tr)
  p["lambda0"] <- .1
  p["lambdaA"] <- .1  # double diversification rate in state A
  p["mu0"]     <- .03 # same as BiSSE paper
  p[sprintf("q%s01.0", LETTERS[1:n.trait])] <- .01
  p[sprintf("q%s10.0", LETTERS[1:n.trait])] <- .01

  p2 <- drop(tr %*% p)

  n.trees <- 100
  lambda <- p2[seq_len(2^n.trait)]
  x0 <- sample(as.integer(which(lambda == .1)), n.trees, TRUE)

  key <- do.call(expand.grid, rep(list(0:1), n.trait))
  names(key) <- LETTERS[1:n.trait]
  key <- as.matrix(key)

  f <- function(x0) {
    phy <- trees(p2, "musse", 1, max.taxa=n.taxa, x0=x0)[[1]]
    states <- as.data.frame(key[phy$tip.state,,drop=FALSE])
    colnames(states) <- LETTERS[seq_len(n.trait)]
    rownames(states) <- phy$tip.label
    phy$tip.state <- states
    phy
  }

  set.seed(seed)
  trees <- lapply(x0, f)
  save(trees, file=filename)
  trees
}

fit.ml <- function(trees, ..., regenerate=FALSE) {
  n.trait <- ncol(trees[[1]]$tip.state)
  n.taxa <- length(trees[[1]]$tip.label)
  filename <- sprintf("output/ml-%d-%d.Rdata", n.trait, n.taxa)
  if ( file.exists(filename) && !regenerate )
    return(load.local(filename))

  fits <- mclapply(trees, fit.ml.one, ...)
  save(fits, file=filename)
  fits
}


## Models to fit:

## (a) Constant extinction:
##   No SSD
##   SSD(A), ..., SSD(n)
##   SSD(A) & ... & SSD(n)
fit.ml.one <- function(x) {
  states <- x$tip.state
  n.trait <- ncol(states)

  ## Constrain lambdas in turn, and fit constrained models:
  lik0 <- make.musse.multitrait(x, states, depth=c(1, 0, 0),
                                strict=FALSE, # daring...
                                control=list(backend="CVODES"))
  p0 <- starting.point.musse.multitrait(x, lik0)

  formulae <- lapply(sprintf("lambda%s ~ 0", LETTERS[1:n.trait]),
                     as.formula, .GlobalEnv)

  ## No SDD:
  lik0.0 <- constrain(lik0, formulae=formulae)
  fit0.0 <- find.mle(lik0.0, p0[argnames(lik0.0)])

  ## Include only one trait:
  f <- function(i) {
    lik <- constrain(lik0, formulae=formulae[-i])
    find.mle(lik, p0[argnames(lik)])
  }
  fit0.1 <- lapply(1:n.trait, f)

  ## Include all traits:
  fit0.2 <- find.mle(lik0, p0)

  list(fit0.0, fit0.1, fit0.2)
}
 
fit.mcmc <- function(trees, ..., incorrect=FALSE, regenerate=FALSE) {
  n.trait <- ncol(trees[[1]]$tip.state)
  n.taxa <- length(trees[[1]]$tip.label)
  if ( incorrect )
    filename <- sprintf("output/mcmc-incorrect-%d-%d.Rdata", n.trait, n.taxa)
  else
    filename <- sprintf("output/mcmc-%d-%d.Rdata", n.trait, n.taxa)
  
  if ( file.exists(filename) && !regenerate )
    return(load.local(filename))

  fits <- mclapply(trees, fit.mcmc.one, ..., incorrect=incorrect)
  save(fits, file=filename)
  fits
}

fit.mcmc.one <- function(x, print.every=100, incorrect=FALSE) {
  states <- x$tip.state
  if ( incorrect )
    states <- states[,-1,drop=FALSE]
  n.trait <- ncol(states)

  ## Constrain lambdas in turn, and fit constrained models:
  lik <- make.musse.multitrait(x, states, depth=c(1, 0, 0),
                               strict=FALSE, # daring...
                               control=list(backend="CVODES"))
  p <- starting.point.musse.multitrait(x, lik)
  prior <- make.prior.multitrait(x, lik)

  set.seed(1)
  tmp <- mcmc(lik, p, 100, w=.1, prior=prior, print.every=print.every)
  w <- diff(apply(tmp[-c(1, ncol(tmp))], 2, quantile, c(1,19)/20))[1,]
  w <- pmax(w, .1) # never drop below .1
  samples <- mcmc(lik, p, 10000, w=.1, prior=prior,
                  print.every=print.every)
  ## Discard index and p to save a little space, also drop first 500
  ## points as burn-in (should be more than sufficient).
  samples[-(1:500),-c(1,ncol(samples))]
}

## Wrappers to make running these easiy.
run.ml <- function(n.trait, n.taxa, regenerate=FALSE) {
  cat(sprintf("Running %d/%d\n", n.trait, n.taxa))
  cat("\tTrees\n")
  trees <- sim.trees(n.trait, n.taxa, regenerate)
  cat("\tFits\n")
  fit.ml(trees, mc.preschedule=FALSE, regenerate=regenerate)
}

run.mcmc <- function(n.trait, n.taxa, regenerate=FALSE,
                     incorrect=FALSE) {
  cat(sprintf("Running %d/%d\n", n.trait, n.taxa))
  cat("\tTrees\n")
  trees <- sim.trees(n.trait, n.taxa, regenerate)
  cat("\tFits\n")
  fit.mcmc(trees, mc.preschedule=FALSE,
           regenerate=regenerate, incorrect=incorrect)
}

## Plot the results:
fig.musse.mcmc <- function(nn, obj) {
  if ( nrow(nn) != length(obj) )
    stop("Length mismatch")

  obj.A <- lapply(obj, sapply, function(x) x[,"lambdaA"])
  obj.B <- lapply(obj, sapply, function(x)
                  if ("lambdaB" %in% colnames(x)) x[,"lambdaB"] else
                  c(NA, NA, NA))

  cols <- c("#1f78b4", "#b2df8a")
  cols.tr <- add.alpha(cols, .25)
  cols.dk <- mix(cols, "black", .5)
  ylim <- c(-.1, .2)
  n.taxa <- sort(unique(nn$n.taxa))  
  par(mfrow=c(2,2), oma=c(4.1, 4.1, 0, 0), mar=rep(.5, 4))
  for ( nt in sort(unique(nn$n.trait)) ) {
    i <- nn$n.trait == nt
    
    sub.means.A <- sapply(obj.A[i], function(x) mean(x[1,]))
    sub.means.B <- sapply(obj.B[i], function(x) mean(x[1,]))
    env.A <- t(sapply(obj.A[i], function(x) apply(x[-1,], 1, median)))
    env.B <- t(sapply(obj.B[i], function(x) apply(x[-1,], 1, median)))

    plot(NA, xlim=range(n.taxa), ylim=ylim, las=1, log="x",
         xaxt="n", yaxt="n")
    if ( nt %in% c(3, 4) ) axis(1, n.taxa, col=NA, col.ticks="black")
    if ( nt %in% c(1, 3) ) axis(2, las=1, col=NA, col.ticks="black")
    abline(h=c(0, .1), col="darkgrey", lty=3, lwd=2)
    str <- sprintf("%s) %d %s", letters[nt], nt,
                   if (nt == 1) "trait" else "traits")
    text(n.taxa[1], ylim[2], str, adj=c(0, 1))

    polygon(c(n.taxa, rev(n.taxa)), c(env.B[,1], rev(env.B[,2])),
            col=cols.tr[2], border=NA)
    polygon(c(n.taxa, rev(n.taxa)), c(env.A[,1], rev(env.A[,2])),
            col=cols.tr[1], border=NA)
    if ( nt > 1 )
      lines(n.taxa, sub.means.B, lwd=2, col=cols.dk[2])
    lines(n.taxa, sub.means.A, lwd=2, col=cols.dk[1])
  }
  mtext("Number of species", 1, outer=TRUE, line=2.25)
  mtext("Speciation rate main effect", 2, outer=TRUE, line=3)
}

fig.power <- function(nn, obj, obj.inc) { 
  ## Significant fits have a ci that has the same sign, so
  ## sign(prod(ci)) is positive.
  is.sig <- function(x, v) if (v %in% colnames(x))
    sign(prod(x[2:3,v])) > 0 else NA
  ## Corectly identify A as significant:
  p.cor <- sapply(obj, function(x)
                  mean(sapply(x, is.sig, "lambdaA")))
  ## Incorrectly also identify B as significant:
  p.inc1 <- sapply(obj, function(x)
                   mean(sapply(x, is.sig, "lambdaB")))
  p.inc2 <- c(rep(NA, 4),
              sapply(obj.inc, function(x)
                     mean(sapply(x, is.sig, "lambdaB"))))

  s.b <- sapply(obj.inc, function(x)
                sapply(x, is.sig, "lambdaB"))
  s.c <- sapply(obj.inc, function(x)
                sapply(x, is.sig, "lambdaC"))
  s.d <- sapply(obj.inc, function(x)
                sapply(x, is.sig, "lambdaD"))
  s.b[is.na(s.b)] <- s.c[is.na(s.c)] <- s.d[is.na(s.d)] <- FALSE
  p.inc3 <- c(rep(NA, 4), colMeans(s.b | s.c | s.d))
  
  cols <- c("#1f78b4", "#33a02c", "#ff7f00")
  cols.tr <- add.alpha(cols, .25)
  cols.dk <- mix(cols, "black", .5)

  par(mfrow=c(2,2), oma=c(4.1, 4.1, 0, 0), mar=rep(.5, 4))
  ylim <- c(0, 1)
  n.taxa <- sort(unique(nn$n.taxa))
  for ( nt in sort(unique(nn$n.trait)) ) {
    plot(NA, xlim=range(n.taxa), ylim=ylim, las=1, log="x",
         xaxt="n", yaxt="n")
    abline(h=1/20, col="darkgrey", lty=3)
    if ( nt %in% c(3, 4) ) axis(1, n.taxa, col=NA, col.ticks="black")
    if ( nt %in% c(1, 3) ) axis(2, las=1, col=NA, col.ticks="black")
    str <- sprintf("%s) %d %s", letters[nt], nt,
                   if (nt == 1) "trait" else "traits")
    text(n.taxa[1], ylim[2], str, adj=c(0, 1))

    i <- nn$n.trait == nt  
    if ( nt > 1 ) {
      lines(n.taxa, p.inc1[i], lwd=2, col=cols[2])
      lines(n.taxa, p.inc2[i], lwd=2, col=cols[2], lty=2)
    }
    if ( nt > 2 )
      lines(n.taxa, p.inc3[i], lwd=2, col=cols[3], lty=3)
    lines(n.taxa, p.cor[i], lwd=2, col=cols[1])
  }
  mtext("Number of species", 1, outer=TRUE, line=2.25)
  mtext("Proportion of trees significant", 2, outer=TRUE, line=3)
}

## Utility functions:
load.local <- function(x) {
  v <- load(x)
  if ( length(v) != 1 )
    stop(".Rdata file must contain exactly one object")
  get(v)
}

make.prior.multitrait <- function(x, lik) {
  p <- starting.point.bd(x)
  r <- -diff(p)[[1]]
  prior1 <- make.prior.exponential(1/(2*r))  
  function(pars)
    prior1(lik(pars, pars.only=TRUE))
}

## Plotting utilities
hdr.uniroot <- diversitree:::hdr.uniroot
add.alpha <- diversitree:::add.alpha

mix <- function(cols, col2, p) {
  m <- col2rgb(cols)
  m2 <- col2rgb(rep(col2, length=length(cols)))
  m3 <- (m * p + m2 * (1-p))/255
  rgb(m3[1,], m3[2,], m3[3,])
}

to.pdf <- function(filename, width, height, expr,
                   ..., pointsize=10, verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}
