## In contrast with the other functions in diversitree, we're going to
## use a formula interface here; this will be a little more cumbersome
## for two parameters.

## The "contrasts" based approach comes from Freckleton 2012; thanks
## to Rob Freckleton for clarifying the approach during implementation.

make.pgls <- function(tree, formula, data, control=list()) {
  control <- check.control.pgls(control)
  cache <- make.cache.pgls(tree, formula, data, control)
  if (control$method == "vcv")
    all.branches <- make.all.branches.pgls.vcv(cache, control)
  else if (control$method == "contrasts")
    all.branches <- make.all.branches.pgls.contrasts(cache, control)
  check.pars <- make.check.pars.pgls(cache$np)

  ll <- function(pars) {
    check.pars(pars)
    all.branches(pars)
  }
  ## NOTE: using pgls.dt here to distinguish this from caper's pgls
  ## class so that 'print' works.
  class(ll) <- c("pgls.dt", "dtlik", "function")
  ll
}

make.cache.pgls <- function(tree, formula, data, control) {
  tree <- check.tree(tree, ultrametric=FALSE)

  cache <- check.data.pgls(tree, formula, data)
  cache$n  <- length(tree$tip.label)
  cache$np <- ncol(cache$predictors) + 1 # plus s2

  if (control$method == "vcv") {
    cache$vcv <- vcv(tree)
    cache$vcv.inverse <- solve(cache$vcv)
    cache$vcv.logdet  <- determinant(cache$vcv, logarithm=TRUE)$modulus
    ## Rename in order to match up with the equations better:
    cache$X <- cache$predictors
    cache$Y <- cache$response
  } else {
    pred <- cache$predictors[,-1,drop=FALSE]
    cache$u.x <- apply(pred, 2, pic, tree)
    cache$u.y <- pic(cache$response, tree)
    cache$V <- pic(cache$predictors[,1], tree,
                   var.contrasts=TRUE)[,"variance"]
    cache$V0 <- pgls.root.var.bm(tree)
    cache$root.x <- apply(pred, 2, pgls.root.mean.bm, tree=tree)
    cache$root.y <- pgls.root.mean.bm(tree, cache$response)
  }

  cache$info <- make.info.pgls(tree, cache$predictors)
  cache
}

make.info.pgls <- function(phy, predictors) {
  list(name="pgls",
       name.pretty="PGLS",
       ## Parameters:
       np=3L,
       argnames=default.argnames.pgls(predictors),
       ## Variables:
       ny=NA,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=FALSE,
       ## These are optional
       doc=NULL,
       ## Need some other references in here, too.
       reference=c(
         "Freckleton R.P. 2012 Methods in Ecology and Evolution 3:940-947"))
}
default.argnames.pgls <- function(predictors)
  c(colnames(predictors), "s2")

make.all.branches.pgls.vcv <- function(cache, control) {
  X <- cache$X
  Y <- cache$Y
  n <- cache$n
  vcv.logdet  <- as.numeric(cache$vcv.logdet)
  vcv.inverse <- cache$vcv.inverse
  ## This is a direct implementation of equation 15 in Freckleton
  ## 2012; the only modification is that we make the translation
  ##   log |s2 V| = n log (s2) + log |V|
  ## so that we don't recompute the determinant, and
  ##   (s2 V)^{-1} = V^{-1} / s2
  ## so that we don't recompute the inverse.
  ## It's possible that the
  ##   t(z) %*% (vcv.inverse / s2)
  ## could be replaced by
  ##   crossprod(z, vcv.inverse/s2)
  ## for a small speed saving, but not bothering right now.
  ##
  ## There are also constants (nlog(2pi) + vcv.logdet) that can be
  ## precomputed easily enough.
  function(pars) {
    b <- pars[-length(pars)]
    s2 <- pars[[length(pars)]]
    z <- as.numeric(Y - X %*% b)
    -(n * log(2 * pi) +
      n * log(s2) + vcv.logdet + 
      as.numeric(t(z) %*% (vcv.inverse / s2) %*% (z))) / 2
  }
}

make.all.branches.pgls.contrasts <- function(cache, control) {
  u.x    <- cache$u.x
  u.y    <- cache$u.y
  np     <- cache$np
  n      <- cache$n
  V      <- cache$V
  V0     <- cache$V0
  root.x <- cache$root.x
  root.y <- cache$root.y

  ## Note that n*log(2 * pi) + sum(log(V)) + log(V0) can be factored
  ## out; this is constant every step.  So we'll end up with something
  ## that looks like:
  ##   const <- - 0.5 * (n * log(2 * pi + sum(log(V)) + log(V0)))
  ##   const +
  ##     -0.5 * (n * log(s2) +
  ##             sum((u.y    - u.x %*% b)^2) / s2 +
  ##             (intercept - (root.y - b %*% root.x))^2 / (s2 * V0))
  function(pars) {
    intercept <- pars[[1]]
    b <- pars[-c(1, np)]
    s2 <- pars[[np]]
    ll <- -(n * log(2 * pi) +
            n * log(s2) +
            sum(log(V)) +
            log(V0) +
            sum((u.y - u.x %*% b)^2) / s2) / 2
    dll <- -0.5 * (intercept - (root.y - b %*% root.x))^2 / (s2 * V0)
    as.numeric(ll + dll)
  }
}

check.data.pgls <- function(tree, formula, data, allow.unnamed=FALSE) {
  if (!inherits(formula, "formula"))
    stop("'formula' must be a formua object")
  if (!is.data.frame(data))
    stop("'data' must be a data.frame")

  if (is.null(rownames(data))) {
    if (allow.unnamed) {
      if (nrow(data) == length(tree$tip.label)) {
        rownames(data) <- tree$tip.label
        warning("Assuming data rows are in tree$tip.label order")
      } else {
        stop(sprintf("Invalid data length (expected %d rows)",
                     length(tree$tip.label)))
      }
    } else {
      stop("'data' must contain row names")
    }
  }

  if (!all(tree$tip.label %in% rownames(data)))
    stop("Not all species have state information")
  ## Reorder the data so it matches the tree
  data <- data[match(rownames(data), tree$tip.label),]
  
  ## For now, missing data will straight up cause a fail.
  m <- model.frame(formula, data, na.action=na.fail)
  predictors <- model.matrix(formula, m)
  response <- structure(m[[1]], names=rownames(m))

  list(tree=tree, formula=formula, data=data,
       response=response, predictors=predictors)
}

make.check.pars.pgls <- function(k) {
  function(pars) {
    if (length(pars) != k)
      stop("Incorrect paramter length")
    if (pars[[k]] <= 0)
      stop("Diffusion coefficient must be greater than zero")
    invisible(TRUE)
  }
}

## Simple control checking.
##
## TODO: Allow a switch that would give us REML estimates?  What would
## that mean in the context of Bayesian inference?  See Rob's email
## for some direction here.
check.control.pgls <- function(control) {
  defaults <- list(method="vcv", REML=FALSE)
  control <- modifyList(defaults, control)

  if ( length(control$method) != 1 )
    stop("control$method must be a scalar")
  methods <- c("vcv", "contrasts")
  if (!(control$method %in% methods))
    stop(sprintf("control$method must be in %s",
                 paste(methods, collapse=", ")))
  control
}

pgls.root.var.bm <- function(tree) {
  n.spp <- length(tree$tip.label)
  x <- rep(1, n.spp)
  tmp <- pic(x, tree, var.contrasts=TRUE, rescaled.tree=TRUE)
  idx <- tmp$rescaled.tree$edge[,1] == n.spp + 1
  root.v <- tmp$rescaled.tree$edge.len[idx]
  prod(root.v)/sum(root.v)
}

## Could probably do this better; this is a pretty heavyweight way of
## doing this calculation.  A better way might be to get this in much
## the same way that pgls.root.var.bm works.  However, that uses
## internal features of ape's pic function that I'd rather not go
## near.
##
## TODO: This is going to be worse when we are trying to do this on a
## rescaled tree.  I wonder if we can pull this from the contrasts
## instead?
pgls.root.mean.bm <- function(tree, states) {
  lik <- make.bm(tree, states, control=list(method="pruning", backend="C"))
  attr(lik(1, intermediates=TRUE), "vals")[[1]]
}

fitted.pgls.dt <- function(object, p, ...) {
  b <- p[-length(p)]
  cache <- get.cache(object)
  X <- cache$predictors
  X %*% b
}

fitted.fit.mle.pgls <- function(object, ...) {
  fitted(get.likelihood(object), coef(object, ...))
}

fitted.mcmcsamples.pgls <- function(object, ...) {
  p <- coef(object, ...)
  lik <- get.likelihood(object)
  ret <- apply(p, 1, function(x) fitted(lik, x))
  rownames(ret) <- rownames(diversitree:::get.cache(lik)$predictors)
  ret
}

## Residuals follow from the fitted:
residuals.pgls.dt <- function(object, p, ...) {
  get.cache(object)$response - fitted(object, p, ...)
}
residuals.fit.mle.pgls <- function(object, ...) {
  get.cache(get.likelihood(object))$response - fitted(object, ...)
}
residuals.mcmcsamples.pgls <- residuals.fit.mle.pgls
