## My simple-minded OU calculator, direct from the SDE:
## 1: make
make.ou <- function(tree, states, states.sd=0, control=list()) {
  if ( missing(control) )
    control <- list(method="vcv")
  control <- check.control.continuous(control)
  if ( control$method == "vcv" )
    stop("VCV-based not yet implemented for OU (consider geiger)")
  cache <- make.cache.ou(tree, states, states.sd)

  ll.ou <- function(pars, root=ROOT.MAX, root.x=NA,
                    intermediates=FALSE) {
    if ( length(pars) != 3 )
      stop("Incorrect parameter length")
    if ( pars[1] < 0 )
      stop("Negative diffusion coefficient")
    if ( pars[2] < 0 )
      stop("Negative alpha coefficient")
    ans <- all.branches.matrix(pars, cache, initial.conditions.ou,
                               branches.ou)
    vals <- ans$init[,cache$root]

    loglik <- root.bm.direct(vals, ans$lq, root, root.x)
    if ( intermediates ) {
      attr(loglik, "intermediates") <- intermediates
      attr(loglik, "vals") <- vals
    }
    loglik
  }

  class(ll.ou) <- c("ou", "function")
  ll.ou
}

## 2: print
print.ou <- function(x, ...) {
  cat("Ornstein Uhlenbeck (OU) likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
argnames.ou <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("diffusion", "alpha", "theta")
  else
    ret
}
`argnames<-.ou` <- function(x, value) {
  if ( length(value) != 3 )
    stop("Incorrect argument length: expected 3")
  attr(x, "argnames") <- value
  x
}

## 4: find.mle
find.mle.ou <- function(func, x.init, method,
                           fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.ou")
}

## 5: make.cache
make.cache.ou <- make.cache.bm.direct

## 6: ll

## 7: initial.conditions
initial.conditions.ou <- initial.conditions.bm.direct

## 8: branches
branches.ou <- function(y, len, pars, t0, idx) {
  m <- y[1]
  v <- y[2]
  z <- y[3]

  sigma2 <- pars[1]
  alpha  <- pars[2]
  theta  <- pars[3]
  ## In ou, for a branch that ends with mean m (y[1]), the mean moves
  ## to
  ##   exp(len * alpha) * (m - theta) + theta
  ## And the variance becomes
  ##   (exp(2 * len * alpha) - 1) / (2 * alpha) + exp(2 * len * alpha) * v
  ## With normalising constant
  ##   exp(t * alpha)

  ## The last line of the second option comes from the limit of
  ##   (exp(2*len*alpha) - 1) * sigma2 / (2*alpha)
  ## as alpha -> 0 being len * sigma2
  ##   if ( alpha > 0 )
  ##     c(len * alpha + z,
  ##       exp(len * alpha) * (m - theta) + theta,
  ##       (exp(2*len*alpha) - 1) * sigma2 / (2*alpha) + exp(2*len*alpha) * v)
  ##   else
  ##     c(len * alpha + z,
  ##       exp(len * alpha) * (m - theta) + theta,
  ##       len * sigma2 + exp(2*len*alpha) * v)

  if ( alpha > 0 )
    list(len * alpha + z,
         c(exp(len * alpha) * (m - theta) + theta,
           (exp(2*len*alpha) - 1) * sigma2 / (2*alpha) +
           exp(2*len*alpha) * v,
           0))
  else
    list(len * alpha + z,
         c(exp(len * alpha) * (m - theta) + theta,
           len * sigma2 + v,
           0))
}
