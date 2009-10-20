## I should use stats4:::mle, but this does not work well for me with
## computing the variance-covariance matrix.
find.mle <- function(func, x.init, ...) {
  UseMethod("find.mle")
}

## TODO: Bug here (see 20091019 email with Andy Wilson) where
##   find.mle(f, x, nosucharg=1)
## will return nonsense, because the fail.value kicks in.  Check that
## ... arguments!
## I think that I have fixed this now, but need to check.
find.mle.bisse <- function(func, x.init,
                           method=c("L-BFGS-B", "Nelder-Mead", "subplex"),
                           control=list(), lower=NULL, upper=NULL,
                           fail.value=NULL, hessian=FALSE, ...) {
  method <- match.arg(method)
  names <- argnames(func)
  npar <- length(names)

  if ( inherits(func, c("constrained", "fixed")) ) {
    ## Identify the parameters we do have:
    arg.idx <- match(names, argnames(environment(func)$f))
    if ( length(x.init) == 6 ) x.init <- x.init[arg.idx]
    if ( length(lower) == 6 )  lower <- lower[arg.idx]
    if ( length(upper) == 6 )  upper <- upper[arg.idx]
  } else {
    arg.idx <- 1:6
  }

  if ( is.null(names(x.init)) )
    names(x.init) <- names

  if ( method == "subplex" ) {
    if ( !require(subplex) )
      stop("The subplex package is required")
    if ( is.null(lower) ) lower <- rep(0, npar)
    if ( is.null(upper) ) upper <- rep(Inf, npar)
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- -Inf
    control <- modifyList(list(reltol=.Machine$double.eps^0.25),
                          control)

    func2 <- function(x) {
      if ( any(x < lower | x > upper) )
        -fail.value
      else
        -func(x, ..., fail.value=fail.value)
    }
    ans <- subplex(x.init, func2, control=control)
    ans$value <- -ans$value
    ans$hessian <- NULL
  } else {
    if ( is.null(lower) )
      lower <- c(1e-4, 1e-4, 0, 0, 1e-4, 1e-4)[arg.idx]
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- func(x.init, ...) - 1000
    dx <- 1e-5
    control <- modifyList(list(fnscale=-1, ndeps=rep(dx, npar)),
                          control)
    ans <- optim(x.init, func, control=control, lower=lower,
                 method="L-BFGS-B", fail.value=fail.value,
                 ...)
  }
  
  names(ans)[names(ans) == "value"] <- "lnLik"
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle: ",
            tolower(ans$message))
  class(ans) <- "mle.bisse"

  at.edge <- ans$par - lower == 0 & lower > 0
  if ( any(at.edge) )
    warning(sprintf("Parameter(s) %s at edge of parameter space",
                    paste(names[at.edge], collapse=", ")))
  if ( hessian ) {
    if ( !require(numDeriv) )
      stop("The package numDeriv is required to compute the hessian")
    ans$hessian <- hessian(func, ans$par, ...)
  }

  ans
}

logLik.mle.bisse <- function(object, ...) {
  ll <- object$lnLik
  attr(ll, "df") <- length(object$par)
  class(ll) <- "logLik"
  ll
}

## Code based on MASS:::anova.negbin and ape:::anova.ace
anova.mle.bisse <- function(object, ...) {
  mlist <- c(list(object), list(...))
  if ( length(mlist) == 1L )
    stop("Need to specify more than one model")
  if ( is.null(names(mlist)) )
    names(mlist) <-
      c("full", model=sprintf("model %d", seq_len(length(mlist)-1)))
  else
    names(mlist)[1] <- "full"

  ll <- lapply(mlist, logLik)
  ll.val <- sapply(ll, as.numeric)
  chisq <- c(NA, abs(2*(ll.val[1] - ll.val[-1])))
  df <- sapply(ll, attr, "df")
  ddf <- c(NA, abs(df[1] - df[-1]))
  
  out <- data.frame(Df=df,
                    lnLik=sapply(ll, as.numeric),
                    AIC=sapply(mlist, AIC),
                    ChiSq=chisq,
                    "Pr(>|Chi|)"=1 - pchisq(chisq, ddf),
                    check.names=FALSE)
  rownames(out) <- names(mlist)
    
  class(out) <- c("anova", "data.frame")
  out
}
