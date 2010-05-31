## Currently, I have added interfaces to several different optimisers:
##   'optim' (L-BFGS-B, Nelder-Mead, BFGS, CG, SANN)
##   'subplex'
##   'nlminb'
##   'nlm'
## These now have a similar interface, though not identical.

## TODO: consider making stats4 mle objects, rather than those below:
## new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef), 
##     vcov = vcov, min = min, details = oout, minuslogl = minuslogl, 
##     method = method)

## TODO: These all assume that the function is "protectable" by
## passing in 'fail.value' - this is not always the case.  However, if
## I move the "protect" call *inside* of find.mle/mcmc this infelicity
## goes away.

find.mle <- function(func, x.init, method, ...) {
  UseMethod("find.mle")
}

find.mle.default <- function(func, x.init, method,
                             fail.value=NA, class.append=NULL, ...) {
  if ( inherits(func, "constrained") )
    x.init <- guess.constrained.start(func, x.init)

  ans <- do.mle.search(func, x.init, method, fail.value, ...)
  class(ans) <- c(class.append, class(ans))
  ans
}

do.mle.search <- function(func, x.init, method, fail.value=NA,
                          hessian=FALSE, verbose=0,
                          ...) {
  method <- match.arg(method, c("optim", "subplex", "nlminb", "nlm"))

  if ( verbose > 0 )
    func2 <- big.brother(protect(func), verbose)
  else
    func2 <- protect(func)

  if ( is.null(names(x.init)) ) {
    names.v <- try(argnames(func), silent=TRUE)
    if ( !inherits(names.v, "try-error") ) {
      if ( length(names.v) != length(x.init) )
        stop("Invalid parameter names length: expected ",
             length(names.v), ", got ", length(x.init))
      names(x.init) <- names.v
    }
  }

  ## TODO: check starting position like so.  However, doing this
  ## fails because things like lower are passed through "...".
  ## y0 <- try(func(x.init, ...))
  ## if ( inherits(y0, "try-error") )
  ##   stop("The above error occured when testing the starting position")
  ## if ( !is.finite(y0) )
  ##   stop("Starting point must have finite probability")

  mle.search <- get(sprintf("do.mle.search.%s", method))
  ans <- mle.search(func2, x.init, fail.value, ...)
  if ( verbose )
    cat("\n")

  if ( hessian ) {
    if ( !require(numDeriv) )
      stop("The package numDeriv is required to compute the hessian")
    ans$hessian <- hessian(func, ans$par, ...)
  }

  ans$func <- func
  class(ans) <- "fit.mle"
  ans
}

do.mle.search.optim <- function(func, x.init, fail.value=NA,
                                control=list(), lower=-Inf, upper=Inf,
                                dx=1e-5, optim.method="L-BFGS-B",
                                ...) {
  if ( is.null(fail.value) || is.na(fail.value) )
    fail.value <- func(x.init, ...) - 1000
  control <- modifyList(list(fnscale=-1,
                             ndeps=rep(dx, length(x.init))), control)
  ans <- optim(x.init, func, fail.value=fail.value,
               control=control, lower=lower, upper=upper,
               method=optim.method, ...)
  names(ans)[names(ans) == "value"] <- "lnLik"

  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (optim): ",
            tolower(ans$message))

  ans$method <- "optim"
  ans$optim.method <- optim.method
  ans
}

do.mle.search.subplex <- function(func, x.init, fail.value=NA,
                                  control=list(), lower=-Inf,
                                  upper=Inf, ...) {
  if ( is.null(fail.value) || is.na(fail.value) )
    fail.value <- -Inf
  ## By default, a less agressive tolerance that is more likely to be
  ## met.
  control <- modifyList(list(reltol=.Machine$double.eps^0.25),
                        control)

  ## This is fairly simple minded treatment of box constraints
  ## compared with optim.  In particular, no check is made to see if
  ## the best solution falls at the edge of parameter space.
  ## Realistically, this would not be too difficult to do (though
  ## fairly expensive); we could check at the end of the optimisation
  ## if the best point was close to a boundary (say, 10 * control
  ## within boundaries) and set those points to their boundary
  ## values.  Constrain the problem and start again.  I'm not sure if
  ## this is worth the effort though.
  if ( any(is.finite(lower)) || any(is.finite(upper)) ) {
    if ( any(x.init < lower | x.init > upper) )
      stop("Starting point falls outside of box constraints")
    func2 <- invert(boxconstrain(func, lower, upper))
  } else {
    func2 <- invert(func)
  }
  ans <- subplex(x.init, func2, fail.value=fail.value,
                 control=control, ...)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"

  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (subplex): ",
            tolower(ans$message))
  
  ans$method <- "subplex"  
  ans
}

do.mle.search.nlminb <- function(func, x.init, fail.value=NA,
                                 control=list(), lower=-Inf,
                                 upper=Inf, ...) {
  if ( is.null(fail.value) || is.na(fail.value) )
    fail.value <- -Inf
  ans <- nlminb(x.init, invert(func), fail.value=fail.value,
                control=control, lower=lower, upper=upper, ...)
  names(ans)[names(ans) == "objective"] <- "lnLik"
  names(ans)[names(ans) == "evaluations"] <- "counts"
  ans$lnLik <- -ans$lnLik
  ans <- ans[c("par", "lnLik", "counts", "convergence", "message",
               "iterations")]
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (nlminb): ",
            tolower(ans$message))
  ans$method <- "nlminb"
  ans
}

do.mle.search.nlm <- function(func, x.init, fail.value=NA,
                              control=list(), lower=-Inf,
                              upper=Inf, ...) {
  if ( is.null(fail.value) || is.na(fail.value) )
    fail.value <- func(x.init, ...) - 1000    

  ans <- nlm(invert(func), x.init, fail.value=fail.value, ...)

  names(ans)[names(ans) == "estimate"] <- "par"
  names(ans)[names(ans) == "minimum"] <- "lnLik"
  names(ans)[names(ans) == "iterations"] <- "counts"
  ans$lnLik <- -ans$lnLik 
  names(ans$par) <- names(x.init)
  ans <- ans[c("par", "lnLik", "counts", "code", "gradient")]
  
  if ( ans$code > 2 )
    warning("Convergence problems in find.mle: code = ",
            ans$code, " (see ?nlm for details)")  

  ans$method <- "nlm"
  ans
}

## For want of a better name, this does the initial parameter
## guessing.
## This can 
guess.constrained.start <- function(func, x.init, warn=TRUE) {
  f.orig <- environment(func)$f
  names.orig <- argnames(f.orig)
  names.cons <- argnames(func)
  n.orig <- length(names.orig)
  n.cons <- length(names.cons)
  arg.idx <- match(names.cons, names.orig)

  if ( length(x.init) == n.cons ) {
    x.init
  } else if  ( length(x.init) == n.orig && !any(is.na(arg.idx)) ) {
    if ( warn )
      warning("Guessing parameters while constraining model - may do badly")
    x.init <- x.init[arg.idx]
  } else {
    stop("Could not guess reduced parameter set from those given")
  }

  x.init
}

logLik.fit.mle <- function(object, ...) {
  ll <- object$lnLik
  attr(ll, "df") <- length(object$par)
  class(ll) <- "logLik"
  ll
}

coef.fit.mle <- function(object, full=FALSE, extra=FALSE, ...) {
  func <- object$func
  if ( full && inherits(func, "constrained") )
    if ( extra && !is.null(extra.v <- attr(func, "extra")) )
      c(object$par[extra.v], func(object$par, pars.only=TRUE))
    else
      func(object$par, pars.only=TRUE)
  else
    object$par
}

extractAIC.fit.mle <- function(fit, scale, k=2, ...)
  c(length(coef(fit)), AIC(fit))


## Code based on MASS:::anova.negbin and ape:::anova.ace
anova.fit.mle <- function(object, ..., sequential=FALSE) {
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
  df <- sapply(ll, attr, "df")

  if ( sequential ) {
    ddf <- c(NA, diff(df))
    if ( any(ddf[-1] < 1) )
      stop("Models are not ordered correctly for sequential anova")
    chisq <- c(NA, 2*diff(ll.val))
    if ( any(chisq[-1] < 0 ))
      warning("Impossible chi-square values - convergence failures?")
  } else {
    chisq <- c(NA, abs(2*(ll.val[1] - ll.val[-1])))
    ddf <- c(NA, abs(df[1] - df[-1]))
  }
  
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
