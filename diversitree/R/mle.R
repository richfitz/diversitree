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

do.mle.search <- function(func, x.init, method, fail.value=-Inf,
                          control=list(), lower=-Inf, upper=Inf,
                          hessian=FALSE, verbose=0, ...) {
  method <- match.arg(method, c("optim", "subplex", "nlminb", "nlm",
                                "minqa", "tgp", "optimize", "int1d",
                                "mixed", "subplexR"))

  control$y.init <- y.init <- func(x.init, ...)
  if ( inherits(y.init, "try-error") )
    stop("The above error occured when testing the starting position")
  if ( !is.finite(y.init) )
    stop("Starting point must have finite probability")

  if ( is.null(control$fail.penalty) )
    control$fail.penalty <- 1000

  ## Can handle -Inf: subplex, nlminb, int1d, mixed, tgp
  ## Require finite values: optim, nlm, minqa, optimize
  if ( is.na(fail.value) ) {
    if ( method %in% c("optim", "nlm", "minqa", "optimize") )
      fail.value <- y.init - control$fail.penalty
    else if ( method %in% c("subplex", "nlminb", "int1d", "mixed",
                            "tgp", "subplexR") )
      fail.value <- -Inf
  }

  ## Protect the function, and combine in the function arguments:
  func2 <- protect(function(x) func(x, ...), fail.value)
  
  ## Add in verbosity, if needed:
  control$verbose <- verbose > 0
  if ( control$verbose )
    func2 <- func2 <- big.brother(func2, verbose)

  ## Remember the names of the input vector, or try and get it from
  ## the supplied function.
  if ( is.null(names(x.init)) && !is.null(x.init) ) {
    names.v <- try(argnames(func), silent=TRUE)
    if ( inherits(names.v, "try-error") ) {
      names.v <- NULL
    } else {
      if ( length(names.v) != length(x.init) )
        stop("Invalid parameter names length: expected ",
             length(names.v), ", got ", length(x.init))
      names(x.init) <- names.v
    }
  } else {
    names.v <- names(x.init)
  }

  mle.search <- get(sprintf("do.mle.search.%s", method))
  ans <- mle.search(func2, x.init, control, lower, upper)

  if ( verbose )
    cat("\n")

  if ( hessian ) {
    if ( !require(numDeriv) )
      stop("The package numDeriv is required to compute the hessian")
    ans$hessian <- hessian(func, ans$par)
  }

  names(ans$par) <- names.v
  ans$func <- func # drop or make option - this can be annoying
  ans$method <- method
  class(ans) <- "fit.mle"
  ans
}

do.mle.search.optim <- function(func, x.init, control, lower, upper) {
  control <- modifyList(list(fnscale=-1,
                             ndeps=rep(1e-5, length(x.init)),
                             optim.method="L-BFGS-B"), control)

  optim.method <- control$optim.method
  control.optim <- control[c("fnscale", "ndeps")]

  ans <- optim(x.init, func, method=optim.method,
               control=control.optim, lower=lower, upper=upper)
  names(ans)[names(ans) == "value"] <- "lnLik"  
  ans$optim.method <- optim.method
  
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (optim): ",
            tolower(ans$message))

  ans
}

do.mle.search.subplex <- function(func, x.init, control, lower, upper) {
  ## By default, lower tolerance-- more likely to be met
  control <- modifyList(list(reltol=.Machine$double.eps^0.25,
                             parscale=rep(.1, length(x.init))),
                        control)

  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)

  ans <- subplex(x.init, func2, control)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"

  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (subplex): ",
            tolower(ans$message))
  
  ans
}

do.mle.search.nlminb <- function(func, x.init, control, lower, upper) {
  control.nlminb.ok <- c("eval.max", "iter.max", "trace", "abs.tol",
                         "rel.tol", "x.tol", "step.min")
  control.nlminb <- control[names(control) %in% control.nlminb.ok]
  ans <- nlminb(x.init, invert(func), control=control.nlminb,
                lower=lower, upper=upper)
  names(ans)[names(ans) == "objective"] <- "lnLik"
  names(ans)[names(ans) == "evaluations"] <- "counts"
  ans$lnLik <- -ans$lnLik
  ans <- ans[c("par", "lnLik", "counts", "convergence", "message",
               "iterations")]
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (nlminb): ",
            tolower(ans$message))
  ans
}

do.mle.search.nlm <- function(func, x.init, control, lower, upper) {
  nlm.defaults <-
    list(typsize=rep(1, length(x.init)), print.level=0, ndigit=12,
         gradtol=1e-06, steptol=1e-06, iterlim=100,
         check.analyticals=TRUE)
  control <- modifyList(nlm.defaults, control)
  if ( is.null(control$stepmax) )
    control$stepmax <- max(1000*sqrt(sum((x.init/control$typsize)^2)), 1000)

  ans <- nlm(invert(func), x.init, typsize=control$typsize,
             print.level=control$print.level, ndigit=control$ndigit,
             gradtol=control$gradtol, stepmax=control$stepmax,
             steptol=control$steptol, iterlim=control$iterlim,
             check.analyticals=control$check.analyticals)

  names(ans)[names(ans) == "estimate"] <- "par"
  names(ans)[names(ans) == "minimum"] <- "lnLik"
  names(ans)[names(ans) == "iterations"] <- "counts"
  ans$lnLik <- -ans$lnLik 
  ans <- ans[c("par", "lnLik", "counts", "code", "gradient")]
  
  if ( ans$code > 2 )
    warning("Convergence problems in find.mle: code = ",
            ans$code, " (see ?nlm for details)")  

  ans
}

do.mle.search.minqa <- function(func, x.init, control, lower, upper) {
  if ( !require(minqa) )
    stop("This method requires the minqa package")
  
  control <- modifyList(list(minqa.method="newuoa"), control)
  minqa.method <- match.arg(control$minqa.method,
                            c("bobyqa", "newuoa", "uobyqa"))

  control.minqa.ok <- c("npt", "rhobeg", "rhoend", "iprint", "rho",
                        "maxfun")
  control.minqa <- control[names(control) %in% control.minqa.ok]

  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)

  opt <- get(minqa.method, "package:minqa")
  if ( minqa.method == "bobyqa" )
    ans <- opt(x.init, func2, lower=lower, upper=upper,
               control=control.minqa)
  else
    ans <- opt(x.init, func2, control=control.minqa)
  
  list(par=ans$par, lnLik=ans$fval, counts=ans$feval,
       minqa.method=minqa.method)
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

