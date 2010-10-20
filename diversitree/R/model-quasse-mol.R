## Some jiggling around will be needed, as the entire "padding" thing
## goes away with the MOL approach.

## There appears to be an issue with the "high" resolution version:
##   Error in lsodes(y = y, times = times, func = func, parms,
##     sparsetype = "1D",  : 
##   REAL() can only be applied to a 'numeric', not a 'list'
## Which is a bit cryptic, but possibly means that the memory
## allocation failed?
make.branches.quasse.mol <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  tc <- control$tc
  
  f.hi <- make.pde.quasse.mol(nx*r, dx/r, 2L)
  f.lo <- make.pde.quasse.mol(nx,   dx,   2L)
  make.branches.quasse(f.hi, f.lo, control)
}

make.pde.quasse.mol <- function(nx, dx, nd) {
  function(y, len, pars, t0) {
    ans <- ode.1D(y, c(t0, t0+len), "derivs_quasse_mol", pars,
                  initfunc="initmod_quasse_mol",
                  nspec=nd, dimens=nx,
                  method="lsodes", dllname="diversitree")[-1,-1]
    matrix(ans, nx, nd)
  }
}

## Stripped down integrator from deSolve 1.7:
## ode.1D.2 <- function(y, times, func, parms, nspec = NULL, dimens = NULL,
##                      ...) {
##   N <- length(y)
##   nspec = N/dimens
##   if (is.null(dimens) )
##     dimens <- N/nspec
##   if ( N%%nspec != 0 )
##     stop ("cannot run ode.1D: nspec is not an integer fraction of number of state variables")

##   out <- lsodes2(y=y, times=times, func=func, parms, sparsetype="1D",
##                  nnz=c(nspec,dimens), ...)

##   attr (out,"dimens") <- dimens
##   attr (out,"nspec") <- nspec

##   out
## }

## lsodes2 <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
##                     jacvec=NULL, sparsetype="sparseint", nnz = NULL,
##                     inz = NULL, verbose=FALSE, tcrit = NULL, hmin=0,
##                     hmax=NULL, hini=0, ynames=TRUE, maxord=NULL,
##                     maxsteps=5000, lrw=NULL, liw=NULL, dllname=NULL,
##                     initfunc=dllname, initpar=parms, rpar=NULL,
##                     ipar=NULL, nout=0, outnames=NULL, forcings=NULL,
##                     initforc = NULL, fcontrol=NULL, events=NULL, lags
##                     = NULL, ...)  {
##   hmax <- max(abs(diff(times)))
##   n <- length(y)
##   maxord <- 5

##   # imp = method flag as used in lsodes
##   imp <- 22   # inz supplied,jac not supplied  

##   ## Special-purpose sparsity structures: 1-D and 2-D
##   ## reaction-transport problems

##   ## Typically these applications are called via ode.1D, ode.2D and
##   ## ode.3D Here the sparsity is specified in the C-code; this needs
##   ## extra input: the number of components *nspec* and the
##   ## dimensionality of the problem (number of boxes in each
##   ## direction).  This information is passed by ode.1D, ode.2D and
##   ## ode.3D in parameter nnz (a vector).  nnz is altered to include
##   ## the number of nonzero elements (element 1).  'Type' contains the
##   ## type of sparsity + nspec + num boxes + cyclicBnd

##   nspec  <- nnz[1]
##   Type   <- c(2,nnz)    #type=2
##   nnz    <- n*(2+nspec)-2*nspec

##   ## model and Jacobian function 
##   JacFunc   <- NULL
##   Ynames    <- attr(y, "names")
##   flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
##   Eventfunc <- NULL
##   events <- list()

##   DLL <- deSolve:::checkDLL(func,jacvec,dllname,
##                             initfunc,verbose,nout, outnames, JT=2)

##   ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE=dllname)$address
##   Func    <- getNativeSymbolInfo(func, PACKAGE=dllname)$address
##   Nglobal <- 0
##   Nmtot   <- NULL

##   rho <- NULL
##   ipar <- 0
##   rpar <- 0

##   ## work arrays iwork, rwork
##   # 1. Estimate length of rwork and iwork if not provided via arguments lrw, liw
##   moss  <- imp%/%100         # method to be used to obtain sparsity
##   meth  <- imp%%100%/%10     # basic linear multistep method
##   miter <- imp%%10           # corrector iteration method
##   lenr = 2     # real to integer wordlength ratio (2 due to double precision)

##   if (is.null(lrw)) {         # make a guess of real work space needed
##     lrw <- 20+n*(maxord+1)+3*n +20  #extra 20 to make sure
##     lrw <- lrw + 2*nnz + 2*n + (nnz+10*n)/lenr
##     lrw <- lrw*1.1 # increase to be sure it is enough...
##   }
  
##   if (is.null(liw)) {         # make a guess of integer work space needed
##     liw <- 31+n+nnz +30  # extra 30
##   }

##   # 2. Allocate and set values
##   # only first 20 elements of rwork passed to solver;
##   # other elements will be allocated in C-code
##   # for iwork: only first 30 elements, except when sparsity imposed
##   rwork <- numeric(20)

##   # iwork will contain sparsity structure (ian,jan)
##   # See documentation of DLSODES how this is done
##   ## sparsity not imposed; only 30 element of iwork allocated.
##   iwork <- numeric(30) # TODO - should be integer?

##   # other elements of iwork, rwork
##   iwork[5] <- maxord
##   iwork[6] <- maxsteps
  
##   rwork[5] <- hini
##   rwork[6] <- hmax
##   rwork[7] <- hmin

##   # the task to be performed.
##   itask <- 1

##   ## calling solver
##   storage.mode(y) <- storage.mode(times) <- "double"
##   IN <- 3

##   lags <- list(islag=0L)

##   out <- .Call("call_lsoda",y,times,Func,initpar,
##                rtol, atol, rho, tcrit, JacFunc, ModelInit, Eventfunc,
##                as.integer(verbose), as.integer(itask), as.double(rwork),
##                as.integer(iwork), as.integer(imp),as.integer(Nglobal),
##                as.integer(lrw),as.integer(liw),as.integer(IN),
##                NULL, as.integer(0), as.double (rpar), as.integer(ipar),
##                as.integer(Type),flist, events, lags, PACKAGE="deSolve")
  
  
## ### saving results    
##   out <- deSolve:::saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
##                  iin=c(1,12:20), iout=c(1:3,14,5:9,17))

##   attr(out, "type") <- "lsodes"
##   if (verbose) diagnostics(out)
##   out
## }
