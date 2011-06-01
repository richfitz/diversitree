/* All of the R interface code goes here */
#include <R.h>
#include <Rinternals.h>

#include "config.h"
#ifdef WITH_CVODES

/* Sundials Header Files */
#include "cvodes/include/cvodes/cvodes.h"
#include "cvodes/include/cvodes/cvodes_dense.h"
#include "cvodes/include/nvector/nvector_serial.h"
#include "cvodes/include/sundials/sundials_dense.h"
#include "cvodes/include/sundials/sundials_types.h"

/* My includes */
#include "cvodes_obj.h"

static void cvodes_finalize(SEXP extPtr);
typedef int rhs_func(realtype t, N_Vector y, N_Vector ydot, void *f_data);

SEXP r_make_cvodes(SEXP r_neq, SEXP r_np, SEXP r_rhs, 
		   SEXP r_rtol, SEXP r_atol) {
  RCvodesObj *obj;
  SEXP extPtr;

  CVRhsFn rhs = (CVRhsFn) R_ExternalPtrAddr(r_rhs);
  int neq = INTEGER(r_neq)[0];

  if ( LENGTH(r_atol) != neq ) 
    error("Incorrect length error");
  
  obj = make_cvodes(neq, INTEGER(r_np)[0], rhs, 
		    REAL(r_rtol)[0], REAL(r_atol));
  if ( obj == NULL )
    error("Error building CVODES integrator");

  PROTECT(extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(extPtr, cvodes_finalize);
  UNPROTECT(1);

  return extPtr;
}

static void cvodes_finalize(SEXP extPtr) {
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);
  cvodes_cleanup(obj);
}

SEXP r_cvodes_set_pars(SEXP extPtr, SEXP r_pars) {
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);

  if ( LENGTH(r_pars) != obj->np )
    error("Incorrect length parameters");

  cvodes_set_pars(obj, REAL(r_pars));

  return(R_NilValue);
}

SEXP r_cvodes_run(SEXP extPtr, SEXP r_y0, SEXP r_times) {
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);
  int nt = LENGTH(r_times), flag;
  SEXP r_ret;

  PROTECT(r_ret = allocMatrix(REALSXP, obj->neq, nt));  

  flag = cvodes_run(obj, REAL(r_y0), REAL(r_times), nt, REAL(r_ret));
  if ( flag != 0 )
    error("Solver failed...");

  UNPROTECT(1);

  return r_ret;
}

/* Debug: return contents of several slots */
SEXP r_test(SEXP extPtr) {
  int i;
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);
  Rprintf("neq = %d, np = %d\n", obj->neq, obj->np);

  Rprintf("Contents of 'y':\n");
  for ( i = 0; i < obj->neq; i++ )
    Rprintf("\t%d: %2.5f\n", i+1, Ith(obj->y, i+1));

  Rprintf("Contents of 'p':\n");
  for ( i = 0; i < obj->np; i++ )
    Rprintf("\tp[%d]: %2.5f\n", i, obj->p[i]);

  Rprintf("Tolerance: %2.10f (rel).  Abs: ", obj->reltol);
  for ( i = 0; i < obj->neq; i++ )
    Rprintf("%2.10f  ", Ith(obj->abstol, i+1));
  Rprintf("\n");

  return R_NilValue;
}

#endif
