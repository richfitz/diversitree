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
#include "cvodes_obj_fwd.h"

static void cvodes_fwd_finalize(SEXP extPtr);

SEXP r_make_cvodes_fwd(SEXP r_neq, SEXP r_np, 
		       SEXP r_rhs, SEXP r_sens1,
		       SEXP r_rtol, SEXP r_atol) {
  RCvodesObj *obj;
  SEXP extPtr;
  CVRhsFn rhs = (CVRhsFn) R_ExternalPtrAddr(r_rhs);
  CVSensRhs1Fn sens1 = (CVSensRhs1Fn) R_ExternalPtrAddr(r_sens1);
  int neq = INTEGER(r_neq)[0];

  if ( LENGTH(r_atol) != neq ) 
    error("Incorrect length error");

  obj = make_cvodes_fwd(neq, INTEGER(r_np)[0], rhs, sens1,
			REAL(r_rtol)[0], REAL(r_atol));
  if ( obj == NULL )
    error("Error building CVODES[fwd] integrator");

  PROTECT(extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(extPtr, cvodes_fwd_finalize);
  UNPROTECT(1);
  return extPtr;
}

static void cvodes_fwd_finalize(SEXP extPtr) {
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);
  Rprintf("Cleaning up (sens)\n");
  cvodes_fwd_cleanup(obj);
}

SEXP r_cvodes_fwd_run(SEXP extPtr, SEXP r_y0, SEXP r_ys0, 
		      SEXP r_times) {
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);
  int neq = obj->neq, nS = obj->np, nt = LENGTH(r_times), flag;
  SEXP r_ret, r_retY, r_retS;

  PROTECT(r_ret  = allocVector(VECSXP, 2));
  PROTECT(r_retY = allocMatrix(REALSXP, nt, neq));
  PROTECT(r_retS = allocMatrix(REALSXP, nt, neq*nS));
  SET_VECTOR_ELT(r_ret, 0, r_retY);
  SET_VECTOR_ELT(r_ret, 1, r_retS);

  flag = cvodes_fwd_run(obj, REAL(r_y0), REAL(r_ys0),
			REAL(r_times), LENGTH(r_times),
			REAL(r_retY), REAL(r_retS));

  if ( flag != 0 )
    error("Solver failed...");

  UNPROTECT(3);

  return r_ret;
}
  
#endif
