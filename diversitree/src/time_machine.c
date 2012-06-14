#include <R.h>
#include <Rinternals.h>
#include "time_machine.h"
#include "util.h"

/* 1. Here are possible types of functions, and their constants */
#define T_CONSTANT 0
#define T_LINEAR   1
#define T_STEPF    2
#define T_SIGMOID  3
#define T_SPLINE   4

double t_linear(double t, double *p) {
  return p[0] + t * p[1];
}
double t_stepf(double t, double *p) {
  return t <= p[2] ? p[0] : p[1];
}
double t_sigmoid(double t, double *p) {
  //y0    y1   - p0             r       tmid
  return p[0] + (p[1] - p[0])/(1 + exp(p[4] * (p[3] - t)));
}

/* A little helper function to order these in R so I don't have to
   rely on remembering */
SEXP r_get_time_machine_types() {
  SEXP ret;
  PROTECT(ret = allocVector(STRSXP, 5));
  SET_STRING_ELT(ret, T_CONSTANT, mkChar("constant.t"));
  SET_STRING_ELT(ret, T_LINEAR,   mkChar("linear.t"));
  SET_STRING_ELT(ret, T_STEPF,    mkChar("stepf.t"));
  SET_STRING_ELT(ret, T_SIGMOID,  mkChar("sigmoid.t"));
  SET_STRING_ELT(ret, T_SPLINE,   mkChar("spline.t"));
  UNPROTECT(1);
  return ret;
}

/* 2. The time machine itself */

/* Make/cleanup function */
static void dt_time_machine_finalize(SEXP extPtr);
SEXP r_time_machine(SEXP obj) {
  SEXP extPtr;
  /* Extract stored R objects to build the time machine */
  SEXP
    types           = getListElement(obj, "types"),
    start           = getListElement(obj, "start");
  int np_in = INTEGER(getListElement(obj, "np.in"))[0],
    np_out  = INTEGER(getListElement(obj, "np.out"))[0],
    nf = LENGTH(types);
  dt_time_machine *ret;

  ret = (dt_time_machine *)Calloc(1, dt_time_machine);
  ret->np_in  = np_in;
  ret->p_in   = Calloc(np_in, double);
  ret->np_out = np_out;
  ret->p_out   = Calloc(np_out, double);
  ret->nf      = nf;

  ret->types  = Calloc(nf, int);
  ret->start  = Calloc(nf, int);
  memcpy(ret->types, INTEGER(types), nf*sizeof(int));
  memcpy(ret->start, INTEGER(start), nf*sizeof(int));

  /* And, done */
  extPtr = R_MakeExternalPtr(ret, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, dt_time_machine_finalize);
  return extPtr;
}

void cleanup_time_machine(dt_time_machine *obj) {
  Free(obj->p_in);
  Free(obj->p_out);
  Free(obj->types);
  Free(obj->start);
  Free(obj);
}

static void dt_time_machine_finalize(SEXP extPtr) {
  dt_time_machine *obj = (dt_time_machine*)R_ExternalPtrAddr(extPtr);
  cleanup_time_machine(obj);
}

/* When initialising a time machine, we set all constant arguments to
   their current value 
   
   It's possible that we can actually do better here; where (say) the
   function is linear and the slope is exactly zero, we can skip this.
   There are a limited enough range of options that this might be
   simple enough to implement and fast enough to warrent implementing.
   
   Furthermore looping over an array of indices
   (which(type!=constant)) will possibly be faster when most things
   aren't variable.
   
   ** TODO **
   check to see if pars is the same as previous
   memcmp(pars, obj->p_in)
   and if so just skip updating?
*/
void init_time_machine(dt_time_machine *obj, double *pars) {
  double *p_out = obj->p_out;
  const int np_in = obj->np_in, nf = obj->nf;
  const int *types = obj->types, *start = obj->start;
  int i;

  /* First, store all parameters */
  memcpy(obj->p_in, pars, np_in*sizeof(double));

  /* Then copy in all constant values */
  for ( i = 0; i < nf; i++ )
    if ( types[i] == T_CONSTANT )
      p_out[i] = pars[start[i]];
}

/* Compute values at a time 't' */
void run_time_machine(dt_time_machine *obj, double t) {
  double *p_in = obj->p_in, *p_out = obj->p_out;
  const int nf = obj->nf, *types = obj->types, *start = obj->start;
  int i;
  for ( i = 0; i < nf; i++ )
    switch(types[i]) {
    case T_LINEAR: 
      p_out[i] = t_linear(t, p_in + start[i]);
      break;
    case T_STEPF:
      p_out[i] = t_stepf(t, p_in + start[i]);
      break;
    case T_SIGMOID:
      p_out[i] = t_sigmoid(t, p_in + start[i]);
      break;
    }
}

/* R interface to the above functions */
SEXP r_init_time_machine(SEXP extPtr, SEXP pars) {
  init_time_machine((dt_time_machine*)R_ExternalPtrAddr(extPtr),
		    REAL(pars));
  return R_NilValue;
}

SEXP r_run_time_machine(SEXP extPtr, SEXP t) {
  dt_time_machine *obj = (dt_time_machine*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  run_time_machine(obj, REAL(t)[0]);
  PROTECT(ret = allocVector(REALSXP, obj->np_out));
  memcpy(REAL(ret), obj->p_out, obj->np_out*sizeof(double));
  UNPROTECT(1);
  return ret;
}

/* Not sure if this one is ever useful */
SEXP r_get_time_machine_pars(SEXP extPtr) {
  dt_time_machine *obj = (dt_time_machine*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, obj->np_out));
  memcpy(REAL(ret), obj->p_out, obj->np_out*sizeof(double));
  UNPROTECT(1);
  return ret;
}
