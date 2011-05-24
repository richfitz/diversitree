/* Second shot at an object-based integration solution */
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

#include "cvodes_obj.h"
#include "cvodes_obj_fwd.h"

RCvodesObj* make_cvodes_fwd(int neq, int np, CVRhsFn rhs, 
			    CVSensRhs1Fn sens1, double rtol,
			    double *atol) {
  RCvodesObj *obj = make_cvodes(neq, np, rhs, rtol, atol);
  N_Vector *yS;
  int i, flag;
  void *cvode_mem = obj->cvode_mem;
  int nS = np; /* TODO: Hard coded for now */

  /*
   * In the CV SIMULTANEOUS approach, the state and sensitivity
     variables are corrected at the same time. If CV NEWTON was
     selected as the non- linear system solution method, this amounts
     to performing a modified Newton iteration on the combined
     nonlinear system;

   * In the CV STAGGERED approach, the correction step for the
     sensitivity variables takes place at the same time for all
     sensitivity equations, but only after the correction of the state
     variables has converged and the state variables have passed the
     local error test;

   * In the CV STAGGERED1 approach, all corrections are done
     sequentially, first for the state variables and then for the
     sensitivity variables, one parameter at a time. If the
     sensitivity variables are not included in the error control, this
     approach is equivalent to CV STAGGERED. Note that the CV
     STAGGERED1 approach can be used only if the user-provided
     sensitivity right-hand side function is of type CVSensRhs1Fn (see
     S5.3).
  */
  int ism = CV_SIMULTANEOUS;

  /* Explain */
  booleantype err_con = FALSE;

  yS = NULL;

  /*** 10: Define the sensitivity problem ***/
  /* Problem parameters (optional):

     Because I am going to provide (for now at least) a RHS, I don't
     have to muck around setting a parameter vector data->p.  However,
     I would like to try this at some point. */

  /* Parameter list (optional):

     Array of nS integers to specify the parameters 'p' with respect
     to which solution sensitivities are required.  Otherwise
     p[0:(nS-1)] are used, which is fine for now. */

  /* Parameter scaling values (optional)

     This is needed if internal difference-quotient function is used.
     |p| for |p| > 0 is good, but not needed here */
  
  /*** 11: Set sensitivity initial conditions ***/
  /* Here, the 'y' vector that points at the initial conditions
     serves to provide the N_Vector type for cloning.  These are then
     zeroed below.  However, it would be best if these initial
     conditions were specifiable.  

     This is basically a matrix of derivatives though, so setting the
     initial conditions will be slightly tricky.  Ideally, we will
     have things in appropriate data structures already. */
  yS = N_VCloneVectorArray_Serial(nS, obj->y);
  if ( cvodes_check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0) )
    return(NULL);
  for ( i = 0; i < nS; i++ )
    N_VConst(RCONST(0.0), yS[i]);

  /*** 12: Activate sensitivity calculations ***/
  /* TODO: fS=NULL, at least within CVodeSensInit() activates internal
     difference quotient RHS routine */
  flag = CVodeSensInit1(cvode_mem, nS, ism, sens1, yS);
  if ( cvodes_check_flag(&flag, "CVodeSensInit", 1) )
    return(NULL);

  /*** 13: Set sensitivity tolerances ***/
  /* CVodeSens{XX}tolerances
     where {XX} is SS, SV or EE
     SS = Scalar relative and absolute tolerances
     SV = Scalar relative and vector absolute
     EE = Estimate tolerances based on scaling factor pbar and
          tolerances for state variables.
  */
  flag = CVodeSensEEtolerances(cvode_mem);
  if(cvodes_check_flag(&flag, "CVodeSensEEtolerances", 1)) 
    return(NULL);

  /*** 14: Set sensitivity analysis optional inputs ***/
  /* Activate, or not, error control.  If err_con is FALSE, then the
     sensitivity calculations are not included in the error control
     mechanisms */
  flag = CVodeSetSensErrCon(cvode_mem, err_con);
  if (cvodes_check_flag(&flag, "CVodeSetSensErrCon", 1)) 
    return(NULL);

  /* Here, the three NULL values are
     p: 'realtype' vector of parameters user_data->p
     pbar: 'realtype' vector of scaling values >0.9
     plist: 'int' vector of indices to check
  */
  flag = CVodeSetSensParams(cvode_mem, NULL, NULL, NULL);
  if (cvodes_check_flag(&flag, "CVodeSetSensParams", 1)) 
    return(NULL);

  /*** 15(10): Specify rootfinding problem ***/
  obj->nS  = nS;
  obj->yS  = yS;
  obj->ism = ism; /* no need for ism or err_con? */
  obj->err_con = err_con;

  return obj;
}

void cvodes_fwd_cleanup(RCvodesObj *obj) {
  /*** 19: Deallocate memory for sensitivity vectors ***/
  N_VDestroyVectorArray_Serial(obj->yS, obj->nS);

  cvodes_cleanup(obj);
}

int cvodes_fwd_run(RCvodesObj *obj, double *y0, double *ys0, 
		   double *times, int nt,
		   double *retY, double *retS) {
  int i, j, flag;
  int neq = obj->neq, nS=obj->nS;
  realtype t, tout;
  void *cvode_mem = obj->cvode_mem;
  N_Vector y = obj->y, *yS = obj->yS;

  /* Initialise y */
  for ( i = 0; i < neq; i++ )
    retY[nt*i] = Ith(y,i+1) = y0[i];
  copy_yS_in(neq, nS, nt, yS, ys0);

  CVodeReInit(cvode_mem, times[0], y);
  CVodeSensReInit(cvode_mem, obj->ism, yS);

  /* Copy into the return matrix */
  copy_yS_out(neq, nS, nt, yS, retS);
  /* TODO: more efficient direct r_y0 -> retS? */

  for ( i = 1; i < nt; i++ ) {
    tout = times[i];
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if (cvodes_check_flag(&flag, "CVode", 1))
      return(-1);
    
    /* Copy solution into array */
    for ( j = 0; j < neq; j++ )
      retY[i + nt*j] = Ith(y, j+1);

    /*** 17 Extract sensitivity solution ***/
    flag = CVodeGetSens(cvode_mem, &t, yS);
    if (cvodes_check_flag(&flag, "CVodeGetSens", 1))
      return(-2);
    copy_yS_out(neq, nS, nt, yS, retS + i);
  }

  return 0;
}

void copy_yS_out(int neq, int nS, int nt, N_Vector *yS, double *out) {
  int i, j, k = 0;
  realtype *ySi;
  for ( i = 0; i < neq; i++ ) {
    ySi = NV_DATA_S(yS[i]);
    for ( j = 0; j < nS; j++, k += nt )
      out[k] = ySi[j];
  }
}
void copy_yS_in(int neq, int nS, int nt, N_Vector *yS, double *in) {
  int i, j, k = 0;
  realtype *ySi;
  for ( i = 0; i < neq; i++ ) {
    ySi = NV_DATA_S(yS[i]);
    for ( j = 0; j < nS; j++, k++ )
      ySi[j] = in[k];
  }
}

#endif
