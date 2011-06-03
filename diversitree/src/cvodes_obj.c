#include <R.h>

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

RCvodesObj* make_cvodes(int neq, int np, CVRhsFn rhs, double rtol,
			double *atol) {
  RCvodesObj *obj;
  void *cvode_mem;
  int i, flag;
  N_Vector y0, abstol;
  UserData *data;
  realtype reltol;
  int maxsteps = 10000; /* TODO: Configurable */
  
  /*** 3: Set vector of initial values ***/
  /* Create serial vector of length neq for the initial conditions.
     This is just a dummy holding pen for now, and we set the values
     all equal to zero.  This should not hurt, but it might be better
     to accept a dummy vector of plausible initial conditions? */
  y0 = N_VNew_Serial(neq);
  if (cvodes_check_flag((void *)y0, "N_VNew_Serial", 0)) 
    return(NULL);
  for( i = 0; i < neq; i++ )
    Ith(y0,i+1) = 0.0;

  /*** 4: Create CVODES object */
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton
   iteration */
  /* TODO: Other options here:
     CV_BDF     specifies the Backward Differentiation Formula (or CV_ADAMS)
     CV_NEWTON  specifies a Newton iteration (or CV_FUNCTIONAL)
     A pointer to the integrator problem memory is returned and stored
     in cvode_mem.     
   */
  /* For stiff problems, this is recommended */
  /* cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);*/
  /* But I am assuming nonstiff and going with */
  cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  if (cvodes_check_flag((void *)cvode_mem, "CVodeCreate", 0)) 
    return(NULL);

  /*** 5: Initialize CVODES solver ***/
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, rhs, RCONST(0.0), y0);
  if (cvodes_check_flag(&flag, "CVodeInit", 1)) 
    return(NULL);

  /*** 6: Specify integration tolerances ***/
  /* Set the scalar relative tolerance */
  reltol = RCONST(rtol);

  /* Set the vector absolute tolerance */
  abstol = N_VNew_Serial(neq); 
  if (cvodes_check_flag((void *)abstol, "N_VNew_Serial", 0)) 
    return(NULL);
  for ( i = 0; i < neq; i++ )
    Ith(abstol,i+1) = atol[i];

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (cvodes_check_flag(&flag, "CVodeSVtolerances", 1)) 
    return(NULL);

  /*** 7: Set optional inputs ***/
  /* This would be where we would use CVodeSet to change any optional
     inputs.  Useful ones, used in Rsundials, are:
     
     CVodeSetUserData(cvode_mem, userdata)
     CVodeSetMaxStep(cvode_mem, maxstepsize)
  */
  CVodeSetMaxNumSteps(cvode_mem, maxsteps);

  /* Here, create some space for the user data, and initialise a
     zeroed vector of parameters */
  data = calloc(1, sizeof(UserData*));
  if (cvodes_check_flag((void *)data, "calloc", 2)) 
    return(NULL);
  data->np  = np;
  data->neq = neq;
  data->p = calloc(np, sizeof(realtype));
  for ( i = 0; i < np; i++ )
    data->p[i] = RCONST(0.0);
  CVodeSetUserData(cvode_mem, data);

  /*** 8: Attach linear solver module ***/
  /* Other options include
     CVLapackDense(cvode_mem, N)
     CVBand(cvode_mem, N, mupper, mlower)
     CVLapackBand(cvode_mem, N, mupper, mlower)
     CVDiag(cvode_mem);
     CVSpgmr(cvode_mem, pretype, maxl)
     CVSptfqmr(cvode_mem, pretype, maxl)
   */
  flag = CVDense(cvode_mem, neq);
  if ( cvodes_check_flag(&flag, "CVDense", 1) )
    return(NULL);

  /*** 9: Set linear solver optional inputs */

  /*** 10: Specify rootfinding problem ***/

  /* At this point, all the setting up we can do is done - let's
     package it up and return */
  obj = calloc(1, sizeof(RCvodesObj));
  obj->neq  = neq;
  obj->np   = np;
  obj->y    = y0;
  obj->data = data;
  obj->p    = data->p;
  obj->cvode_mem = cvode_mem;
  obj->reltol = reltol;
  obj->abstol = abstol;

  return obj;
}

void cvodes_cleanup(RCvodesObj *obj) {
  /*** 13: Deallocate memory for solution vector ***/
  N_VDestroy_Serial(obj->y);
  N_VDestroy_Serial(obj->abstol);
  free(obj->p);
  free(obj->data);

  /*** 14: Free solver memory ***/
  CVodeFree(&(obj->cvode_mem));
  free(obj);
}

void cvodes_set_pars(RCvodesObj* obj, double *pars) {
  int i;
  realtype *p = obj->p;
  for ( i = 0; i < obj->np; i++ )
    p[i] = pars[i];
}

int cvodes_run(RCvodesObj *obj, double *y0, double *times, int nt, 
	       double *ret) {
  int i, flag;
  int neq = obj->neq;
  realtype t, tout;
  void *cvode_mem = obj->cvode_mem;
  N_Vector y = obj->y;
  double *ydat = NV_DATA_S(y);

  /* Initialise y */
  for( i = 0; i < neq; i++ )
    ret[i] = ydat[i] = y0[i];

  CVodeReInit(cvode_mem, times[0], y);

  /*** 11: Advance solution in time ***/
  for ( i = 1; i < nt; i++ ) {
    tout = times[i];
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if ( flag != CV_TOO_CLOSE && cvodes_check_flag(&flag, "CVode", 1) )
      return -1;

    memcpy(ret + i*neq, ydat, neq * sizeof(double));
  }
  return 0;
}

int cvodes_check_flag(void *flagvalue, char *funcname, int opt) {
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

#endif
