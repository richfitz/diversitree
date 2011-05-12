/* TODO: I need to be careful where realtype is not double, i.e.
     ifndef SUNDIALS_DOUBLE_PRECISION // trouble
*/
#include "cvodes/include/user_data.h"

typedef struct {
  int neq;
  int np;

  N_Vector y;

  UserData *data;
  realtype *p;

  void *cvode_mem;

  realtype reltol;
  N_Vector abstol;

  /* ----------------------------
     Forward sensitivity analysis
     ---------------------------- */
  int nS;
  N_Vector *yS;
  int ism;
  booleantype err_con;

} RCvodesObj;

/* 
   Description from CVODES/examples/cvodes/serial/cvsRoberts_dns.c

   These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

RCvodesObj* make_cvodes(int neq, int np, CVRhsFn rhs, double rtol,
			double *atol);
void cvodes_cleanup(RCvodesObj *obj);
void cvodes_set_pars(RCvodesObj* obj, double *pars);
int cvodes_run(RCvodesObj *obj, double *y0, double *times, int nt, 
		double *ret);
void cvodes_set_pars(RCvodesObj* obj, double *pars);

int cvodes_check_flag(void *flagvalue, char *funcname, int opt);
