/*
 * The ClaSSE equations, implemented in C
 * By Emma Goldberg <eeg@uic.edu>
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include "util.h"

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_classe(int n, double *pars, double *y, double *ydot, 
                      int jk_array[][2]) {
  /* note: n = num states is called k elsewhere, but k is used as 
   * an index here */

  double *E = y, *D = y + n;
  double *dEdt = ydot, *dDdt = ydot + n;
  int len_lam_i = n * (n + 1) / 2;
  int len_lam = n * len_lam_i;
  double *lambda = pars, *lambda_i;
  double *mu = pars + len_lam, *Q = pars + len_lam + n;
  double Ei, Di;
  int i, j, k, m;

  for (i=0; i<n; i++)
  {
    Ei = E[i];
    Di = D[i];

    /* extinction (start ydot) */
    dEdt[i] = mu[i] * (1 - Ei);
    dDdt[i] = -mu[i] * Di;

    /* speciation (continue ydot) */
    lambda_i = lambda + i * len_lam_i;
    for (m=0; m<len_lam_i; m++)
    {
      j = jk_array[m][0];
      k = jk_array[m][1];
      dEdt[i] += lambda_i[m] * (-Ei + E[j] * E[k]);
      dDdt[i] += lambda_i[m] * (-Di + D[j] * E[k] + D[k] * E[j]);
    }
  }

  /* transitions (complete ydot by adding Q y to it) */
  /* RGF TODO: Replace with mult_mv() from util-matrix.c */
  /* mult_mv(k, Q, y, 1.0, ydot); */
  do_gemm2(Q, n, n, y, n, 2, ydot);
}

/* Plain ClaSSE */
/* deSolve / LSODA */
static double *parms_classe;

void initmod_classe(void (* odeparms)(int *, double *)) {
  /* parameter length checking will be done in R */
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_classe = REAL(get_deSolve_gparms());
} 

/* compute indices for the lambda's */
void fill_jk_array(int jk_array[][2], int n) {
  int j, k, m;
  m = 0;
  for (j=0; j<n; j++)
  {
    for (k=j; k<n; k++)
    {
      jk_array[m][0] = j;
      jk_array[m][1] = k;
      m++;
    }
  }
}

void derivs_classe(int *neq, double *t, double *y, double *ydot,
                     double *yout, int *ip) {
  int n = *neq / 2;

  /* saves a little time to pre-compute indices outside of do_derivs_classe;
   * but do it outside of here, too? */
  int len_lam_i = n * (n + 1) / 2;
  int jk_array[500][2]; /* was: int jk_array[len_lam_i][2]; */
  fill_jk_array(jk_array, n);

  do_derivs_classe(n, parms_classe, y, ydot, jk_array);
}

/* CVODES */
int derivs_classe_cvode(realtype t, N_Vector y, N_Vector ydot,
                        void *user_data) {
  const UserData *data = (UserData*) user_data;

  int n = data->neq/2;

  int len_lam_i = n * (n + 1) / 2;
  int jk_array[500][2]; /* was: int jk_array[len_lam_i][2]; */
  fill_jk_array(jk_array, n);

  do_derivs_classe(data->neq/2,
                   data->p,
                   NV_DATA_S(y),
                   NV_DATA_S(ydot), jk_array);
  return 0;
}

void initial_conditions_classe(int neq, double *vars_l, double *vars_r,
                              double *pars, double t, double *vars_out) {
  /* note: n = num states is called k elsewhere, but k is used as 
   * an index here */
  const int n = neq/2;
  double *lambda = pars, *lambda_i;
  int i, j, k, m;

  int len_lam_i = n * (n + 1) / 2;
  int jk_array[500][2]; /* was: int jk_array[len_lam_i][2]; */
  fill_jk_array(jk_array, n);

  /* E: */
  memcpy(vars_out, vars_l, n * sizeof(double));

  /* D: */
  for(i=0; i<n; i++)
  {
    vars_out[i+n] = 0;
    lambda_i = lambda + i * len_lam_i;
    for (m=0; m<len_lam_i; m++)
    {
      j = jk_array[m][0] + n;
      k = jk_array[m][1] + n;
      vars_out[i+n] += 0.5 * lambda_i[m] * (vars_l[j] * vars_r[k] 
                       + vars_l[k] * vars_r[j]);
    }
  }
}


/* no time-dependence yet */
