/*
 * The GeoSSE equations, implemented in C
 * By Emma Goldberg <eeg@uic.edu>
 */

#include <R.h>

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

void do_derivs_geosse(double *pars, double *y, double *ydot) {
  /* states:
   * 1 = in both regions
   * 2 = region A endemic
   * 3 = region B endemic
   */
  double E_1 = y[0];
  double E_2 = y[1];
  double E_3 = y[2];

  double D_N1 = y[3];
  double D_N2 = y[4];
  double D_N3 = y[5];

  double 
    sA  = pars[0],     /* speciation within region A */
    sB  = pars[1],     /* speciation within region B */
    sAB = pars[2],     /* between-region speciation  */
    xA  = pars[3],     /* extinction from region A   */
    xB  = pars[4],     /* extinction from region B   */
    dA  = pars[5],     /* dispersal from A to B      */
    dB  = pars[6];     /* dispersal from B to A      */

  /*  dE_1 / dt  */
  ydot[0] = -(sA + sB + xA + xB + sAB) * E_1 
            + xA * E_3 + xB * E_2 
      + sA * E_1 * E_2 + sB * E_1 * E_3 + sAB * E_2 * E_3;

  /*  dE_2 / dt  */
  ydot[1] = -(sA + dA + xA) * E_2 
            + xA + dA * E_1 + sA * E_2 * E_2;

  /*  E_3 / dt  */
  ydot[2] = -(sB + dB + xB) * E_3 
            + xB + dB * E_1 + sB * E_3 * E_3;

  /*  dD_N1 / dt  */
  ydot[3] = -(sA + sB + sAB + xA + xB) * D_N1 
            + xA * D_N3 + xB * D_N2 
      + sA * (E_2 * D_N1 + E_1 * D_N2) 
      + sB * (E_3 * D_N1 + E_1 * D_N3)
      + sAB * (E_2 * D_N3 + E_3 * D_N2);

  /*  dD_N2 / dt  */
  ydot[4] = -(sA + dA + xA) * D_N2 
            + dA * D_N1 + 2 * sA * D_N2 * E_2;

  /*  dD_N3 / dt  */
  ydot[5] = -(sB + dB + xB) * D_N3 
            + dB * D_N1 + 2 * sB * D_N3 * E_3;
}

/* deSolve / LSODA */
static double parms_geosse[7];

void initmod_geosse(void (* odeparms)(int *, double *))
{
  int N = 7;
  odeparms(&N, parms_geosse);
}

void derivs_geosse(int *neq, double *t, double *y, double *ydot, 
		   double *yout, int *ip) {
  do_derivs_geosse(parms_geosse, y, ydot);
}

/* CVODES */
int derivs_geosse_cvode(realtype t, N_Vector y, N_Vector ydot,
			void *user_data) {
  do_derivs_geosse(((UserData*) user_data)->p,
		   NV_DATA_S(y), NV_DATA_S(ydot));
  return 0;
}

void initial_conditions_geosse(int neq, double *vars_l, double *vars_r,
			       double *pars, double t, 
			       double *vars_out) {
  /* E.AB, E.A, E.B */
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1];
  vars_out[2] = vars_l[2];

  /* D.CAB (Eq. 2c) */
  vars_out[3] =
    0.5 * ((vars_l[3] * vars_r[4] + vars_l[4] * vars_r[3]) * pars[0] +
	   (vars_l[3] * vars_r[5] + vars_l[5] * vars_r[3]) * pars[1] +
	   (vars_l[4] * vars_r[5] + vars_l[5] * vars_r[4]) * pars[2]);

  /* D.CA, D.CB  (Eq. 2ab) */
  vars_out[4] = vars_l[4] * vars_r[4] * pars[0];
  vars_out[5] = vars_l[5] * vars_r[5] * pars[1];
}


/*** For split ***/

/* Auxilliary (just compute E) */

void initmod_geosse_aux(void (* odeparms)(int *, double *))
{
  int N = 7;
  odeparms(&N, parms_geosse);
}

void do_derivs_geosse_aux(double *pars, double *y, double *ydot)
{
  double E_1 = y[0];
  double E_2 = y[1];
  double E_3 = y[2];

  double sA = pars[0], sB = pars[1], sAB = pars[2], 
         xA = pars[3], xB = pars[4], dA = pars[5], dB = pars[6]; 

  /*  dE_1 / dt  */
  ydot[0] = -(sA + sB + xA + xB + sAB) * E_1 
            + xA * E_3 + xB * E_2 
      + sA * E_1 * E_2 + sB * E_1 * E_3 + sAB * E_2 * E_3;

  /*  dE_2 / dt  */
  ydot[1] = -(sA + dA + xA) * E_2 
            + xA + dA * E_1 + sA * E_2 * E_2;

  /*  dE_3 / dt  */
  ydot[2] = -(sB + dB + xB) * E_3 
            + xB + dB * E_1 + sB * E_3 * E_3;
}

/* deSolve */
void derivs_geosse_aux(int *neq, double *t, double *y, double *ydot, 
                       double *yout, int *ip) {
  do_derivs_geosse_aux(parms_geosse, y, ydot);
}

/* CVODES */
int derivs_geosse_aux_cvode(realtype t, N_Vector y, N_Vector ydot,
                           void *user_data) {
  do_derivs_geosse_aux(((UserData*) user_data)->p,
                      NV_DATA_S(y), NV_DATA_S(ydot));
  return 0;
}
