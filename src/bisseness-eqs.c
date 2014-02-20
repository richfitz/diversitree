/*
 * 
 *
 * bisse_ness-eqs.c
 *  
 *
 * Created by Karen Magnuson-Ford on February 16, 2012.
 * Copyright 2012. All rights reserved.
 *
 * Modified by RGF 4 March 2013.
 *
 */
#include <R.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bisseness(double *pars, const double *y, double *ydot) {
  double E0 = y[0], E1 = y[1], D0 = y[2], D1 = y[3];
  double la0 = pars[0], la1 = pars[1], mu0 = pars[2], mu1 = pars[3],
    q01 = pars[4], q10 = pars[5],
    p0c = pars[6], p0a = pars[7], p1c = pars[8], p1a = pars[9];

  ydot[0] = -(mu0 + q01 + la0) * E0 + 
    la0 * E0 * E0 * (1 - p0c) + mu0 + q01 * E1 + 
    la0 * E0 * E1 * p0c * p0a + la0 * E1 * E1 * p0c * (1 - p0a);
  ydot[1] = -(mu1 + q10 + la1) * E1 + 
    la1 * E1 * E1 * (1 - p1c) + mu1 + q10 * E0 + 
    la1 * E0 * E1 * p1c * p1a + la1 * E0 * E0 * p1c * (1 - p1a);
  ydot[2] = -(mu0 + q01 + la0) * D0 + 
    2 * la0 * E0 * D0 * (1 - p0c) + q01 * D1 + 
    (E0 * D1 + E1 * D0) * la0 * p0c * p0a + 
    2 * la0 * E1 * D1 * p0c * (1 - p0a);
  ydot[3] = -(mu1 + q10 + la1) * D1 + 
    2 * la1 * E1 * D1 * (1 - p1c) + q10 * D0 + 
    (E1 * D0 + E0 * D1) * la1 * p1c * p1a + 
    2 * la1 * E0 * D0 * p1c * (1 - p1a);
}

/* 
   Wrap these up for GslOde, which allows time and neq to be here too.
*/
void derivs_bisseness_gslode(int neqs, double t, double *pars, 
			 const double *y, double *dydt) {
  do_derivs_bisseness(pars, y, dydt);
}
