#include <R.h>
#include <Rinternals.h>

void logistic(int neqs, double t, double *pars, const double *y, 
	      double *dydt) {
  const double r = pars[0], K = pars[1];
  int i;
  for ( i = 0; i < neqs; ++i )
    dydt[i] = r * y[i] * (1 - y[i]/K);
}
