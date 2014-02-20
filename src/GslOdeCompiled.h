// -*-c++-*-

#ifndef __GSLODECOMPILED_H_
#define __GSLODECOMPILED_H_

#include <vector>
#include "GslOdeBase.h"

class GslOdeCompiled : public GslOdeBase {
public:
  GslOdeCompiled(SEXP extPtr, size_t size);

  // Our obligations:
  void derivs(double t, const double y[], double dydt[]);
  void set_pars(SEXP pars_);
  void clear_pars();

private:
  // Our internal machinery.

  // This is the typedef of the core derivative function.  This should
  // capture all of the ways that I currently use these functions.
  //                      neq, t,      pars,     y,              dydt
  typedef void deriv_func(int, double, double *, const double *, double *);
  // And here is the function itself.
  deriv_func *derivs_f;
  // Assumption is real parameters
  double *pars;
};

#endif
