// -*-c++-*-

#ifndef __GSLODETIME_H_
#define __GSLODETIME_H_

#include <vector>

#include "TimeMachine.h"
#include "GslOdeBase.h"

// There is some overlap here with GslOdeCompiled, which is
// undesirable.
class GslOdeTime : public GslOdeBase {
public:
  GslOdeTime(SEXP extPtr, size_t size, TimeMachine tm_);
  
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
  // The time machine!
  TimeMachine tm;
};

#endif
