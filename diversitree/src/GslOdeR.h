// -*-c++-*-
#ifndef __GSLODER_H_
#define __GSLODER_H_

#include "GslOdeBase.h"

class GslOdeR : public GslOdeBase {
public:
  GslOdeR(SEXP fun_, SEXP env_, size_t size);

  // Our obligations:
  void derivs(double t, const double y[], double dydt[]);
  void set_pars(SEXP pars_);
  void clear_pars();

private:
  // Our internal machinery.
  // TODO: fun should also take t as argument.
  SEXP fun;  // R derivatives function (signature: function(t, y, pars))
  SEXP env;  // Environment to evaluate fun within
  SEXP pars; // Transient parameters object, passed to fun

  // Later (for derivs, so could just be done there?)
  SEXP target(double t, SEXP y);
};

#endif
