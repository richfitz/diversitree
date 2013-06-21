#include "GslOdeCompiled.h"
#include <R.h>
#include "util.h"

GslOdeCompiled::GslOdeCompiled(SEXP extPtr, size_t size) :
  GslOdeBase(size) {
  derivs_f = get_deriv_func_from_ptr(extPtr);
}

void GslOdeCompiled::set_pars(SEXP pars_) {
  // Implicit short lifespan, as we don't copy (probably we
  // should...), but these are passed as read only elsewhere.
  // Declaring parameters as 
  //   double const * pars
  // should give us a mutable pointer to a constant vector (so that
  // pars can be reassigned, but the contents cannot be).  This
  // *should* work well here.
  pars = REAL(pars_);
}

void GslOdeCompiled::clear_pars() {
  pars = NULL;
}

void GslOdeCompiled::derivs(double t, const double y[], 
			    double dydt[]) {
  derivs_f(static_cast<int>(size()), t, pars, y, dydt);
}

