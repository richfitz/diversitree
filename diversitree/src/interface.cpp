#include <Rcpp.h>

#include "GslOdeBase.h"

#include "GslOdeR.h"
#include "GslOdeCompiled.h"
#include "GslOdeTime.h"

#include "TimeMachine.h"

// Need to do this.
RCPP_EXPOSED_CLASS(TimeMachine)

RCPP_MODULE(diversitree) {
  using namespace Rcpp;

  class_<GslOdeBase>("GslOdeBase")
    .property("size",      &GslOdeBase::r_size)
    .method("derivs",      &GslOdeBase::r_derivs)
    .method("run",         &GslOdeBase::r_run)
    .method("set_control", &GslOdeBase::r_set_control)
    ;

  class_<GslOdeR>("GslOdeR")
    .constructor<SEXP, SEXP, int>()
    .derives<GslOdeBase>("GslOdeBase")
    ;

  class_<GslOdeCompiled>("GslOdeCompiled")
    .constructor<SEXP, int>()
    .derives<GslOdeBase>("GslOdeBase")
    ;

  class_<GslOdeTime>("GslOdeTime")
    .constructor<SEXP, int, TimeMachine>()
    .derives<GslOdeBase>("GslOdeBase")
    ;

  // Or with typedefs
  // typedef vec_str  = std::vector<std::string>;
  // typedef vec_bool = std::vector<bool>;
  // .constructor<vec_str, vec_str, vec_bool, vec_bool>()
  class_<TimeMachine>("TimeMachine")
    .constructor< std::vector<std::string>,
		  std::vector<std::string>,
		  std::vector<bool>,
		  std::vector<bool>,
		  int,
		  std::vector<double>,
		  std::vector<double> >()
    .property("size",  &TimeMachine::r_size)
    .property("names", &TimeMachine::names)
    .method("set",     &TimeMachine::set)
    .method("get",     &TimeMachine::get)
    .method("getv",    &TimeMachine::getv)
    ;
}
