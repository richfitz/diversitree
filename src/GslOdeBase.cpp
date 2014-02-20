#include <Rcpp.h>

#include "GslOdeBase.h"

GslOdeBase::GslOdeBase(size_t size) : neq(size) {
  if ( size == 0 )
    Rf_error("Cannot create zero-sized problem");
  is_initialised = false;

  s = NULL;
  c = NULL;
  e = NULL;

  // Guess at sensible defaults, however, these will depend on the
  // scale.  However, something has to be set here, and these aren't
  // that different to what is used by deSolve's lsoda.
  atol = 0.0;
  rtol = 1e-6;
  hini = 1e-6;
  hmax = 1.0;
  hmin = 1e-10; // Arbitrary for now...

  // Runge-Kutta Cash-Karp as default step type (?)
  step_fun = gsl_odeiv2_step_rkck;

  // Fill GSL object with bits and bobs.
  sys.function  = helper_gsl_ode;
  sys.dimension = size;
  sys.params    = this;
}

GslOdeBase::~GslOdeBase() {
  free_gsl();
}

// [R] Given a set of parameters, time, and parameters, compute
// derivatives.  For testing, not for speed.  Note that we go through
// the virtual function 'derivs' here, defined by subclasses.  The
// set_pars() / clear_pars() functions are also defined by the
// subclasses, and are constrained only in that they take R objects
// (SEXPs).
std::vector<double>
GslOdeBase::r_derivs(double t, std::vector<double> y, SEXP pars_) {
  // TODO: The call to derivs() should be wrapped in a try (or
  // equivalent) so that pars are reset if the R call fails.
  if (y.size() != size())
    Rf_error("Incorrect input length (expected %d, got %d)",
	     size(), y.size());
  set_pars(pars_);

  std::vector<double> ret(size());
  derivs(t, &y[0], &ret[0]);
  clear_pars();

  return ret;
}

// [R] Run the system given initial conditions 'y', returning values
// at all times in 'times'.  The first value in times is taken as the
// starting time, and the initial conditions are returned for this
// time.  The matrix is size() rows and times.size() columns (each row
// representing a variable, each column representing a time).
Rcpp::NumericMatrix GslOdeBase::r_run(std::vector<double> times, 
				      std::vector<double> y_,
				      SEXP pars_) {
  set_pars(pars_);

  Rcpp::NumericMatrix ret(static_cast<int>(size()), 
			  static_cast<int>(times.size()-1));
  Rcpp::NumericMatrix::iterator out = ret.begin();

  std::vector<double>::iterator t = times.begin();
  set_state(*t++, y_); // This makes 'y' contains mutable state.

  while ( t != times.end() ) {
    advance(*t++);
    std::copy(y.begin(), y.end(), out);
    out += size();
  }

  clear_pars();

  return ret;
}

void GslOdeBase::r_set_control(Rcpp::List control) {
  std::vector<std::string> names = control.names();
  for (size_t i = 0; i < static_cast<size_t>(control.size()); i++)
    set_control1(names[i], control[static_cast<int>(i)]);
}

///////////////////////////////////////////////////////////////////////////////
// Private below here...

// Get ready to integrate a problem from time 't0' with initial
// conditions 'y_'
void GslOdeBase::set_state(double t0, const std::vector<double> y_) {
  init();
  must_be_initialised();

  if ( y_.size() != size() )
    Rf_error("Expected 'y' of size %d (recieved %d)", 
	     size(), y_.size());
  t = t0;
  y = y_;

  step_size = hini;

  reset_gsl();
}

// Possibly reallocate gsl space and set up the problem.  
void GslOdeBase::init() {
  if ( !is_initialised ) {
    alloc_gsl();
    is_initialised = true;
  }
}

// Simple check that everything that should be set up before anything
// can touch memory that has been allocated.
void GslOdeBase::must_be_initialised() const {
  if ( !is_initialised )
    Rf_error("Integrator not initialised");
}

// Control setting:
void GslOdeBase::set_control1(std::string key, SEXP value) {
  if ( key == "atol" )
    atol = Rcpp::as<double>(value);
  else if ( key == "rtol" )
    rtol = Rcpp::as<double>(value);
  else if ( key == "hini" )
    hini = Rcpp::as<double>(value);
  else if ( key == "hmax" )
    hmax = Rcpp::as<double>(value);
  else if ( key == "algorithm" )
    set_stepper_type(Rcpp::as<std::string>(value));
  else
    Rf_error("Unknown key `%s'", key.c_str());

  // Everything but hini will need gsl structures reallocated.
  if ( key != "hini" )
    free_gsl();
}

// Control setting: more work to set integrator type:
void GslOdeBase::set_stepper_type(const std::string type) {
  if ( type == "rk2" )
    step_fun = gsl_odeiv2_step_rk2;
  else if ( type == "rk4" )
    step_fun = gsl_odeiv2_step_rk4;
  else if ( type == "rkf45" )
    step_fun = gsl_odeiv2_step_rkf45;
  else if ( type == "rkck" )
    step_fun = gsl_odeiv2_step_rkck;
  else if ( type == "rk8pd" )
    step_fun = gsl_odeiv2_step_rk8pd;
  else
    Rf_error("Invalid stepper type specified");
}

// Stepping
// Advances problem 
void GslOdeBase::advance(double tmax) {
  must_be_initialised();

  // TODO: should check that tmax is finite and greater than t.
  while ( t < tmax ) {
    step(tmax);
  }
}

void GslOdeBase::step(double tmax) {
  int status = gsl_odeiv2_evolve_apply(e, c, s, 
				       &sys, &t, tmax,
				       &step_size, &y[0]);
  if(status != GSL_SUCCESS)
    Rf_error("GSL failure, return value = %d", status);
  if ( step_size > hmax ) // assuming positive sign step size here.
    step_size = hmax;
  else if ( step_size < hmin )
    Rf_error("Step size too small");
}

// Gsl control:
void GslOdeBase::alloc_gsl() {
  s = gsl_odeiv2_step_alloc(step_fun, size());
  if ( s == NULL )
    Rf_error("failed to allocate step object");

  c = gsl_odeiv2_control_y_new(atol, rtol);
  if ( c == NULL ) {
    gsl_odeiv2_step_free(s);
    Rf_error("failed to allocate control object");
  }

  e = gsl_odeiv2_evolve_alloc(size());
  if ( e == NULL ) {
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    Rf_error("failed to allocate evolve object");
  }
}

void GslOdeBase::free_gsl() {
  if ( is_initialised ) { // TODO: Not sure about logic?
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    s = NULL;
    c = NULL;
    e = NULL;
    is_initialised = false;
  }
}

// TODO: Not clear what should happen if either error is trigger.  Set
// flag broken and have must_be_initialised() fail there?
void GslOdeBase::reset_gsl() {
  must_be_initialised();

  const int flag_s = gsl_odeiv2_step_reset(s);
  if ( flag_s != GSL_SUCCESS )
    Rf_error("Error %d in resetting stepper", flag_s);

  const int flag_e = gsl_odeiv2_evolve_reset(e);
  if ( flag_e != GSL_SUCCESS )
    Rf_error("Error %d in resetting stepper", flag_e);
}

int helper_gsl_ode(double t, const double y[], double dydt[],
		   void *params) {
  GslOdeBase *pr = static_cast<GslOdeBase*>(params);
  pr->derivs(t, y, dydt);
  return GSL_SUCCESS;
}
