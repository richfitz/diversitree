// -*-c++-*-

// TODO: There is no facilities for checking parameter length yet.
// This is hard to do at the moment because different models might
// have different requirements (e.g., what is the parameter
// requirement for the R versions?  Possibly supply a function?)

// TODO: Don't much like the variation in argument lists:
//   derivs -- (t, y, dydt)
//   R function -- (x, t, pars)
//   C function -- (t, y, dydt)
//   deSolve -- (t, y, pars)
//   exposed derivs -- (y, t, pars)
//   
//   Settle on t, y, [pars], [dydt]

// There are a couple of bits that are strongly tied to R: the
// handling of infinities and the Matrix return (from the less useful
// vector<vector<double>>)

#ifndef __GslOdeBase_H_
#define __GslOdeBase_H_

#include <vector>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_nan.h>

#include <Rcpp.h>

class GslOdeBase {
  // Need to take problem (R function) and size (integer) as
  // arguments.  It might be nice to delay argument size until later
  // when the initial conditions are provided.  But that relies on
  // very late allocation of data structures.
  // 
  // It might be nicer to have this inherit from the base class, but
  // then it's all templated and that's going to be confusing to pass
  // along.
  
public:
  // TODO: Should have (or prevent) copy constructor and assignment
  // (Rule of Three).
  GslOdeBase(size_t size);
  // Destructor should be declared virtual -- see:
  // http://stackoverflow.com/questions/461203/when-to-use-virtual-destructors
  // http://www.gotw.ca/publications/mill18.htm
  virtual ~GslOdeBase();

  // Dimension of the problem.
  // TODO: Because readonly, reference size directly, not via getter?
  size_t size() const { return neq; }
  // This is because CRAN won't allow Rcpp to wrap size_t on Windows
  // http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2012-July/004016.html
  int r_size() const { return static_cast<int>(size()); }

  // For GSL; expected name and arguments for computing derivatives.
  virtual void derivs(double t, const double y[], double dydt[]) = 0;

  // [R] Given a set of parameters, time, and parameters, compute
  // derivatives.  For testing, not for speed.
  std::vector<double> r_derivs(double t, std::vector<double> y, SEXP pars_);

  // [R] Run the system given initial conditions 'y', returning values
  // at all times in 'times'.  The first value in times is taken as
  // the starting time, and the initial conditions are returned for
  // this time.
  Rcpp::NumericMatrix r_run(std::vector<double> times, 
			    std::vector<double> y, 
			    SEXP pars_);

  // [R] Set control parameters (tolerances, integrator type, etc).
  void r_set_control(Rcpp::List control);

protected:
  virtual void set_pars(SEXP pars) = 0;
  virtual void clear_pars() = 0;
  
private:
  // Construction and control
  void set_state(double t0, const std::vector<double> y_);
  // Ensure all gsl memory allocated.
  void init();
  // Assert (or throw R error) that gsl ready to go
  void must_be_initialised() const;

  // Control setting
  void set_control1(std::string key, SEXP value);
  void set_stepper_type(std::string type);

  // Stepping
  // Advance time from t to tmax:
  void advance(double tmax);
  // Advance one step, up to maximum of tmax:
  void step(double tmax);

  // GSL internals
  void alloc_gsl();
  void free_gsl();
  void reset_gsl();

  // Problem dimension (number of equations)
  size_t neq;

  // Variables (length neq)
  std::vector<double> y;

  // Keep track of if gsl structures are allocated.
  bool is_initialised;

  // Integrator control:
  double hini, hmax, hmin, atol, rtol;

  // Suggested size of the next step, same as the actual size of the
  // last step once a step has been taken.
  double step_size;

  // Time, reset after init(), incremented by stepping or running.
  double t;

  // GSL internals
  gsl_odeiv2_system sys;
  gsl_odeiv2_step    *s;
  gsl_odeiv2_control *c;
  gsl_odeiv2_evolve  *e;
  const gsl_odeiv2_step_type *step_fun;
};

int helper_gsl_ode(double t, const double y[], double dydt[],
		   void *params);

#endif
