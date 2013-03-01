// -*-c++-*-

#ifndef __TIME_MACHINE_H
#define __TIME_MACHINE_H

#include <string>
#include <vector>

#include "Spline.h"

// The basic idea is pretty easy for non-Q based models, and that is
// what I'll work on first.  For models involving Q matrices, we need
// to skip a bunch of output variables on the first pass, and then
// come back and normalise the diagonals.


class TimeMachineFunction {
public:
  TimeMachineFunction(std::string name_, std::string func_,
		      bool nonnegative_, bool truncate_);
  void set(std::vector<double>::iterator p);
  double get(double t);
  bool is_constant;
  int np;

private:
  double check_ok(double x);

  typedef double tm_fun(double, 
			std::vector<double>::const_iterator,
			Spline *);
  bool nonnegative;
  bool truncate;
  std::string func_name;
  std::string variable_name;
  Spline *spline;
  tm_fun *f;
  std::vector<double> p_in;
};

class TimeMachine {
  // Default constructor.  
  //   We need to a vector of functions (identified by name)
  //   We need to know the number of target parameters
public: 
  TimeMachine(std::vector<std::string> names,
	      std::vector<std::string> funcs,
	      std::vector<bool> nonnegative,
	      std::vector<bool> truncate);

  void set(std::vector<double> pars);
  std::vector<double> get(double t);

private:
  // Number of high level (hyper) parameters (in).
  int np_in;
  // Number of underlying model parameters (out)
  int np_out;
  std::vector<double> p_out;
  // Number of functions.  Often equal to np_out, except for the case
  // of Q matrices.
  int nf;

  // Actual "functions":
  std::vector<TimeMachineFunction> functions;

  // Index of the starting parameter position in p_in (i.e., idx[i] is
  // the index within p_in for the first parameter of the ith
  // function).
  std::vector<int> idx;
};


double tm_fun_constant(double t, 
		       std::vector<double>::const_iterator pars_in,
		       Spline *spline);
double tm_fun_linear(double t, 
		     std::vector<double>::const_iterator pars_in,
		     Spline *spline);
double tm_fun_stepf(double t, 
		    std::vector<double>::const_iterator pars_in,
		    Spline *spline);
double tm_fun_sigmoid(double t,
		      std::vector<double>::const_iterator pars_in,
		      Spline *spline);
double tm_fun_spline(double t,
		     std::vector<double>::const_iterator pars_in,
		     Spline *spline);

#endif
