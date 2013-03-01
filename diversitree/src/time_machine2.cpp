#include <R.h>

#include "time_machine2.h"

TimeMachine::TimeMachine(std::vector<std::string> names,
			 std::vector<std::string> funcs,
			 std::vector<bool> nonnegative,
			 std::vector<bool> truncate) {
  nf = funcs.size();
  idx.resize(nf);
  np_in = 0;

  // For now -- will change when dealing with Q
  np_out = nf;
  p_out.resize(nf);

  for ( int i = 0; i < nf; i++ ) {
    idx[i] = np_in;
    TimeMachineFunction f = TimeMachineFunction(names[i], funcs[i], 
						nonnegative[i],
						truncate[i]);
    functions.push_back(f); // causes splines to be copied...
    np_in += f.np;
  }
}

void TimeMachine::set(std::vector<double> pars) {
  if ( (int)pars.size() != np_in )
    error("Expected %d parameters, recieved %d", np_in, pars.size());

  std::vector<TimeMachineFunction>::iterator f = functions.begin();
  std::vector<double>::iterator p_iter = pars.begin();
  for ( int i = 0; i < nf; i++ ) {
    f->set(p_iter + idx[i]);
    if ( f->is_constant )
      p_out[i] = f->get(0); // any time would work here.
    f++;
  }
}

// Probably also offer an iterator filling version for use within
// branches, eventually.
// 
// Make the iterator version of this an operator?  Not sure.
std::vector<double> TimeMachine::get(double t) {
  std::vector<TimeMachineFunction>::iterator f = functions.begin();
  for ( int i = 0; i < nf; i++ ) {
    if ( !f->is_constant )
      p_out[i] = f->get(t);
    f++;
  }

  return p_out;
}

// void normalise_q(dt_time_machine *obj, int is_const) {
//   for ( int i = 0; i < k; i++ ) {
//     if ( obj->const_q[i] == is_const ) {
//       int offset = idx_q + i;
//       double tmp = 0;
//       for ( int j = 0; j < k; j++ )
// 	if ( j != i )
// 	  tmp += p_out[offset + j*k];
//       p_out[offset + i*k] = -tmp;
//     }
//   }
// }


void TimeMachine::normalise_q(bool is_const) {
  for ( int i = 0; i < k; i++ ) {
    if ( const_q[i] == is_const ) {
      qi_out = q_out + i;
      double tmp = 0.0;
      for ( int j = 0; j < k; j++ )
	if ( j != i )
	  tmp += qi_out[j*k];
      qi_out[i*k] = -tmp;
    }
  }
}

//////////////////////////////////////////////////////////////////////

TimeMachineFunction::TimeMachineFunction(std::string name_,
					 std::string func_,
					 bool nonnegative_,
					 bool truncate_):
  variable_name(name_), func_name(func_),
  nonnegative(nonnegative_), truncate(truncate_) {

  // Splines disabled for now
  spline = NULL;

  is_constant = func_name == "constant.t";

  if ( func_name == "constant.t" ) {
    f = &tm_fun_constant;
    np = 1;
  } else if ( func_name == "linear.t" ) {
    f = &tm_fun_linear;
    np = 2;
  } else if ( func_name == "stepf.t" ) {
    f = &tm_fun_stepf;
    np = 3;
  } else if ( func_name == "sigmoid.t" ) {
    f = &tm_fun_sigmoid;
    np = 4;
  } else if ( func_name == "spline.t" ) {
    f = &tm_fun_spline;
    np = 2;
  } else {
    error("Unknown function type %s", func_name.c_str());
  }

  p_in.resize(np);
}

void TimeMachineFunction::set(std::vector<double>::iterator p) {
  std::copy(p, p+np, p_in.begin());
}

double TimeMachineFunction::get(double t) {
  return check_ok(f(t, p_in.begin(), spline));
}

double TimeMachineFunction::check_ok(double x) {
  const bool neg = x < 0;
  if ( neg ) {
    if ( truncate )
      x = 0;
    else if ( nonnegative )
      error("Value of %s (%s) must be nonnegative", 
	    variable_name.c_str(), func_name.c_str());
  }
  return x;
}

//////////////////////////////////////////////////////////////////////

double tm_fun_constant(double t,
		       std::vector<double>::const_iterator pars_in,
		       Spline *spline) {
  return *pars_in++;
}

double tm_fun_linear(double t,
		     std::vector<double>::const_iterator pars_in,
		     Spline *spline) {
  const double c = *pars_in++, m = *pars_in++;
  return c + m * t;
}

double tm_fun_stepf(double t,
		    std::vector<double>::const_iterator pars_in,
		    Spline *spline) {
  const double t0 = *pars_in++, t1 = *pars_in++, tc = *pars_in++;
  return t <= tc ? t0 : t1;
}

double tm_fun_sigmoid(double t,
		      std::vector<double>::const_iterator pars_in,
		      Spline *spline) {
  const double y0 = *pars_in++, y1 = *pars_in++, tmid = *pars_in++, 
    r = *pars_in++;
  return y0 + (y1 - y0)/(1 + exp(r * (tmid - t)));
}

double tm_fun_spline(double t,
		     std::vector<double>::const_iterator pars_in,
		     Spline *spline) {
  const double y0 = *pars_in++, y1 = *pars_in++;
  return y0 + (y1 - y0) * spline->eval(t);
}
