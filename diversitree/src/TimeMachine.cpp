#include <R.h>

#include "TimeMachine.h"

TimeMachine::TimeMachine(std::vector<std::string> names,
			 std::vector<std::string> funcs,
			 std::vector<bool> nonnegative,
			 std::vector<bool> truncate,
			 size_t k,
			 std::vector<double> spline_t,
			 std::vector<double> spline_y) {
  nf = funcs.size();
  idx.resize(nf);
  target.resize(nf);
  np_in = 0;

  np_out = nf;
  p_out.resize(nf);

  if ( spline_t.size() > 0 )
    spline.init(spline_t, spline_y);

  for (size_t i = 0; i < nf; i++) {
    idx[i] = np_in;
    target[i] = i;
    TimeMachineFunction f = TimeMachineFunction(names[i], funcs[i], 
						nonnegative[i],
						truncate[i],
						&spline);
    functions.push_back(f);
    np_in += f.np;
  }

  setup_q(k);
}

void TimeMachine::set(std::vector<double> pars) {
  if (pars.size() != np_in)
    error("Expected %d parameters, recieved %d", np_in, pars.size());
  // Only go through the extra effort below if the parameters differ.
  if ( pars == p_in )
    return;

  p_in = pars;

  std::vector<TimeMachineFunction>::iterator f = functions.begin();
  std::vector<double>::iterator p_iter = pars.begin();
  for (size_t i = 0; i < nf; i++) {
    // TODO: Cast error here
    f->set(p_iter + static_cast<int>(idx[i]));
    if ( f->is_constant )
      p_out[target[i]] = f->get(0); // any time would work here.
    f++;
  }
  if ( k > 0 )
    normalise_q(true);
}

// Probably also offer an iterator filling version for use within
// branches, eventually.
// 
// Make the iterator version of this an operator?  Not sure.
std::vector<double> TimeMachine::get(double t) {
  std::vector<TimeMachineFunction>::iterator f = functions.begin();
  for (size_t i = 0; i < nf; i++) {
    if ( !f->is_constant )
      p_out[target[i]] = f->get(t);
    f++;
  }
  if ( k > 0 )
    normalise_q(false);

  return p_out;
}

// This is not actually done - we need this to work on a vector of
// time, returning a matrix...
std::vector<double> TimeMachine::getv(double t) {
  std::vector<double> p = get(t);
  if ( k > 0 ) {
    std::vector<double> ret;
    std::vector<double>::iterator pi = p.begin();

    for (size_t i = 0; i < idx_q_out; i++)
      ret.push_back(*pi++);
    for (size_t i = 0; i < k; i++)
      for (size_t j = 0; j < k; j++, pi++)
	if ( i != j )
	  ret.push_back(*pi);
    p = ret;
  }
  return p;
}

std::vector<std::string> TimeMachine::names() const {
  std::vector<std::string> ret;
  std::vector<TimeMachineFunction>::const_iterator f;
  for ( f = functions.begin(); f != functions.end(); f++ )
    ret.push_back(f->name());
  return ret;
}

void TimeMachine::setup_q(size_t k_) {
  k = k_;
  if ( k > 0 ) {
    // More space in the output vector
    np_out += k;
    p_out.resize(np_out);

    // Index of start of Q matrix in output parameters
    idx_q_out = np_out - k*k;
    idx_q_f   = nf     - k*(k-1);

    // Space to indicate when rows/columns(?) of Q are constant
    const_q.resize(k, true);
    std::vector<TimeMachineFunction>::iterator f = functions.begin();
    std::advance(f, idx_q_f);
    for (size_t i = 0; i < k; i++) {
      for (size_t j = 0; j < k - 1; j++) {
	const_q[i] = const_q[i] && f->is_constant;
	f++;
      }
    }
      
    // Work out what the actual place things belong (recomputing the
    // target mapping) given the input Q matrix is transposed relative
    // to true Q matrix, and lacks diagonals.  This is pretty opaque
    // :(
    for (size_t i = 0, l = idx_q_f; i < k; i++)
      for (size_t j = 0; j < k; j++)
	if (i != j)
	  target[l++] = i + j * k + idx_q_f;
  }
}

void TimeMachine::normalise_q(bool is_const) {
  for (size_t i = 0; i < k; i++) {
    if ( const_q[i] == is_const) {
      size_t offset = idx_q_out + i;
      double tmp = 0.0;
      for (size_t j = 0; j < k; j++)
	if ( j != i )
	  tmp += p_out[offset + j*k];
      p_out[offset + i*k] = -tmp;
    }
  }
}


//////////////////////////////////////////////////////////////////////

TimeMachineFunction::TimeMachineFunction(std::string name_,
					 std::string func_,
					 bool nonnegative_,
					 bool truncate_,
					 Spline *spline_):
  variable_name(name_), func_name(func_),
  nonnegative(nonnegative_), truncate(truncate_),
  spline(spline_) {
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
    if ( spline == NULL )
      error("Should not be able to get here!");
  } else {
    error("Unknown function type %s", func_name.c_str());
  }

  p_in.resize(np);
}

void TimeMachineFunction::set(std::vector<double>::iterator p) {
  // TODO: Cast error here
  std::copy(p, p+static_cast<int>(np), p_in.begin());
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
