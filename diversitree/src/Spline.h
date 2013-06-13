// -*-c++-*-
#ifndef SPLINE_H_
#define SPLINE_H_

#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class Spline {
public:
  Spline();
  ~Spline();
  Spline(const Spline& x);
  void reset();
  void init_self();
  void init(std::vector<double> x_, std::vector<double> y_);
  double eval(double u) const;
  void add_point(double xi, double yi);
  size_t size() const;

protected:
  std::vector<double> x, y;

  // Specific to GSL splines.
  gsl_interp_accel *acc;
  gsl_spline *spline;  

  void gsl_free_spline();
  void gsl_free_acc();

  void do_gsl_alloc_spline(std::vector<double> x, std::vector<double> y);
};

#endif // SPLINE_H_
