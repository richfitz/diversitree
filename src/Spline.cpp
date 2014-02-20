#include "Spline.h"

Spline::Spline() {
  acc    = NULL;
  spline = NULL;
}

Spline::~Spline() {
  gsl_free_spline();
  gsl_free_acc();
}

Spline::Spline(const Spline& obj) : x(obj.x), y(obj.y) {
  acc = NULL;
  spline = NULL;
  if (size() > 0)
    init_self();
}

void Spline::gsl_free_spline() {
  if ( spline != NULL ) {
    gsl_spline_free(spline);
    spline = NULL;
  }
}

void Spline::gsl_free_acc() {
  if ( acc != NULL ) {
    gsl_interp_accel_free(acc);
    acc = NULL;
  }
}

void Spline::do_gsl_alloc_spline(std::vector<double> x_,
				 std::vector<double> y_) {
  const size_t n = size();
  gsl_free_spline();
  gsl_free_acc();
  acc    = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, &x_.front(), &y_.front(), n);
}

void Spline::init(std::vector<double> x_, std::vector<double> y_) {
  x      = x_;
  y      = y_;
  do_gsl_alloc_spline(x, y);
}

double Spline::eval(double u) const {
  return gsl_spline_eval(spline, u, acc);
}

void Spline::reset() {
  x.clear();
  y.clear();
  gsl_free_spline();
}

void Spline::add_point(double xi, double yi) {
  x.push_back(xi);
  y.push_back(yi);
}

void Spline::init_self() {
  const size_t n = size();
  gsl_free_spline();
  gsl_free_acc();
  spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, &x.front(), &y.front(), n);  
}

size_t Spline::size() const {
  return x.size();
}
