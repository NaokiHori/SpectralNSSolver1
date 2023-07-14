#include "runge_kutta.h"

// for RK4, max is 2.8 according to eigenvalue analysis
// the smaller the more stable and more accurate
const double runge_kutta_cfl = 2.;

// Butcher-tableau, a_ij (b_j are merged)
const double runge_kutta_coef_as[RKSTEPMAX][RKSTEPMAX] = {
  {1. / 2.,      0.,      0.,      0.},
  {     0., 1. / 2.,      0.,      0.},
  {     0.,      0.,      1.,      0.},
  {1. / 6., 1. / 3., 1. / 3., 1. / 6.},
};

// Butcher-tableau, c_i
// since c_i - c_j are needed, n-step contribution is appended
const double runge_kutta_coef_cs[RKSTEPMAX + 1] = {
  0.0,
  0.5,
  0.5,
  1.0,
  1.0,
};

