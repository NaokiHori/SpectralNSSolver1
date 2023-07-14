#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

#define RKSTEPMAX 4

extern const double runge_kutta_cfl;

// Butcher tableau
extern const double runge_kutta_coef_as[RKSTEPMAX][RKSTEPMAX];
extern const double runge_kutta_coef_cs[RKSTEPMAX + 1];

#endif // RUNGE_KUTTA_H
