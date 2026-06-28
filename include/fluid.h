#if !defined(FLUID_H)
#define FLUID_H

#include <stdbool.h>
#include <fftw3.h>
#include "domain.h"
#include "runge_kutta.h"

typedef struct {
  // array in spectral domain, x1 pencil 
  fftw_complex * restrict s_x1_array;
  // array in physical domain, z1 pencil 
  double * restrict p_z1_array;
  // storage to store intermediate field A^{0,1,2,...,RKSTEPMAX-1} of RK scheme
  fftw_complex * restrict s_x1_array_int;
  // storage to store slopes (right-hand-side terms) of RK schemes
  fftw_complex * restrict s_x1_slopes[RKSTEPMAX];
  // diffusivity of this quantity
  double diffusivity;
} field_t;

typedef enum {
  enum_ux,
  enum_uy,
  enum_uz,
  enum_sc,
} field_number_t;

typedef struct {
  // flow fields, momentum in each direction and scalar field
  field_t * fields[NDIMS + 1];
  // 2/3 dealiasing mask
  bool * s_x1_mask;
} fluid_t;

extern int fluid_init(
    const char dirname[],
    const domain_t * domain,
    fluid_t * fluid
);

extern int fluid_load(
    const char dirname[],
    const domain_t * domain,
    fluid_t * fluid
);

extern int fluid_save(
    const char dirname[],
    const domain_t * domain,
    const fluid_t * fluid
);

extern int fluid_integrate(
    const domain_t * domain,
    fluid_t * fluid,
    double * dt
);

#endif // FLUID_H
