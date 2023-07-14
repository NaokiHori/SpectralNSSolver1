#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#define FLUID_INTERNAL
#include "internal.h"

static inline double compute_factor(
    const double diffusivity,
    const double k2,
    const double weight,
    const double dt
){
  return exp(diffusivity * k2 * weight * dt);
}

static int update_field(
    const domain_t * domain,
    const size_t rkstep,
    const double * restrict coef_as,
    const double * restrict coef_cs,
    const double dt,
    const double diffusivity,
    // n-step field
    const fftw_complex * restrict array0,
    fftw_complex * const restrict slopes[RKSTEPMAX],
    fftw_complex * restrict array1
){
  const size_t * mysizes = domain->s_x1_mysizes;
  const double * restrict xfreqs = domain->x1_xfreqs;
  const double * restrict yfreqs = domain->x1_yfreqs;
#if NDIMS == 3
  const double * restrict zfreqs = domain->x1_zfreqs;
#endif
  // u^n contribution
  {
#if NDIMS == 2
    const size_t nitems = mysizes[0] * mysizes[1];
#else
    const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
#endif
    for(size_t index = 0; index < nitems; index++){
      array1[index] = array0[index];
    }
  }
  // append f^k contributions
  for(size_t l = 0; l < rkstep + 1; l++){
    const double coef_a = coef_as[l];
    const double coef_c = coef_cs[l];
    const fftw_complex * restrict slope = slopes[l];
    if(0. == coef_a){
      continue;
    }
#if NDIMS == 2
    for(size_t index = 0, j = 0; j < mysizes[1]; j++){
      const double ky = yfreqs[j];
      for(size_t i = 0; i < mysizes[0]; i++, index++){
        const double kx = xfreqs[i];
        const double k2 =
          + 1. * kx * kx
          + 1. * ky * ky;
        const double e = compute_factor(diffusivity, k2, coef_c, dt);
        array1[index] += coef_a * dt * e * slope[index];
      }
    }
#else
    for(size_t index = 0, k = 0; k < mysizes[2]; k++){
      const double kz = zfreqs[k];
      for(size_t j = 0; j < mysizes[1]; j++){
        const double ky = yfreqs[j];
        for(size_t i = 0; i < mysizes[0]; i++, index++){
          const double kx = xfreqs[i];
          const double k2 =
            + 1. * kx * kx
            + 1. * ky * ky
            + 1. * kz * kz;
          const double e = compute_factor(diffusivity, k2, coef_c, dt);
          array1[index] += coef_a * dt * e * slope[index];
        }
      }
    }
#endif
  }
  // compute new field
  {
    const double coef_c = coef_cs[rkstep + 1];
#if NDIMS == 2
    for(size_t index = 0, j = 0; j < mysizes[1]; j++){
      const double ky = yfreqs[j];
      for(size_t i = 0; i < mysizes[0]; i++, index++){
        const double kx = xfreqs[i];
        const double k2 =
          + 1. * kx * kx
          + 1. * ky * ky;
        const double e = compute_factor(diffusivity, k2, coef_c, dt);
        array1[index] /= e;
      }
    }
#else
    for(size_t index = 0, k = 0; k < mysizes[2]; k++){
      const double kz = zfreqs[k];
      for(size_t j = 0; j < mysizes[1]; j++){
        const double ky = yfreqs[j];
        for(size_t i = 0; i < mysizes[0]; i++, index++){
          const double kx = xfreqs[i];
          const double k2 =
            + 1. * kx * kx
            + 1. * ky * ky
            + 1. * kz * kz;
          const double e = compute_factor(diffusivity, k2, coef_c, dt);
          array1[index] /= e;
        }
      }
    }
#endif
  }
  return 0;
}

int update_fields(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  // update fields using computed slopes
  for(size_t n = 0; n < NDIMS + 1; n++){
    field_t * field = fluid->fields[n];
    const double diffusivity = field->diffusivity;
    if(0 != update_field(
        domain,
        rkstep,
        runge_kutta_coef_as[rkstep],
        runge_kutta_coef_cs,
        dt,
        diffusivity,
        // n-step field
        field->s_x1_array,
        // slopes
        field->s_x1_slopes,
        // intermediate field
        field->s_x1_array_int
    )) return 1;
  }
  return 0;
}

