#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#include "memory.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"

static inline size_t iabs(
    const int val
){
  return val < 0 ? - val : val;
}

static int allocate_and_init_mask(
    const domain_t * domain,
    bool ** mask
){
  // allocate (this is x1 pencil in spectral domain)
  const size_t * mysizes = domain->s_x1_mysizes;
  const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
  *mask = memory_fftw_calloc(nitems, sizeof(bool));
  // create mask
  const size_t * glsizes = domain->s_glsizes;
  const int * xwaves = domain->x1_xwaves;
  const int * ywaves = domain->x1_ywaves;
  const int * zwaves = domain->x1_zwaves;
  for(size_t index = 0, k = 0; k < mysizes[2]; k++){
    const int kz = zwaves[k];
    for(size_t j = 0; j < mysizes[1]; j++){
      const int ky = ywaves[j];
      for(size_t i = 0; i < mysizes[0]; i++, index++){
        const int kx = xwaves[i];
        if(
               iabs(kx) >= glsizes[0] / 3
            || iabs(ky) >= glsizes[1] / 3
            || iabs(kz) >= glsizes[2] / 3
        ){
          (*mask)[index] = false;
        }else{
          (*mask)[index] = true;
        }
      }
    }
  }
  return 0;
}

static int allocate_and_init_field(
    const domain_t * domain,
    field_t ** field,
    const double diffusivity
){
  // structure itself
  *field = memory_calloc(1, sizeof(field_t));
  // main and sub fields in spectral domain
  {
    const size_t * mysizes = domain->s_x1_mysizes;
    const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
    (*field)->s_x1_array     = memory_fftw_calloc(nitems, sizeof(fftw_complex));
    (*field)->s_x1_array_int = memory_fftw_calloc(nitems, sizeof(fftw_complex));
    for(size_t rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
      (*field)->s_x1_slopes[rkstep] = memory_fftw_calloc(nitems, sizeof(fftw_complex));
    }
  }
  // auxiliary field in physical domain to compute convolution sum
  {
    const size_t * mysizes = domain->p_z1_mysizes;
    const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
    (*field)->p_z1_array = memory_fftw_calloc(nitems, sizeof(double));
  }
  (*field)->diffusivity = diffusivity;
  return 0;
}

int fluid_init(
    const char dirname[],
    const domain_t * domain,
    fluid_t * fluid
){
  // load non-dimensional parameters
  double Re = 0.;
  double Sc = 0.;
  if(0 != config.get_double("Re", &Re)){
    return 1;
  }
  if(0 != config.get_double("Sc", &Sc)){
    return 1;
  }
  // allocate and prepare 2/3 dealiasing mask
  if(0 != allocate_and_init_mask(domain, &fluid->s_x1_mask)){
    return 1;
  }
  // allocate buffers for flow each field and set diffusivity
  if(0 != allocate_and_init_field(domain, &fluid->fields[enum_ux], 1. / Re     )){
    return 1;
  }
  if(0 != allocate_and_init_field(domain, &fluid->fields[enum_uy], 1. / Re     )){
    return 1;
  }
  if(0 != allocate_and_init_field(domain, &fluid->fields[enum_uz], 1. / Re     )){
    return 1;
  }
  if(0 != allocate_and_init_field(domain, &fluid->fields[enum_sc], 1. / Re / Sc)){
    return 1;
  }
  // load initial condition from files
  if(0 != fluid_load(dirname, domain, fluid)){
    return 1;
  }
  return 0;
}

