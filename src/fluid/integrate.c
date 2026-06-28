#include <string.h>
#include <stdbool.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#define FLUID_INTERNAL
#include "internal.h"

static int copy_fields(
    const domain_t * domain,
    fluid_t * fluid,
    const bool is_deep_copy
){
  // two fields exist in fluid_t:
  //   1. s_x1_array:     store main field
  //   2. s_x1_array_int: store intermediate field
  // this function performs their memory copies,
  //   which are needed in the RK iterations
  // in particular, copy from the main to the intermediate is
  //   to be done at the beginning to set the initial RK field (for simplicity),
  //   while the other way should be done at the end
  //   to extract the RK result to the main field
  // since the n-step (main) field is used throughout the RK iteration,
  //   a deep copy is needed for the former, whereas a shallow copy is enough
  //   for the latter since the intermediate field is not used
  //   after the whole RK iterations
  for(size_t n = 0; n < NDIMS + 1; n++){
    field_t * field = fluid->fields[n];
    if(is_deep_copy){
      // use memcpy
      const size_t * mysizes = domain->s_x1_mysizes;
      const fftw_complex * restrict buf0 = field->s_x1_array;
      fftw_complex * restrict buf1 = field->s_x1_array_int;
      const size_t nitems = mysizes[0] * mysizes[1];
      memcpy(buf1, buf0, nitems * sizeof(fftw_complex));
    }else{
      // swap buffers
      fftw_complex * restrict tmp = field->s_x1_array;
      field->s_x1_array = field->s_x1_array_int;
      field->s_x1_array_int = tmp;
    }
  }
  return 0;
}

int fluid_integrate(
    const domain_t * domain,
    fluid_t * fluid,
    double * restrict dt
){
  // set 0-step RK values
  if(0 != copy_fields(domain, fluid, true)){
    return 1;
  }
  // RK iteration
  for(size_t rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
    // compute fields in physical space
    //   1. to decide time step size
    //   2. to compute convolution sum
    if(0 != compute_physical_fields(domain, fluid)){
      return 1;
    }
    // at the begining of RK,
    //   decide time step size using the physical velocity
    // NOTE: this exists inside the RK loop
    //   since decide_dt uses physical velocity
    if(0 == rkstep && 0 != decide_dt(domain, fluid, dt)){
      return 1;
    }
    // compute right-hand side of RK scheme: slopes
    //   i.e. "f" of dy/dt = f
    // NOTE: the only contribution is the advection
    if(0 != compute_slopes(domain, rkstep, fluid)){
      return 1;
    }
    // update fields: u^1, u^2, ..., u^{RKSTEPMAX-1}
    if(0 != update_fields(domain, rkstep, *dt, fluid)){
      return 1;
    }
  }
  // extract result
  if(0 != copy_fields(domain, fluid, false)){
    return 1;
  }
  return 0;
}

