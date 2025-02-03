#include <stdbool.h>
#include <complex.h>
#include <fftw3.h>
#include "memory.h"
#include "domain.h"
#include "fluid.h"
#include "transform.h"
#define FLUID_INTERNAL
#include "internal.h"

// internal buffers
typedef struct {
  bool initialised;
  fftw_complex * masked;
} st_t;

static st_t st = {
  .initialised = false,
};

int compute_physical_fields(
    const domain_t * domain,
    fluid_t * fluid
){
  // compute array size,
  //   which is same in this function
  const size_t * mysizes = domain->s_x1_mysizes;
  const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
  // compute velocities in the physical space,
  //   which is done by transforming spectral velocity (iDFT)
  if(!st.initialised){
    st.masked = memory_fftw_calloc(nitems, sizeof(fftw_complex));
    st.initialised = true;
  }
  // for each field (momentum + scalar)
  for(size_t n = 0; n < NDIMS + 1; n++){
    field_t * field = fluid->fields[n];
    // input spectral field
    const fftw_complex * restrict iarray = field->s_x1_array_int;
    // output physical field
    double * restrict oarray = field->p_z1_array;
    const bool * restrict mask = fluid->s_x1_mask;
    fftw_complex * restrict masked = st.masked;
    // mask input array
    for(size_t index = 0; index < nitems; index++){
      masked[index] = mask[index] ? iarray[index] : 0.;
    }
    // iDFT, from spectral to physical
    if(0 != transform_s2p(domain, masked, oarray)){
      return 1;
    }
  }
  return 0;
}

