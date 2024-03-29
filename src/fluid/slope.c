#include <stdbool.h>
#include <string.h>
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
  // arrays used inside the function "convolute"
  // store product of two arrays in the physical domain
  double * restrict p_y1_buf;
  // store product in the spectral domain,
  //   i.e. iDFT(pbuf)
  fftw_complex * restrict s_x1_buf;
} st_t;
static st_t st = {
  .initialised = false,
};

static int convolute(
    const domain_t * domain,
    const double * restrict parr0,
    const double * restrict parr1,
    fftw_complex * restrict sbuf
){
  // compute convolution sum of two arrays
  // NOTE: two arrays should already be in the physical space (i.e. after iDFT-ed)
  const size_t * mysizes = domain->p_y1_mysizes;
  const size_t nitems = mysizes[0] * mysizes[1];
  double * restrict pbuf = st.p_y1_buf;
  // compute product in the physical domain
  for(size_t index = 0; index < nitems; index++){
    pbuf[index] = parr0[index] * parr1[index];
  }
  // go back to the spectral domain
  if(0 != transform_p2s(domain, pbuf, sbuf)){
    return 1;
  }
  return 0;
}

static int compute_adv(
    const domain_t * domain,
    const fluid_t * fluid,
    const double * restrict q,
    fftw_complex * restrict slope
){
  // evaluate advective terms of the given field "q",
  //   i.e. - d(u_j q) / dx_j
  const size_t * mysizes = domain->s_x1_mysizes;
  const double * restrict xfreqs = domain->x1_xfreqs;
  const double * restrict yfreqs = domain->x1_yfreqs;
  const double * restrict ux = fluid->fields[enum_ux]->p_y1_array;
  const double * restrict uy = fluid->fields[enum_uy]->p_y1_array;
  fftw_complex * restrict buf = st.s_x1_buf;
  // zero-clear buffer
  // although not necessary, I do this
  //   just to treat all directions consistently below
  memset(slope, 0, sizeof(fftw_complex) * mysizes[0] * mysizes[1]);
  // - d(ux q)/dx
  if(0 != convolute(domain, ux, q, buf)){
    return 1;
  }
  for(size_t index = 0, j = 0; j < mysizes[1]; j++){
    for(size_t i = 0; i < mysizes[0]; i++, index++){
      const double kx = xfreqs[i];
      slope[index] -= I * kx * buf[index];
    }
  }
  // - d(uy q)/dy
  if(0 != convolute(domain, uy, q, buf)){
    return 1;
  }
  for(size_t index = 0, j = 0; j < mysizes[1]; j++){
    const double ky = yfreqs[j];
    for(size_t i = 0; i < mysizes[0]; i++, index++){
      slope[index] -= I * ky * buf[index];
    }
  }
  return 0;
}

static int project_velocity(
    const domain_t * domain,
    const size_t rkstep,
    fluid_t * fluid
){
  const size_t * mysizes = domain->s_x1_mysizes;
  const double * restrict xfreqs = domain->x1_xfreqs;
  const double * restrict yfreqs = domain->x1_yfreqs;
  fftw_complex * restrict slopeux = fluid->fields[enum_ux]->s_x1_slopes[rkstep];
  fftw_complex * restrict slopeuy = fluid->fields[enum_uy]->s_x1_slopes[rkstep];
  for(size_t index = 0, j = 0; j < mysizes[1]; j++){
    const double ky = yfreqs[j];
    for(size_t i = 0; i < mysizes[0]; i++, index++){
      const double kx = xfreqs[i];
      if(0 == i && 0 == j){
        continue;
      }
      const double k2 =
        + 1. * kx * kx
        + 1. * ky * ky;
      const fftw_complex ip =
        + 1. * kx * slopeux[index]
        + 1. * ky * slopeuy[index];
      slopeux[index] -= kx / k2 * ip;
      slopeuy[index] -= ky / k2 * ip;
    }
  }
  return 0;
}

int compute_slopes(
    const domain_t * domain,
    const size_t rkstep,
    fluid_t * fluid
){
  if(!st.initialised){
    // allocate internal buffers
    const size_t s_x1_nitems = domain->s_x1_mysizes[0] * domain->s_x1_mysizes[1];
    const size_t p_y1_nitems = domain->p_y1_mysizes[0] * domain->p_y1_mysizes[1];
    st.s_x1_buf = memory_fftw_calloc(s_x1_nitems, sizeof(fftw_complex));
    st.p_y1_buf = memory_fftw_calloc(p_y1_nitems, sizeof(double));
    st.initialised = true;
  }
  // repeat the same thing for each field 
  // NOTE: velocity in each direction and one scalar field
  for(size_t n = 0; n < NDIMS + 1; n++){
    const double * restrict iarray = fluid->fields[n]->p_y1_array;
    fftw_complex * restrict oarray = fluid->fields[n]->s_x1_slopes[rkstep];
    if(0 != compute_adv(domain, fluid, iarray, oarray)){
      return 1;
    }
  }
  // evaluate correction terms to make the velocity field non-solenoidal
  if(0 != project_velocity(domain, rkstep, fluid)){
    return 1;
  }
  return 0;
}

