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
  // store product in physical domain
  double * restrict p_y1_buf;
  // store product in spectral domain
  fftw_complex * restrict s_x1_buf;
  // arrays to store output
  fftw_complex * restrict s_x1_advs[NDIMS + 1];
} st_t;
static st_t st = {
  .initialised = false,
};

static int convolute(
    const domain_t * domain,
    const double * restrict arr0,
    const double * restrict arr1
){
  // compute convolution sum of two arrays
  // NOTE: two arrays should already be in the physical space (i.e. after iDFT-ed)
  const size_t * mysizes = domain->p_y1_mysizes;
  const size_t nitems = mysizes[0] * mysizes[1];
  double * restrict pbuf = st.p_y1_buf;
  fftw_complex * restrict sbuf = st.s_x1_buf;
  // compute product in the physical domain
  for(size_t index = 0; index < nitems; index++){
    pbuf[index] = arr0[index] * arr1[index];
  }
  // go back to the spectral domain
  if(0 != transform_p2s(domain, pbuf, sbuf)) return 1;
  return 0;
}

static int kernel(
    const domain_t * domain,
    const fluid_t * fluid,
    const double * restrict q,
    fftw_complex * restrict adv
){
  // evaluate advective terms of the given field "q"
  // NOTE: do not restrict-qualify
  const size_t * mysizes = domain->s_x1_mysizes;
  const double * restrict xfreqs = domain->x1_xfreqs;
  const double * restrict yfreqs = domain->x1_yfreqs;
  const double * restrict ux = fluid->fields[enum_ux]->p_y1_array;
  const double * restrict uy = fluid->fields[enum_uy]->p_y1_array;
  fftw_complex * restrict buf = st.s_x1_buf;
  // zero-clear buffer
  // although not necessary, I do this
  //   just to treat all directions consistently below
  memset(adv, 0, sizeof(fftw_complex) * mysizes[0] * mysizes[1]);
  // d(ux q)/dx
  if(0 != convolute(domain, ux, q)) return 1;
  for(size_t index = 0, j = 0; j < mysizes[1]; j++){
    for(size_t i = 0; i < mysizes[0]; i++, index++){
      const double kx = xfreqs[i];
      adv[index] += I * kx * buf[index];
    }
  }
  // d(uy q)/dy
  if(0 != convolute(domain, uy, q)) return 1;
  for(size_t index = 0, j = 0; j < mysizes[1]; j++){
    const double ky = yfreqs[j];
    for(size_t i = 0; i < mysizes[0]; i++, index++){
      adv[index] += I * ky * buf[index];
    }
  }
  return 0;
}

static int compute_rhs(
    const domain_t * domain,
    const size_t rkstep,
    fluid_t * fluid
){
  const size_t * mysizes = domain->s_x1_mysizes;
  const double * restrict xfreqs = domain->x1_xfreqs;
  const double * restrict yfreqs = domain->x1_yfreqs;
  const fftw_complex * restrict advux = st.s_x1_advs[enum_ux];
  const fftw_complex * restrict advuy = st.s_x1_advs[enum_uy];
  const fftw_complex * restrict advsc = st.s_x1_advs[enum_sc];
  fftw_complex * restrict slopeux = fluid->fields[enum_ux]->s_x1_slopes[rkstep];
  fftw_complex * restrict slopeuy = fluid->fields[enum_uy]->s_x1_slopes[rkstep];
  fftw_complex * restrict slopesc = fluid->fields[enum_sc]->s_x1_slopes[rkstep];
  for(size_t index = 0, j = 0; j < mysizes[1]; j++){
    const double ky = yfreqs[j];
    for(size_t i = 0; i < mysizes[0]; i++, index++){
      const double kx = xfreqs[i];
      if(0 == i && 0 == j) continue;
      const double k2 = + 1. * kx * kx + 1. * ky * ky;
      const fftw_complex ip =
        + 1. * kx * advux[index]
        + 1. * ky * advuy[index];
      slopeux[index] = -1. * advux[index] + kx / k2 * ip;
      slopeuy[index] = -1. * advuy[index] + ky / k2 * ip;
    }
  }
  for(size_t index = 0; index < mysizes[0] * mysizes[1]; index++){
    slopesc[index] = -1. * advsc[index];
  }
  return 0;
}

int compute_slopes(
    const domain_t * domain,
    const size_t rkstep,
    fluid_t * fluid
){
  if(!st.initialised){
    const size_t s_x1_nitems = domain->s_x1_mysizes[0] * domain->s_x1_mysizes[1];
    const size_t p_y1_nitems = domain->p_y1_mysizes[0] * domain->p_y1_mysizes[1];
    st.s_x1_buf = memory_fftw_calloc(s_x1_nitems, sizeof(fftw_complex));
    for(size_t n = 0; n < NDIMS + 1; n++){
      st.s_x1_advs[n] = memory_fftw_calloc(s_x1_nitems, sizeof(fftw_complex));
    }
    st.p_y1_buf = memory_fftw_calloc(p_y1_nitems, sizeof(double));
    st.initialised = true;
  }
  // repeat the same thing for each field 
  for(size_t n = 0; n < NDIMS + 1; n++){
    const double * restrict iarray = fluid->fields[n]->p_y1_array;
    fftw_complex * restrict oarray = st.s_x1_advs[n];
    if(0 != kernel(domain, fluid, iarray, oarray)) return 1;
  }
  // compute right-hand-side non-linear terms using the results by kernel
  if(0 != compute_rhs(domain, rkstep, fluid)) return 1;
  return 0;
}

