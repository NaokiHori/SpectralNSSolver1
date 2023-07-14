#include <stdbool.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"

// internal buffers and plans
typedef struct {
  bool initialised;
  size_t p_glsizes[NDIMS];
  size_t s_glsizes[NDIMS];
  size_t s_x1_mysizes[NDIMS];
  size_t s_y1_mysizes[NDIMS];
#if NDIMS == 3
  size_t s_z1_mysizes[NDIMS];
#endif
#if NDIMS == 2
  size_t p_y1_mysizes[NDIMS];
#else
  size_t p_z1_mysizes[NDIMS];
#endif
  fftw_complex * restrict s_x1_pencil_s;
  fftw_complex * restrict s_x1_pencil_p;
  fftw_complex * restrict s_y1_pencil_s;
#if NDIMS == 2
  double       * restrict p_y1_pencil_p;
#else
  fftw_complex * restrict s_y1_pencil_p;
  fftw_complex * restrict s_z1_pencil_s;
  double       * restrict p_z1_pencil_p;
#endif
  fftw_plan s2p[NDIMS];
  fftw_plan p2s[NDIMS];
  sdecomp_transpose_plan_t * x1_to_y1;
  sdecomp_transpose_plan_t * y1_to_x1;
#if NDIMS == 3
  sdecomp_transpose_plan_t * y1_to_z1;
  sdecomp_transpose_plan_t * z1_to_y1;
#endif
} st_t;
static st_t st = {
  .initialised = false,
};

static int init(
    const domain_t * domain
){
  const sdecomp_info_t * info = domain->info;
  for(size_t dim = 0; dim < NDIMS; dim++){
    st.p_glsizes[dim] = domain->p_glsizes[dim];
    st.s_glsizes[dim] = domain->s_glsizes[dim];
  }
  // local array size, x1 pencil
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, st.s_glsizes[dim], &st.s_x1_mysizes[dim]);
  }
  // local array size, y1 pencil
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, st.s_glsizes[dim], &st.s_y1_mysizes[dim]);
#if NDIMS == 2
    sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, st.p_glsizes[dim], &st.p_y1_mysizes[dim]);
#endif
  }
#if NDIMS == 3
  // local array size, z1 pencil
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, st.s_glsizes[dim], &st.s_z1_mysizes[dim]);
    sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, st.p_glsizes[dim], &st.p_z1_mysizes[dim]);
  }
#endif
  // buffers
  fftw_complex * restrict * s_x1_pencil_s = &st.s_x1_pencil_s;
  fftw_complex * restrict * s_x1_pencil_p = &st.s_x1_pencil_p;
  fftw_complex * restrict * s_y1_pencil_s = &st.s_y1_pencil_s;
#if NDIMS == 2
  double       * restrict * p_y1_pencil_p = &st.p_y1_pencil_p;
#else
  fftw_complex * restrict * s_y1_pencil_p = &st.s_y1_pencil_p;
  fftw_complex * restrict * s_z1_pencil_s = &st.s_z1_pencil_s;
  double       * restrict * p_z1_pencil_p = &st.p_z1_pencil_p;
#endif
#if NDIMS == 2
  const size_t s_x1_pencil_s_nitems = st.s_x1_mysizes[0] * st.s_x1_mysizes[1];
  const size_t s_x1_pencil_p_nitems = st.s_x1_mysizes[0] * st.s_x1_mysizes[1];
  const size_t s_y1_pencil_s_nitems = st.s_y1_mysizes[0] * st.s_y1_mysizes[1];
  const size_t p_y1_pencil_p_nitems = st.p_y1_mysizes[0] * st.p_y1_mysizes[1];
#else
  const size_t s_x1_pencil_s_nitems = st.s_x1_mysizes[0] * st.s_x1_mysizes[1] * st.s_x1_mysizes[2];
  const size_t s_x1_pencil_p_nitems = st.s_x1_mysizes[0] * st.s_x1_mysizes[1] * st.s_x1_mysizes[2];
  const size_t s_y1_pencil_s_nitems = st.s_y1_mysizes[0] * st.s_y1_mysizes[1] * st.s_y1_mysizes[2];
  const size_t s_y1_pencil_p_nitems = st.s_y1_mysizes[0] * st.s_y1_mysizes[1] * st.s_y1_mysizes[2];
  const size_t s_z1_pencil_s_nitems = st.s_z1_mysizes[0] * st.s_z1_mysizes[1] * st.s_z1_mysizes[2];
  const size_t p_z1_pencil_p_nitems = st.p_z1_mysizes[0] * st.p_z1_mysizes[1] * st.p_z1_mysizes[2];
#endif
  *s_x1_pencil_s = memory_fftw_calloc(s_x1_pencil_s_nitems, sizeof(fftw_complex));
  *s_x1_pencil_p = memory_fftw_calloc(s_x1_pencil_p_nitems, sizeof(fftw_complex));
  *s_y1_pencil_s = memory_fftw_calloc(s_y1_pencil_s_nitems, sizeof(fftw_complex));
#if NDIMS == 2
  *p_y1_pencil_p = memory_fftw_calloc(p_y1_pencil_p_nitems, sizeof(      double));
#else
  *s_y1_pencil_p = memory_fftw_calloc(s_y1_pencil_p_nitems, sizeof(fftw_complex));
  *s_z1_pencil_s = memory_fftw_calloc(s_z1_pencil_s_nitems, sizeof(fftw_complex));
  *p_z1_pencil_p = memory_fftw_calloc(p_z1_pencil_p_nitems, sizeof(      double));
#endif
  // fftw plans
  fftw_plan * s2p = st.s2p;
  fftw_plan * p2s = st.p2s;
  // x iDFT
  s2p[0] = fftw_plan_many_dft(
      1, (int [1]){st.p_glsizes[0]},
#if NDIMS == 2
      st.s_x1_mysizes[1],
#else
      st.s_x1_mysizes[1] * st.s_x1_mysizes[2],
#endif
      *s_x1_pencil_s, NULL, 1, st.s_x1_mysizes[0],
      *s_x1_pencil_p, NULL, 1, st.s_x1_mysizes[0],
      FFTW_BACKWARD, FFTW_MEASURE
  );
  // x DFT
  p2s[0] = fftw_plan_many_dft(
      1, (int [1]){st.p_glsizes[0]},
#if NDIMS == 2
      st.s_x1_mysizes[1],
#else
      st.s_x1_mysizes[1] * st.s_x1_mysizes[2],
#endif
      *s_x1_pencil_p, NULL, 1, st.s_x1_mysizes[0],
      *s_x1_pencil_s, NULL, 1, st.s_x1_mysizes[0],
      FFTW_FORWARD, FFTW_MEASURE
  );
#if NDIMS == 2
  // y iRDFT
  s2p[1] = fftw_plan_many_dft_c2r(
      1, (int [1]){st.p_glsizes[1]},
      st.p_y1_mysizes[0],
      *s_y1_pencil_s, NULL, 1, st.s_y1_mysizes[1],
      *p_y1_pencil_p, NULL, 1, st.p_y1_mysizes[1],
      FFTW_MEASURE
  );
  // y RDFT
  p2s[1] = fftw_plan_many_dft_r2c(
      1, (int [1]){st.p_glsizes[1]},
      st.p_y1_mysizes[0],
      *p_y1_pencil_p, NULL, 1, st.p_y1_mysizes[1],
      *s_y1_pencil_s, NULL, 1, st.s_y1_mysizes[1],
      FFTW_MEASURE
  );
#else
  // y iDFT
  s2p[1] = fftw_plan_many_dft(
      1, (int [1]){st.p_glsizes[1]},
      st.s_y1_mysizes[2] * st.s_y1_mysizes[0],
      *s_y1_pencil_s, NULL, 1, st.s_y1_mysizes[1],
      *s_y1_pencil_p, NULL, 1, st.s_y1_mysizes[1],
      FFTW_BACKWARD, FFTW_MEASURE
  );
  // y DFT
  p2s[1] = fftw_plan_many_dft(
      1, (int [1]){st.p_glsizes[1]},
      st.s_y1_mysizes[2] * st.s_y1_mysizes[0],
      *s_y1_pencil_p, NULL, 1, st.s_y1_mysizes[1],
      *s_y1_pencil_s, NULL, 1, st.s_y1_mysizes[1],
      FFTW_FORWARD, FFTW_MEASURE
  );
  // z iRDFT
  s2p[2] = fftw_plan_many_dft_c2r(
      1, (int [1]){st.p_glsizes[2]},
      st.p_z1_mysizes[0] * st.p_z1_mysizes[1],
      *s_z1_pencil_s, NULL, 1, st.s_z1_mysizes[2],
      *p_z1_pencil_p, NULL, 1, st.p_z1_mysizes[2],
      FFTW_MEASURE
  );
  // z RDFT
  p2s[2] = fftw_plan_many_dft_r2c(
      1, (int [1]){st.p_glsizes[2]},
      st.p_z1_mysizes[0] * st.p_z1_mysizes[1],
      *p_z1_pencil_p, NULL, 1, st.p_z1_mysizes[2],
      *s_z1_pencil_s, NULL, 1, st.s_z1_mysizes[2],
      FFTW_MEASURE
  );
#endif
  // pencil rotations
#if NDIMS == 2
  if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, st.s_glsizes, sizeof(fftw_complex), &st.x1_to_y1)){
    printf("x1 to y1 plan creation failed\n");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, st.s_glsizes, sizeof(fftw_complex), &st.y1_to_x1)){
    printf("y1 to x1 plan creation failed\n");
    return 1;
  }
#else
  if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, st.s_glsizes, sizeof(fftw_complex), &st.x1_to_y1)){
    printf("x1 to y1 plan creation failed\n");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, st.s_glsizes, sizeof(fftw_complex), &st.y1_to_x1)){
    printf("y1 to x1 plan creation failed\n");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL, st.s_glsizes, sizeof(fftw_complex), &st.y1_to_z1)){
    printf("y1 to z1 plan creation failed\n");
    return 1;
  }
  if(0 != sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_Y1PENCIL, st.s_glsizes, sizeof(fftw_complex), &st.z1_to_y1)){
    printf("z1 to y1 plan creation failed\n");
    return 1;
  }
#endif
  // update flag
  st.initialised = true;
  return 0;
}

int transform_s2p(
    const domain_t * domain,
    const fftw_complex * restrict bef,
    double * restrict aft
){
  if(!st.initialised){
    if(0 != init(domain)){
      return 1;
    }
  }
  // inverse Fourier transform from spectral domain to physical domain
#if NDIMS == 2
  // copy buffer
  memcpy(st.s_x1_pencil_s, bef, sizeof(fftw_complex) * st.s_x1_mysizes[0] * st.s_x1_mysizes[1]);
  // iFFT in x
  fftw_execute_dft(st.s2p[0], st.s_x1_pencil_s, st.s_x1_pencil_p);
  // rotate x1 pencil to y1 pencil
  sdecomp.transpose.execute(st.x1_to_y1, st.s_x1_pencil_p, st.s_y1_pencil_s);
  // iFFT in y
  fftw_execute_dft_c2r(st.s2p[1], st.s_y1_pencil_s, st.p_y1_pencil_p);
  // normalise FFT
  const size_t * glsizes = st.p_glsizes;
  const size_t * mysizes = st.p_y1_mysizes;
  const double norm = 1. / glsizes[0] / glsizes[1];
  for(size_t index = 0; index < mysizes[0] * mysizes[1]; index++){
    aft[index] = st.p_y1_pencil_p[index] * norm;
  }
#else
  // copy buffer
  memcpy(st.s_x1_pencil_s, bef, sizeof(fftw_complex) * st.s_x1_mysizes[0] * st.s_x1_mysizes[1] * st.s_x1_mysizes[2]);
  // iFFT in x
  fftw_execute_dft(st.s2p[0], st.s_x1_pencil_s, st.s_x1_pencil_p);
  // rotate x1 pencil to y1 pencil
  sdecomp.transpose.execute(st.x1_to_y1, st.s_x1_pencil_p, st.s_y1_pencil_s);
  // iFFT in y
  fftw_execute_dft(st.s2p[1], st.s_y1_pencil_s, st.s_y1_pencil_p);
  // rotate y1 pencil to z1 pencil
  sdecomp.transpose.execute(st.y1_to_z1, st.s_y1_pencil_p, st.s_z1_pencil_s);
  // iFFT in z
  fftw_execute_dft_c2r(st.s2p[2], st.s_z1_pencil_s, st.p_z1_pencil_p);
  // normalise FFT
  const size_t * glsizes = st.p_glsizes;
  const size_t * mysizes = st.p_z1_mysizes;
  const double norm = 1. / glsizes[0] / glsizes[1] / glsizes[2];
  for(size_t index = 0; index < mysizes[0] * mysizes[1] * mysizes[2]; index++){
    aft[index] = st.p_z1_pencil_p[index] * norm;
  }
#endif
  return 0;
}

int transform_p2s(
    const domain_t * domain,
    const double * restrict bef,
    fftw_complex * restrict aft
){
  if(!st.initialised){
    if(0 != init(domain)){
      return 1;
    }
  }
  // Fourier transform from physical domain to spectral domain
#if NDIMS == 2
  // copy buffer
  memcpy(st.p_y1_pencil_p, bef, sizeof(double) * st.p_y1_mysizes[0] * st.p_y1_mysizes[1]);
  // FFT in y
  fftw_execute_dft_r2c(st.p2s[1], st.p_y1_pencil_p, st.s_y1_pencil_s);
  // rotate y1 pencil to x1 pencil
  sdecomp.transpose.execute(st.y1_to_x1, st.s_y1_pencil_s, st.s_x1_pencil_p);
  // FFT in x
  fftw_execute_dft(st.p2s[0], st.s_x1_pencil_p, st.s_x1_pencil_s);
  // copy buffer
  memcpy(aft, st.s_x1_pencil_s, sizeof(fftw_complex) * st.s_x1_mysizes[0] * st.s_x1_mysizes[1]);
#else
  // copy buffer
  memcpy(st.p_z1_pencil_p, bef, sizeof(double) * st.p_z1_mysizes[0] * st.p_z1_mysizes[1] * st.p_z1_mysizes[2]);
  // FFT in z
  fftw_execute_dft_r2c(st.p2s[2], st.p_z1_pencil_p, st.s_z1_pencil_s);
  // rotate z1 pencil to y1 pencil
  sdecomp.transpose.execute(st.z1_to_y1, st.s_z1_pencil_s, st.s_y1_pencil_p);
  // FFT in y
  fftw_execute_dft(st.p2s[1], st.s_y1_pencil_p, st.s_y1_pencil_s);
  // rotate y1 pencil to x1 pencil
  sdecomp.transpose.execute(st.y1_to_x1, st.s_y1_pencil_s, st.s_x1_pencil_p);
  // FFT in x
  fftw_execute_dft(st.p2s[0], st.s_x1_pencil_p, st.s_x1_pencil_s);
  // copy buffer
  memcpy(aft, st.s_x1_pencil_s, sizeof(fftw_complex) * st.s_x1_mysizes[0] * st.s_x1_mysizes[1] * st.s_x1_mysizes[2]);
#endif
  return 0;
}

