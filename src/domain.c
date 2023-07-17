#include <stdlib.h>
#include <stdbool.h>
#include "memory.h"
#include "sdecomp.h"
#include "domain.h"
#include "fileio.h"

#if !defined(M_PI)
#define M_PI 3.141592653589793238462
#endif

static int load(
    const char dirname[],
    domain_t * domain
){
  if(0 != fileio_r_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, NPY_SZT, sizeof(size_t), domain->p_glsizes)){
    return 1;
  }
  if(0 != fileio_r_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, NPY_DBL, sizeof(double), domain->  lengths)){
    return 1;
  }
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    for(size_t dim = 0; dim < NDIMS; dim++){
      printf("domain->(p_glsizes, lengths)[%zu]: (%5zu, % .2e)\n", dim, domain->p_glsizes[dim], domain->lengths[dim]);
    }
  }
  return 0;
}

int domain_save(
    const char dirname[],
    const domain_t * domain
){
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    fileio_w_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, NPY_SZT, sizeof(size_t), domain->p_glsizes);
    fileio_w_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, NPY_DBL, sizeof(double), domain->  lengths);
  }
  return 0;
}

static int init_wave_numbers(
    domain_t * domain
){
  int * restrict * xwaves = &domain->x1_xwaves;
  int * restrict * ywaves = &domain->x1_ywaves;
#if NDIMS == 3
  int * restrict * zwaves = &domain->x1_zwaves;
#endif
  *xwaves = memory_calloc(domain->s_x1_mysizes[0], sizeof(int));
  *ywaves = memory_calloc(domain->s_x1_mysizes[1], sizeof(int));
#if NDIMS == 3
  *zwaves = memory_calloc(domain->s_x1_mysizes[2], sizeof(int));
#endif
  int * restrict allwaves[NDIMS] = {
    *xwaves,
    *ywaves,
#if NDIMS == 3
    *zwaves,
#endif
  };
  // compute wave numbers | 17
  // wave numbers
  //   0, 1, ..., N/2-1, -N/2, -N/2+1, ..., -2, -1
  // e.g. glsize = 8
  //   -> +0 +1 +2 +3 -4 -3 -2 -1
  for(size_t dim = 0; dim < NDIMS; dim++){
    int * restrict waves = allwaves[dim];
    const int glsize = (int)domain->p_glsizes[dim];
    const int mysize = (int)domain->s_x1_mysizes[dim];
    const int offset = (int)domain->s_x1_offsets[dim];
    for(int n = 0; n < mysize; n++){
      waves[n] = n + offset;
      if(glsize / 2 <= waves[n]){
        // negative wave numbers
        waves[n] -= glsize;
      }
    }
  }
  return 0;
}

static int init_angular_frequency(
    domain_t * domain
){
  int * restrict xwaves = domain->x1_xwaves;
  int * restrict ywaves = domain->x1_ywaves;
#if NDIMS == 3
  int * restrict zwaves = domain->x1_zwaves;
#endif
  double * restrict * xfreqs = &domain->x1_xfreqs;
  double * restrict * yfreqs = &domain->x1_yfreqs;
#if NDIMS == 3
  double * restrict * zfreqs = &domain->x1_zfreqs;
#endif
  *xfreqs = memory_calloc(domain->s_x1_mysizes[0], sizeof(double));
  *yfreqs = memory_calloc(domain->s_x1_mysizes[1], sizeof(double));
#if NDIMS == 3
  *zfreqs = memory_calloc(domain->s_x1_mysizes[2], sizeof(double));
#endif
  const int * restrict allwaves[NDIMS] = {
    xwaves,
    ywaves,
#if NDIMS == 3
    zwaves,
#endif
  };
  double * restrict allfreqs[NDIMS] = {
    *xfreqs,
    *yfreqs,
#if NDIMS == 3
    *zfreqs,
#endif
  };
  // compute angular frequency | 9
  for(size_t dim = 0; dim < NDIMS; dim++){
    const double length = domain->lengths[dim];
    const size_t mysize = domain->s_x1_mysizes[dim];
    const int * restrict waves = allwaves[dim];
    double * restrict freqs = allfreqs[dim];
    for(size_t n = 0; n < mysize; n++){
      freqs[n] = 2. * M_PI / length * waves[n];
    }
  }
  return 0;
}

int domain_init(
    const char dirname[],
    domain_t * domain
){
  // decompose domain
#if NDIMS == 2
  if(0 != sdecomp.construct(MPI_COMM_WORLD, NDIMS, (size_t [NDIMS]){0, 0}, (bool [NDIMS]){true, true}, &domain->info)){
    printf("%s:%d domain decomposition failed\n", __FILE__, __LINE__);
    return 1;
  }
#else
  if(0 != sdecomp.construct(MPI_COMM_WORLD, NDIMS, (size_t [NDIMS]){0, 0, 0}, (bool [NDIMS]){true, true, true}, &domain->info)){
    printf("%s:%d domain decomposition failed\n", __FILE__, __LINE__);
    return 1;
  }
#endif
  // load parameters
  if(0 != load(dirname, domain)){
    return 1;
  }
  // global domain size, spectral domain
#if NDIMS == 2
  domain->s_glsizes[0] = domain->p_glsizes[0];
  domain->s_glsizes[1] = domain->p_glsizes[1] / 2 + 1;
#else
  domain->s_glsizes[0] = domain->p_glsizes[0];
  domain->s_glsizes[1] = domain->p_glsizes[1];
  domain->s_glsizes[2] = domain->p_glsizes[2] / 2 + 1;
#endif
  // local coordinate
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(domain->info, SDECOMP_X1PENCIL, dim, domain->s_glsizes[dim], &domain->s_x1_mysizes[dim]);
    sdecomp.get_pencil_offset(domain->info, SDECOMP_X1PENCIL, dim, domain->s_glsizes[dim], &domain->s_x1_offsets[dim]);
#if NDIMS == 2
    sdecomp.get_pencil_mysize(domain->info, SDECOMP_Y1PENCIL, dim, domain->p_glsizes[dim], &domain->p_y1_mysizes[dim]);
    sdecomp.get_pencil_offset(domain->info, SDECOMP_Y1PENCIL, dim, domain->p_glsizes[dim], &domain->p_y1_offsets[dim]);
#else
    sdecomp.get_pencil_mysize(domain->info, SDECOMP_Z1PENCIL, dim, domain->p_glsizes[dim], &domain->p_z1_mysizes[dim]);
    sdecomp.get_pencil_offset(domain->info, SDECOMP_Z1PENCIL, dim, domain->p_glsizes[dim], &domain->p_z1_offsets[dim]);
#endif
  }
  init_wave_numbers(domain);
  init_angular_frequency(domain);
  return 0;
}

