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
  if(0 != fileio.r_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), domain->p_glsizes)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), domain->  lengths)){
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
    fileio.w_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), domain->p_glsizes);
    fileio.w_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), domain->  lengths);
  }
  return 0;
}

static int init_wave_numbers(
    domain_t * domain
){
  int * restrict * xwaves = &domain->x1_xwaves;
  int * restrict * ywaves = &domain->x1_ywaves;
  *xwaves = memory_calloc(domain->s_x1_mysizes[0], sizeof(int));
  *ywaves = memory_calloc(domain->s_x1_mysizes[1], sizeof(int));
  int * restrict allwaves[NDIMS] = {
    *xwaves,
    *ywaves,
  };
  // compute wave numbers 
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
  double * restrict * xfreqs = &domain->x1_xfreqs;
  double * restrict * yfreqs = &domain->x1_yfreqs;
  *xfreqs = memory_calloc(domain->s_x1_mysizes[0], sizeof(double));
  *yfreqs = memory_calloc(domain->s_x1_mysizes[1], sizeof(double));
  const int * restrict allwaves[NDIMS] = {
    xwaves,
    ywaves,
  };
  double * restrict allfreqs[NDIMS] = {
    *xfreqs,
    *yfreqs,
  };
  // compute angular frequency 
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
  if(0 != sdecomp.construct(MPI_COMM_WORLD, NDIMS, (size_t [NDIMS]){0, 0}, (bool [NDIMS]){true, true}, &domain->info)){
    printf("%s:%d domain decomposition failed\n", __FILE__, __LINE__);
    return 1;
  }
  // load parameters
  if(0 != load(dirname, domain)){
    return 1;
  }
  // global domain size, spectral domain
  domain->s_glsizes[0] = domain->p_glsizes[0];
  domain->s_glsizes[1] = domain->p_glsizes[1] / 2 + 1;
  // local coordinate
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(domain->info, SDECOMP_X1PENCIL, dim, domain->s_glsizes[dim], &domain->s_x1_mysizes[dim]);
    sdecomp.get_pencil_offset(domain->info, SDECOMP_X1PENCIL, dim, domain->s_glsizes[dim], &domain->s_x1_offsets[dim]);
    sdecomp.get_pencil_mysize(domain->info, SDECOMP_Y1PENCIL, dim, domain->p_glsizes[dim], &domain->p_y1_mysizes[dim]);
    sdecomp.get_pencil_offset(domain->info, SDECOMP_Y1PENCIL, dim, domain->p_glsizes[dim], &domain->p_y1_offsets[dim]);
  }
  init_wave_numbers(domain);
  init_angular_frequency(domain);
  return 0;
}

