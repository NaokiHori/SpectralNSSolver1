#include "domain.h"
#include "fluid.h"
#include "fileio.h"

int fluid_load(
    const char dirname[],
    const domain_t * domain,
    fluid_t * fluid
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
#if NDIMS == 2
  const int glsizes[NDIMS] = {domain->   s_glsizes[1], domain->   s_glsizes[0]};
  const int mysizes[NDIMS] = {domain->s_x1_mysizes[1], domain->s_x1_mysizes[0]};
  const int offsets[NDIMS] = {domain->s_x1_offsets[1], domain->s_x1_offsets[0]};
#else
  const int glsizes[NDIMS] = {domain->   s_glsizes[2], domain->   s_glsizes[1], domain->   s_glsizes[0]};
  const int mysizes[NDIMS] = {domain->s_x1_mysizes[2], domain->s_x1_mysizes[1], domain->s_x1_mysizes[0]};
  const int offsets[NDIMS] = {domain->s_x1_offsets[2], domain->s_x1_offsets[1], domain->s_x1_offsets[0]};
#endif
  void * arrays[] = {
    fluid->fields[enum_ux]->s_x1_array,
    fluid->fields[enum_uy]->s_x1_array,
#if NDIMS == 3
    fluid->fields[enum_uz]->s_x1_array,
#endif
    fluid->fields[enum_sc]->s_x1_array,
  };
  const char * dsetnames[] = {
    "ux",
    "uy",
#if NDIMS == 3
    "uz",
#endif
    "sc",
  };
  for(size_t index = 0; index < sizeof(arrays) / sizeof(arrays[0]); index++){
    if(0 != fileio_r_nd_parallel(
          comm_cart,
          dirname,
          dsetnames[index],
          NDIMS,
          glsizes,
          mysizes,
          offsets,
          NPY_CMP,
          sizeof(fftw_complex),
          arrays[index]
    )) return 1;
  }
  return 0;
}

int fluid_save(
    const char dirname[],
    const domain_t * domain,
    const fluid_t * fluid
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
#if NDIMS == 2
  const int glsizes[NDIMS] = {domain->   s_glsizes[1], domain->   s_glsizes[0]};
  const int mysizes[NDIMS] = {domain->s_x1_mysizes[1], domain->s_x1_mysizes[0]};
  const int offsets[NDIMS] = {domain->s_x1_offsets[1], domain->s_x1_offsets[0]};
#else
  const int glsizes[NDIMS] = {domain->   s_glsizes[2], domain->   s_glsizes[1], domain->   s_glsizes[0]};
  const int mysizes[NDIMS] = {domain->s_x1_mysizes[2], domain->s_x1_mysizes[1], domain->s_x1_mysizes[0]};
  const int offsets[NDIMS] = {domain->s_x1_offsets[2], domain->s_x1_offsets[1], domain->s_x1_offsets[0]};
#endif
  const void * arrays[] = {
    fluid->fields[enum_ux]->s_x1_array,
    fluid->fields[enum_uy]->s_x1_array,
#if NDIMS == 3
    fluid->fields[enum_uz]->s_x1_array,
#endif
    fluid->fields[enum_sc]->s_x1_array,
  };
  const char * dsetnames[] = {
    "ux",
    "uy",
#if NDIMS == 3
    "uz",
#endif
    "sc",
  };
  for(size_t index = 0; index < sizeof(arrays) / sizeof(arrays[0]); index++){
    if(0 != fileio_w_nd_parallel(
        comm_cart,
        dirname,
        dsetnames[index],
        NDIMS,
        glsizes,
        mysizes,
        offsets,
        NPY_CMP,
        sizeof(fftw_complex),
        arrays[index]
    )) return 1;
  }
  return 0;
}

