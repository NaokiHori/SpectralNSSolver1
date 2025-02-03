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
  const int glsizes[NDIMS] = {domain->   s_glsizes[1], domain->   s_glsizes[0]};
  const int mysizes[NDIMS] = {domain->s_x1_mysizes[1], domain->s_x1_mysizes[0]};
  const int offsets[NDIMS] = {domain->s_x1_offsets[1], domain->s_x1_offsets[0]};
  void * arrays[] = {
    fluid->fields[enum_ux]->s_x1_array,
    fluid->fields[enum_uy]->s_x1_array,
    fluid->fields[enum_sc]->s_x1_array,
  };
  const char * dsetnames[] = {
    "ux",
    "uy",
    "sc",
  };
  for(size_t index = 0; index < sizeof(arrays) / sizeof(arrays[0]); index++){
    if(0 != fileio.r_nd_parallel(
          comm_cart,
          dirname,
          dsetnames[index],
          NDIMS,
          glsizes,
          mysizes,
          offsets,
          fileio.npy_complex,
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
  const int glsizes[NDIMS] = {domain->   s_glsizes[1], domain->   s_glsizes[0]};
  const int mysizes[NDIMS] = {domain->s_x1_mysizes[1], domain->s_x1_mysizes[0]};
  const int offsets[NDIMS] = {domain->s_x1_offsets[1], domain->s_x1_offsets[0]};
  const void * arrays[] = {
    fluid->fields[enum_ux]->s_x1_array,
    fluid->fields[enum_uy]->s_x1_array,
    fluid->fields[enum_sc]->s_x1_array,
  };
  const char * dsetnames[] = {
    "ux",
    "uy",
    "sc",
  };
  for(size_t index = 0; index < sizeof(arrays) / sizeof(arrays[0]); index++){
    if(0 != fileio.w_nd_parallel(
        comm_cart,
        dirname,
        dsetnames[index],
        NDIMS,
        glsizes,
        mysizes,
        offsets,
        fileio.npy_complex,
        sizeof(fftw_complex),
        arrays[index]
    )) return 1;
  }
  return 0;
}

