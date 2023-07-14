#include <stdio.h>
#include <stdbool.h>
#include "memory.h"
#include "fileio.h"
#define FILEIO_INTERNAL
#include "internal.h"

/**
 * @brief read data from a npy file, by one process
 * @param[in]  dirname  : name of directory in which a target npy file is contained
 * @param[in]  dsetname : name of dataset
 * @param[in]  ndims    : number of dimensions of dataset
 * @param[in]  shape    : shape of dataset
 * @param[in]  dtype    : datatype, e.g. '<f8'
 * @param[in]  size     : size of each element
 * @param[out] data     : pointer to the data to be loaded
 */
int fileio_r_serial(
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    void * data
){
  char * fname = fileio_internal_create_npyfname(dirname, dsetname);
  const size_t header_size = fileio_internal_r_npy_header(fname, ndims, shape, dtype, false);
  if(0 == header_size){
    fprintf(stderr, "%s: NPY header load failed\n", fname);
    memory_free(fname);
    return 1;
  }
  FILE * fp = fileio_fopen(fname, "r");
  if(NULL == fp){
    memory_free(fname);
    return 1;
  }
  if(0 != fseek(fp, (long)header_size, SEEK_SET)){
    fprintf(stderr, "%s: fseek failed\n", fname);
    fileio_fclose(fp);
    memory_free(fname);
    return 1;
  }
  size_t nitems = 1;
  for(size_t dim = 0; dim < ndims; dim++){
    nitems *= shape[dim];
  }
  const size_t nitems_ = fread(data, size, nitems, fp);
  if(nitems_ != nitems){
    fprintf(stderr, "%s: fread failed\n", fname);
    fileio_fclose(fp);
    memory_free(fname);
    return 1;
  }
  fileio_fclose(fp);
  memory_free(fname);
  return 0;
}

/**
 * @brief write data to a npy file, by one process
 * @param[in] dirname  : name of directory in which a target npy file is contained
 * @param[in] dsetname : name of dataset
 * @param[in] ndims    : number of dimensions of dataset
 * @param[in] shape    : shape of dataset
 * @param[in] dtype    : datatype, e.g. '<f8'
 * @param[in] size     : size of each element
 * @param[in] data     : pointer to the data to be written
 */
int fileio_w_serial(
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    const void * data
){
  char * fname = fileio_internal_create_npyfname(dirname, dsetname);
  const size_t header_size = fileio_internal_w_npy_header(fname, ndims, shape, dtype, false);
  if(0 == header_size){
    memory_free(fname);
    return 1;
  }
  FILE * fp = fileio_fopen(fname, "a");
  if(NULL == fp){
    memory_free(fname);
    return 1;
  }
  if(0 != fseek(fp, (long)header_size, SEEK_SET)){
    fprintf(stderr, "%s: fseek failed\n", fname);
    fileio_fclose(fp);
    memory_free(fname);
    return 1;
  }
  size_t nitems = 1;
  for(size_t dim = 0; dim < ndims; dim++){
    nitems *= shape[dim];
  }
  const size_t nitems_ = fwrite(data, size, nitems, fp);
  if(nitems_ != nitems){
    fprintf(stderr, "%s: fwrite failed\n", fname);
    fileio_fclose(fp);
    memory_free(fname);
    return 1;
  }
  fileio_fclose(fp);
  memory_free(fname);
  return 0;
}

