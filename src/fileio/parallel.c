#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <mpi.h>
#include "memory.h"
#include "fileio.h"
#define FILEIO_INTERNAL
#include "internal.h"

static int mpi_file_open(
    const MPI_Comm comm,
    char * fname,
    int amode,
    MPI_File * fh
){
  const int mpi_error_code = MPI_File_open(comm, fname, amode, MPI_INFO_NULL, fh);
  if(MPI_SUCCESS != mpi_error_code){
    char string[MPI_MAX_ERROR_STRING] = {'\0'};
    int resultlen = 0;
    MPI_Error_string(mpi_error_code, string, &resultlen);
    fprintf(stderr, "%s: %s\n", fname, string);
    memory_free(fname);
    return 1;
  }
  return 0;
}

static int get_count(
    const size_t ndims,
    const int * mysizes
){
  int count = 1;
  for(size_t n = 0; n < ndims; n++){
    count *= mysizes[n];
  }
  return count;
}

static int prepare_view(
    const size_t ndims,
    const int * glsizes,
    const int * mysizes,
    const int * offsets,
    MPI_File fh,
    const size_t header_size,
    const MPI_Datatype basetype,
    MPI_Datatype * filetype
){
  // create data type and set file view
  MPI_Type_create_subarray((int)ndims, glsizes, mysizes, offsets, MPI_ORDER_C, basetype, filetype);
  MPI_Type_commit(filetype);
  MPI_File_set_view(fh, (MPI_Offset)header_size, basetype, *filetype, "native", MPI_INFO_NULL);
  return 0;
}

static int destroy_view(
    MPI_Datatype * filetype
){
  // clean-up datatype
  MPI_Type_free(filetype);
  return 0;
}

/**
 * @brief read N-dimensional data from a npy file, by all processes
 * @param[in]  comm     : communicator to which all processes calling this function belong
 * @param[in]  dirname  : name of directory in which a target npy file is contained
 * @param[in]  dsetname : name of dataset
 * @param[in]  ndims    : number of dimensions of the array
 * @param[in]  glsizes  : global sizes   of the dataset
 * @param[in]  mysizes  : local  sizes   of the dataset
 * @param[in]  offsets  : local  offsets of the dataset
 * @param[in]  dtype    : NPY data type
 * @param[in]  size     : size of each element
 * @param[out] data     : pointer to the data to be loaded
 */
int fileio_r_nd_parallel(
    const MPI_Comm comm,
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const int * glsizes,
    const int * mysizes,
    const int * offsets,
    const char dtype[],
    const size_t size,
    void * data
){
  const int root = 0;
  int myrank = 0;
  MPI_Comm_rank(comm, &myrank);
  char * fname = fileio_internal_create_npyfname(dirname, dsetname);
  // check header by main process
  size_t header_size = 0;
  if(root == myrank){
    // set values which are expected to be in NPY file
    size_t * shape = memory_calloc(ndims, sizeof(size_t));
    for(size_t dim = 0; dim < ndims; dim++){
      shape[dim] = (size_t)glsizes[dim];
    }
    header_size = fileio_internal_r_npy_header(fname, ndims, shape, dtype, false);
    memory_free(shape);
  }
  // share result
  MPI_Bcast(&header_size, sizeof(size_t) / sizeof(uint8_t), MPI_BYTE, root, comm);
  if(0 == header_size){
    memory_free(fname);
    return 1;
  }
  // open file
  MPI_File fh = NULL;
  if(0 != mpi_file_open(comm, fname, MPI_MODE_RDONLY, &fh)){
    return 1;
  }
  // prepare file view
  MPI_Datatype basetype = MPI_BYTE;
  MPI_Datatype filetype = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(size, basetype, &basetype);
  MPI_Type_commit(&basetype);
  prepare_view((int)ndims, glsizes, mysizes, offsets, fh, header_size, basetype, &filetype);
  // get number of elements which are locally read
  const int count = get_count(ndims, mysizes);
  // read
  MPI_File_read_all(fh, data, count, basetype, MPI_STATUS_IGNORE);
  // clean-up file view
  MPI_Type_free(&basetype);
  destroy_view(&filetype);
  // close file
  MPI_File_close(&fh);
  memory_free(fname);
  return 0;
}

/**
 * @brief write N-dimensional data to a npy file, by all processes
 * @param[in] comm     : communicator to which all processes calling this function belong
 * @param[in] dirname  : name of directory in which a target npy file is contained
 * @param[in] dsetname : name of dataset
 * @param[in] ndims    : number of dimensions of the array
 * @param[in] glsizes  : global sizes   of the dataset
 * @param[in] mysizes  : local  sizes   of the dataset
 * @param[in] offsets  : local  offsets of the dataset
 * @param[in] dtype    : NPY data type
 * @param[in] size     : size of each element
 * @param[in] data     : pointer to the data to be written
 */
int fileio_w_nd_parallel(
    const MPI_Comm comm,
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const int * glsizes,
    const int * mysizes,
    const int * offsets,
    const char dtype[],
    const size_t size,
    const void * data
){
  const int root = 0;
  int myrank = 0;
  MPI_Comm_rank(comm, &myrank);
  char * fname = fileio_internal_create_npyfname(dirname, dsetname);
  // check header by main process
  size_t header_size = 0;
  if(root == myrank){
    // set values which are expected to be in NPY file
    size_t * shape = memory_calloc(ndims, sizeof(size_t));
    for(size_t dim = 0; dim < ndims; dim++){
      shape[dim] = (size_t)glsizes[dim];
    }
    header_size = fileio_internal_w_npy_header(fname, ndims, shape, dtype, false);
    memory_free(shape);
  }
  // share result
  MPI_Bcast(&header_size, sizeof(size_t) / sizeof(uint8_t), MPI_BYTE, root, comm);
  if(0 == header_size){
    memory_free(fname);
    return 1;
  }
  // open file
  MPI_File fh = NULL;
  if(0 != mpi_file_open(comm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, &fh)){
    return 1;
  }
  // prepare file view
  MPI_Datatype basetype = MPI_BYTE;
  MPI_Datatype filetype = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(size, basetype, &basetype);
  MPI_Type_commit(&basetype);
  prepare_view((int)ndims, glsizes, mysizes, offsets, fh, header_size, basetype, &filetype);
  // get number of elements which are locally written
  const int count = get_count(ndims, mysizes);
  // write
  MPI_File_write_all(fh, data, count, basetype, MPI_STATUS_IGNORE);
  // clean-up file view
  MPI_Type_free(&basetype);
  destroy_view(&filetype);
  // close file
  MPI_File_close(&fh);
  memory_free(fname);
  return 0;
}

