#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <fftw3.h>
#include "memory.h"

void * memory_calloc(
    const size_t count,
    const size_t size
){
  void * ptr = calloc(count, size);
  if(NULL == ptr){
    fprintf(stderr, "memory allocation error\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return ptr;
}

void memory_free(
    void * ptr
){
  free(ptr);
}

void * memory_fftw_calloc(
    const size_t count,
    const size_t size
){
  void * ptr = fftw_malloc(count * size);
  if(NULL == ptr){
    fprintf(stderr, "memory allocation error\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  memset(ptr, 0, count * size);
  return ptr;
}

void memory_fftw_free(
    void * ptr
){
  fftw_free(ptr);
}

