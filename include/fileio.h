#if !defined(FILEIO_H)
#define FILEIO_H

#include <stdio.h> // FILE, size_t
#include <mpi.h>   // MPI_Datatype

// general-purpose file opener
extern FILE * fileio_fopen(
    const char * path,
    const char * mode
);

// general-purpose file closer
extern int fileio_fclose(
    FILE * stream
);

// prepare directory to be stored
extern int fileio_mkdir(
    const char dirname[]
);

// NPY datatypes, which are embedded in NPY files ("dtype" argument)
// they are declared here and defined in src/fileio/entrypoint.c
// 1-byte boolean
extern const char NPY_BOL[];
// 4-byte little-endian integer
extern const char NPY_INT[];
// 8-byte little-endian floating point
extern const char NPY_DBL[];
// 8-byte little-endian unsigned long (need check)
extern const char NPY_SZT[];
// 16-byte little-endian double complex
extern const char NPY_CMP[];

// NPY serial read (called by one process)
extern int fileio_r_serial(
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    void * data
);

// NPY serial write (called by one process)
extern int fileio_w_serial(
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    const void * data
);

// NPY parallel read of N-dimensional array (called by all processes)
extern int fileio_r_nd_parallel(
    const MPI_Comm comm,
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const int * array_of_sizes,
    const int * array_of_subsizes,
    const int * array_of_starts,
    const char dtype[],
    const size_t size,
    void * data
);

// NPY parallel write of N-dimensional array (called by all processes)
extern int fileio_w_nd_parallel(
    const MPI_Comm comm,
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const int * array_of_sizes,
    const int * array_of_subsizes,
    const int * array_of_starts,
    const char dtype[],
    const size_t size,
    const void * data
);

#endif // FILEIO_H
