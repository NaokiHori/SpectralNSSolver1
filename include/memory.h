#if !defined(MEMORY_H)
#define MEMORY_H

#include <stddef.h> // size_t

extern void * memory_calloc(
    const size_t count,
    const size_t size
);

extern void memory_free(
    void * ptr
);

extern void * memory_fftw_calloc(
    const size_t count,
    const size_t size
);

extern void memory_fftw_free(
    void * ptr
);

#endif // MEMORY_H
