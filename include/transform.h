#if !defined(TRANSFORM_H)
#define TRANSFORM_H

#include <fftw3.h>
#include "domain.h"

extern int transform_s2p(
    const domain_t * domain,
    const fftw_complex * restrict bef,
    double * restrict aft
);

extern int transform_p2s(
    const domain_t * domain,
    const double * restrict bef,
    fftw_complex * restrict aft
);

#endif // TRANSFORM_H
