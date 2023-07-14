#if !defined(FLUID_INTERNAL_H)
#define FLUID_INTERNAL_H

#if !defined(FLUID_INTERNAL)
#error "do not include this header"
#endif

extern int compute_physical_fields(
    const domain_t * domain,
    fluid_t * fluid
);

extern int decide_dt(
    const domain_t * domain,
    const fluid_t * fluid,
    double * restrict dt
);

extern int compute_slopes(
    const domain_t * domain,
    const size_t rkstep,
    fluid_t * fluid
);

extern int update_fields(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#endif // FLUID_INTERNAL_H
