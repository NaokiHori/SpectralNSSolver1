#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "sdecomp.h"
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#define FLUID_INTERNAL
#include "internal.h"

int decide_dt(
    const domain_t * domain,
    const fluid_t * fluid,
    double * dt
){
  // z1 pencil
  const size_t * mysizes = domain->p_z1_mysizes;
  const double dx = domain->lengths[0] / domain->p_glsizes[0];
  const double dy = domain->lengths[1] / domain->p_glsizes[1];
  const double dz = domain->lengths[2] / domain->p_glsizes[2];
  // use physical velocity
  const double * restrict ux = fluid->fields[enum_ux]->p_z1_array;
  const double * restrict uy = fluid->fields[enum_uy]->p_z1_array;
  const double * restrict uz = fluid->fields[enum_uz]->p_z1_array;
  // val: local (physical) velocity / grid size
  // check maximum value in my range
  double maxval = 0.;
  const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
  for(size_t index = 0; index < nitems; index++){
    double val = 0.;
    val += fabs(ux[index]) / dx;
    val += fabs(uy[index]) / dy;
    val += fabs(uz[index]) / dz;
    maxval = fmax(maxval, val);
  }
  // communicate maximum value among all pencils
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &maxval, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  // multiply safety factor to decide the time step size
  const double dt_adv = runge_kutta_cfl / NDIMS / maxval;
  // decide time step size
  // NOTE: limit maximum to avoid abrupt change
  *dt = *dt * 1.2;
  *dt = fmin(*dt, dt_adv);
  return 0;
}

