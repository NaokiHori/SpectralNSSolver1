#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <fftw3.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "logging.h"

static double g_rate = DBL_MAX;
static double g_next = 0.;

static int init(
    const domain_t * domain,
    const double time
){
  if(0 != config.get_double("log_rate", &g_rate)){
    return 1;
  }
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, time) / g_rate
  );
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    printf("LOGGING\n");
    printf("\tnext: % .3e\n", g_next);
    printf("\trate: % .3e\n", g_rate);
    fflush(stdout);
  }
  return 0;
}

static void show_progress(
    const char fname[],
    const domain_t * domain,
    const double time,
    const size_t step,
    const double dt,
    const double wtime
){
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    FILE * const fp = fileio.fopen(fname, "a");
    if(NULL != fp){
      // show progress to standard output and file
      // output to stdout and file
#define MPRINT(...) { \
      fprintf(fp,     __VA_ARGS__); \
      fprintf(stdout, __VA_ARGS__); \
}
      MPRINT("step %zu, time %.1f, dt %.2e, elapsed %.1f [sec]\n", step, time, dt, wtime);
#undef MPRINT
      fileio.fclose(fp);
    }
  }
}

static void check_divergence(
    const char fname[],
    const domain_t * domain,
    const double time,
    const size_t step,
    const fluid_t * fluid
){
  const size_t * mysizes = domain->s_x1_mysizes;
  const double * restrict xfreqs = domain->x1_xfreqs;
  const double * restrict yfreqs = domain->x1_yfreqs;
  const double * restrict zfreqs = domain->x1_zfreqs;
  const fftw_complex * restrict ux = fluid->fields[enum_ux]->s_x1_array;
  const fftw_complex * restrict uy = fluid->fields[enum_uy]->s_x1_array;
  const fftw_complex * restrict uz = fluid->fields[enum_uz]->s_x1_array;
  double maxdiv = 0.;
  for(size_t index = 0, k = 0; k < mysizes[2]; k++){
    const double kz = zfreqs[k];
    for(size_t j = 0; j < mysizes[1]; j++){
      const double ky = yfreqs[j];
      for(size_t i = 0; i < mysizes[0]; i++, index++){
        const double kx = xfreqs[i];
        const double div = cabs(
            + I * kx * ux[index]
            + I * ky * uy[index]
            + I * kz * uz[index]
        );
        maxdiv = fmax(maxdiv, div);
      }
    }
  }
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &maxdiv, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    FILE * const fp = fileio.fopen(fname, "a");
    if(NULL != fp){
      fprintf(fp, "%10zu % 8.2e % .1e\n", step, time, maxdiv);
      fileio.fclose(fp);
    }
  }
}

static void check_extrema(
    const char fname[],
    const domain_t * domain,
    const double time,
    const size_t step,
    const fluid_t * fluid
){
  const size_t * mysizes = domain->p_z1_mysizes;
  const double * restrict ux = fluid->fields[enum_ux]->p_z1_array;
  const double * restrict uy = fluid->fields[enum_uy]->p_z1_array;
  const double * restrict uz = fluid->fields[enum_uz]->p_z1_array;
  const double * restrict sc = fluid->fields[enum_sc]->p_z1_array;
  double maxvals[NDIMS + 1] = {0.};
  const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
  for(size_t index = 0; index < nitems; index++){
    const double vals[NDIMS + 1] = {
      fabs(ux[index]),
      fabs(uy[index]),
      fabs(uz[index]),
      fabs(sc[index]),
    };
    for(size_t dim = 0; dim < NDIMS + 1; dim++){
      maxvals[dim] = fmax(maxvals[dim], vals[dim]);
    }
  }
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, maxvals, NDIMS + 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    FILE * const fp = fileio.fopen(fname, "a");
    if(NULL != fp){
      fprintf(fp, "%10zu % 8.2e ", step, time);
      for(size_t dim = 0; dim < NDIMS + 1; dim++){
        const char del = NDIMS == dim ? '\n' : ' ';
        fprintf(fp, "% .4e%c", maxvals[dim], del);
      }
      fileio.fclose(fp);
    }
  }
}

static void check_energy(
    const char fname[],
    const domain_t * domain,
    const double time,
    const size_t step,
    const fluid_t * fluid
){
  const size_t * mysizes = domain->p_z1_mysizes;
  const double * restrict ux = fluid->fields[enum_ux]->p_z1_array;
  const double * restrict uy = fluid->fields[enum_uy]->p_z1_array;
  const double * restrict uz = fluid->fields[enum_uz]->p_z1_array;
  const double * restrict sc = fluid->fields[enum_sc]->p_z1_array;
  const double cellsize = 1.
    * domain->lengths[0] / domain->p_glsizes[0]
    * domain->lengths[1] / domain->p_glsizes[1]
    * domain->lengths[2] / domain->p_glsizes[2];
  double vals[2] = {0.};
  const size_t nitems = mysizes[0] * mysizes[1] * mysizes[2];
  for(size_t index = 0; index < nitems; index++){
    vals[0] += 0.5 * ux[index] * ux[index] * cellsize;
    vals[0] += 0.5 * uy[index] * uy[index] * cellsize;
    vals[0] += 0.5 * uz[index] * uz[index] * cellsize;
    vals[1] += 0.5 * sc[index] * sc[index] * cellsize;
  }
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_DOUBLE, MPI_SUM, comm_cart);
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL != fp){
      fprintf(fp, "%10zu % 8.2e ", step, time);
      fprintf(fp, "% .15e % .15e\n", vals[0], vals[1]);
      fileio.fclose(fp);
    }
  }
}

static void check_and_output(
    const domain_t * domain,
    const size_t step,
    const double time,
    const double dt,
    const double wtime,
    const fluid_t * fluid
){
  show_progress   ("output/log/progress.dat",   domain, time, step, dt, wtime);
  check_divergence("output/log/divergence.dat", domain, time, step, fluid);
  check_extrema   ("output/log/extrema.dat",    domain, time, step, fluid);
  check_energy    ("output/log/energy.dat",     domain, time, step, fluid);
  g_next += g_rate;
}

static double get_next_time(
    void
){
  return g_next;
}

const logging_t logging = {
  .init             = init,
  .check_and_output = check_and_output,
  .get_next_time    = get_next_time,
};

