#include <stdio.h>
#include <mpi.h>
#include "config.h"
#include "timer.h"
#include "domain.h"
#include "fluid.h"
#include "logging.h"
#include "save.h"
#include "fileio.h"

static int save_entrypoint(
    const domain_t * domain,
    const int step,
    const double time,
    const fluid_t * fluid
){
  char * dirname = NULL;
  save.prepare(domain, step, &dirname);
  fileio_w_serial(dirname, "step", 0, NULL, NPY_INT, sizeof(   int), &step);
  fileio_w_serial(dirname, "time", 0, NULL, NPY_DBL, sizeof(double), &time);
  domain_save(dirname, domain);
  fluid_save(dirname, domain, fluid);
  return 0;
}

int main(
    int argc,
    char * argv[]
){
  // launch MPI
  MPI_Init(NULL, NULL);
  // check my rank to dump logs only from the main process
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const double tic = timer();
  // check name of the initial velocity field is given
  if(2 != argc){
    if(0 == myrank) printf("give directory name: ./a.out <name of directory>\n");
    goto abort;
  }
  const char * dirname_ic = argv[1];
  // initialise structure, domain_t
  domain_t domain = {0};
  if(0 != domain_init(dirname_ic, &domain)){
    goto abort;
  }
  // initialise structure, fluid_t
  fluid_t fluid = {0};
  if(0 != fluid_init(dirname_ic, &domain, &fluid)){
    goto abort;
  }
  // load conditions to terminate the solver from environment variables
  double  timemax = 0.;
  double wtimemax = 0.;
  if(0 != config.get_double("timemax", &timemax)){
    goto abort;
  }
  if(0 != config.get_double("wtimemax", &wtimemax)){
    goto abort;
  }
  // load current time step and simulation time units
  size_t step = 0;
  double time = 0.;
  if(0 != fileio_r_serial(dirname_ic, "step", 0, NULL, NPY_SZT, sizeof(size_t), &step)){
    goto abort;
  }
  if(0 != fileio_r_serial(dirname_ic, "time", 0, NULL, NPY_DBL, sizeof(double), &time)){
    goto abort;
  }
  if(0 == myrank) printf("start from: step %zu, time % .7e\n", step, time);
  // initialise logger
  if(0 != logging.init(&domain, time)){
    goto abort;
  }
  // initialise flow field saver
  if(0 != save.init(&domain, time)){
    goto abort;
  }
  // main loop to integrate NS equations in time
  for(double dt = 1.; ; ){
    // integrate the flow field in time
    if(0 != fluid_integrate(&domain, &fluid, &dt)){
      goto abort;
    }
    // now flow field is updated, increment counter and time
    time += dt;
    step += 1;
    const double toc = timer();
    // terminate if the simulation is done
    if(time > timemax){
      break;
    }
    // terminate if the maximum duration is reached
    if(toc - tic > wtimemax){
      break;
    }
    // dump log files regulary
    if(logging.get_next_time() < time){
      logging.check_and_output(&domain, step, time, dt, toc - tic, &fluid);
    }
    // save flow fields regulary
    if(save.get_next_time() < time){
      save_entrypoint(&domain, step, time, &fluid);
    }
  }
  // save last field
  save_entrypoint(&domain, step, time, &fluid);
abort:
  MPI_Finalize();
  return 0;
}

