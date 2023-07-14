#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "sdecomp.h"
#include "memory.h"
#include "config.h"
#include "domain.h"
#include "save.h"
#include "fileio.h"

// parameters deciding directory name
static const char dirname_prefix[] = {"output/save/step"};
static const int dirname_ndigits = 10;

// name of directory
static char * g_dirname = NULL;
static size_t g_dirname_nchars = 0;

// scheduler
static double g_rate = DBL_MAX;
static double g_next = 0.;

/**
 * @brief constructor - schedule saving flow fields
 * @param[in] domain : MPI communicator
 * @param[in] time   : current time (hereafter in free-fall time units)
 */
static int init(
    const domain_t * domain,
    const double time
){
  if(0 != config.get_double("save_rate", &g_rate)){
    return 1;
  }
  // schedule next event
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, time) / g_rate
  );
  // allocate directory name
  g_dirname_nchars =
    + strlen(dirname_prefix)
    + dirname_ndigits;
  g_dirname = memory_calloc(g_dirname_nchars + 2, sizeof(char));
  // report
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(myrank == 0){
    printf("SAVE\n");
    printf("\tnext: % .3e\n", g_next);
    printf("\trate: % .3e\n", g_rate);
    fflush(stdout);
  }
  return 0;
}

/**
 * @brief prepare place to output flow fields and save auxiliary data
 * @param[in]  domain  : information related to MPI domain decomposition
 * @param[in]  step    : time step
 * @param[out] dirname : name of created directory
 */
static int prepare(
    const domain_t * domain,
    const size_t step,
    char ** dirname
){
  // set directory name
  snprintf(g_dirname, g_dirname_nchars + 1, "%s%0*zu", dirname_prefix, dirname_ndigits, step);
  *dirname = g_dirname;
  // get communicator to identify the main process
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // create directory
  if(0 == myrank){
    // although it may fail, anyway continue, which is designed to be safe
    fileio_mkdir(g_dirname);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // schedule next saving event
  g_next += g_rate;
  return 0;
}

/**
 * @brief getter of a member: g_next
 * @return : g_next
 */
static double get_next_time(
    void
){
  return g_next;
}

const save_t save = {
  .init          = init,
  .prepare       = prepare,
  .get_next_time = get_next_time,
};

