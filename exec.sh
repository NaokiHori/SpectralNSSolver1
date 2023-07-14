#!/bin/bash

## temporal information
# maximum duration (in free-fall time)
export timemax=3.0e+1
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=5.0e-1
# save rate (in free-fall time)
export save_rate=1.0e+0

## physical parameters
export Re=1.0e+2
export Sc=1.0e+1

# give name of the directory in which the initial conditions
#   (incl. domain size etc.) are stored as an argument
dirname_ic=initial_condition/output

mpirun -n 2 --oversubscribe ./a.out ${dirname_ic}
