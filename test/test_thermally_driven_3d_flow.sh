#!/bin/bash

PROG="./test_thermally_driven_flow"

INPUT="./input_files/thermally_driven_3d_flow.in"

PETSC_OPTIONS="-pc_type ilu" 

#-pc_factor_mat_solver_package mumps"

$PROG $INPUT 
