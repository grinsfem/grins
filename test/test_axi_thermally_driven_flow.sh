#!/bin/bash

PROG="./test_thermally_driven_flow"

INPUT="./input_files/axi_thermally_driven_flow.in"

PETSC_OPTIONS="-pc_type ilu"
#PETSC_OPTIONS="-pc_type lu -pc_factor_mat_solver_package mumps"

$PROG $INPUT $PETSC_OPTIONS 
