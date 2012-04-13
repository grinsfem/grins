#!/bin/bash

PROG="./test_axi_ns_con_cyl_flow"

INPUT="./input_files/axi_con_cyl_flow.in"

PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps"
#PETSC_OPTIONS="-pc_type ilu -pc_factor_levels 20"
#PETSC_OPTIONS="-pc_type lu -pc_factor_shift_type NONZERO"

$PROG $INPUT $PETSC_OPTIONS 
