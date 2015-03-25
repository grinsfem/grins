#!/bin/bash

PROG="../test/test_turbulent_channel"

INPUT="../test/input_files/turbulent_channel_input.in"

#PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps"
#PETSC_OPTIONS="-ksp_type gmres -pc_type ilu -pc_factor_levels 4"
PETSC_OPTIONS="-ksp_type preonly -pc_type lu"

$PROG $INPUT $PETSC_OPTIONS
