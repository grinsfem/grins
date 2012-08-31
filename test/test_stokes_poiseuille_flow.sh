#!/bin/bash

PROG="./test_stokes_poiseuille_flow"

INPUT="./input_files/stokes_poiseuille_flow_input.in"

#PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps"
PETSC_OPTIONS="-ksp_type gmres -pc_type ilu -pc_factor_levels 4"

$PROG $INPUT 
