#!/bin/bash

PROG="../src/grins"

INPUT="./input_files/axi_magnetostatics.in"

#PETSC_OPTIONS="-pc_type ilu"
PETSC_OPTIONS="ksp_type preonly -pc_type lu -pc_factor_mat_solver_package superlu"

$PROG $INPUT $PETSC_OPTIONS 
