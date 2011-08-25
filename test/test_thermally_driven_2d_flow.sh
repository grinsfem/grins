#!/bin/bash

PROG="./test_thermally_driven_flow"

INPUT="./input_files/thermally_driven_2d_flow.in"

PETSC_OPTIONS="-pc_type lu -pc_factor_mat_solver_package mumps"

$PROG $INPUT 
