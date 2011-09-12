#!/bin/bash

PROG="./test_ns_couette_flow_2d_x"

INPUT="./input_files/couette_flow_input_2d_x.in"

PETSC_OPTIONS="-pc_type lu -pc_factor_mat_solver_package mumps"
$PROG $INPUT
