#!/bin/bash

PROG="./test_ns_couette_flow"

INPUT="./input_files/couette_flow_input.in"

PETSC_OPTIONS="-pc_type lu -pc_factor_mat_solver_package mumps"
$PROG $INPUT
