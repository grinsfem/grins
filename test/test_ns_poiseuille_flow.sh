#!/bin/bash

PROG="./test_ns_poiseuille_flow"

INPUT="./input_files/poiseuille_flow_input.in"

PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps"

$PROG $INPUT 
