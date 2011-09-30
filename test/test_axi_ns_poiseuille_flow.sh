#!/bin/bash

PROG="./test_axi_ns_poiseuille_flow"

INPUT="./input_files/axi_poiseuille_flow_input.in"

PETSC_OPTIONS="-pc_type lu -pc_factor_mat_solver_package mumps"

$PROG $INPUT $PETSC_OPTIONS 
