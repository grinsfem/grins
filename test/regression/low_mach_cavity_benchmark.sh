#!/bin/bash

PROG="${GRINS_TEST_DIR}/low_mach_cavity_benchmark_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/low_mach_cavity_benchmark_regression_input.in"

# A MOAB preconditioner
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 8 -sub_pc_type ilu -sub_pc_factor_shift_type nonzero -sub_pc_factor_mat_ordering_type 1wd -sub_pc_factor_levels 10 -ksp_rtol 1.0e-8"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
