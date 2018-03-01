#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/axi_thermally_driven_flow.in"
DATA="${GRINS_TEST_DATA_DIR}/axi_thermally_driven.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_mat_ordering_type 1wd -sub_pc_factor_shift_type nonzero"
${LIBMESH_RUN:-} $PROG --input $INPUT soln-data=$DATA vars='u v p T' norms='L2 H1' tol='2.0e-8' $PETSC_OPTIONS
