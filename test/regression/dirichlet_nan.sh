#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/dirichlet_nan.in"
DATA="${GRINS_TEST_DATA_DIR}/dirichlet_nan.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $PETSC_OPTIONS input=$INPUT soln-data=$DATA vars='u v p T' norms='L2 H1' tol='8.0e-9'
