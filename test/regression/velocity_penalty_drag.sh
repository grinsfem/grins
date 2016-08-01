#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/velocity_penalty_drag.in"
DATA="${GRINS_TEST_DATA_DIR}/velocity_penalty_drag.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_type lu -sub_pc_factor_shift_type nonzero"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='u v p' norms='L2 H1' tol='1.0e-9' $PETSC_OPTIONS
