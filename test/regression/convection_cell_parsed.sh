#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/convection_cell_parsed_regression.in"
DATA="${GRINS_TEST_DATA_DIR}/convection_cell_regression.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $PETSC_OPTIONS input=$INPUT soln-data=$DATA vars='u v p T' norms='L2 H1' tol='8.0e-9'
