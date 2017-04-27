#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_mooney_rivlin_square_hookean_stiffeners_rayleigh_damping_newmark_regression.in"
DATA="${GRINS_TEST_DATA_DIR}/elastic_mooney_rivlin_square_hookean_stiffeners_rayleigh_damping_newmark_regression.xdr"

PETSC_OPTIONS="-ksp_type cg -pc_type bjacobi -sub_pc_type icc"

${LIBMESH_RUN:-} $PROG $PETSC_OPTIONS --input $INPUT soln-data=$DATA vars='u v w' norms='L2 H1' tol='1.0e-10'
