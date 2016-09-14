#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/suspended_cable_test.in"
DATA="${GRINS_TEST_DATA_DIR}/suspended_cable_test.xdr"

PETSC_OPTIONS="-ksp_type gmres -pc_type bjacobi -sub_pc_type lu -sub_pc_factor_shift_type nonzero -ksp_converged_reason"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='Ux Uy' norms='L2 H1' tol='1.0e-10' $PETSC_OPTIONS
