#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/suspended_cable_test.in"
DATA="${GRINS_TEST_DATA_DIR}/suspended_cable_test.xdr"

if grep "PETSC_HAVE_MUMPS 1" "$PETSC_DIR/$PETSC_ARCH"/include/*; then
  PETSC_OPTIONS="-pc_type lu -pc_factor_mat_solver_package mumps"
else
  PETSC_OPTIONS="-ksp_type cg -pc_type bjacobi -sub_pc_type icc -sub_pc_factor_shift_type nonzero -ksp_converged_reason"
fi

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='Ux Uy' norms='L2 H1' tol='1.0e-10' $PETSC_OPTIONS
