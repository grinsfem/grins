#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/ozone_flame_antioch_constant_regression.in"
DATA="${GRINS_TEST_DATA_DIR}/ozone_flame_antioch_constant_regression.xdr"

PETSC_OPTIONS="-ksp_type gmres -pc_type bjacobi -sub_pc_type lu -sub_pc_factor_shift_type nonzero"

if [ $GRINS_ANTIOCH_ENABLED == 1 ]; then
   ${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='Ux Uy T p Y_O Y_O2 Y_O3' norms='L2 H1' tol='1.0e-6' $PETSC_OPTIONS
else
   exit 77;
fi
