#!/bin/bash

PROG="${GRINS_TEST_DIR}/axisym_reacting_low_mach_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/axisym_reacting_low_mach_antioch_cea_constant_regression.in"
DATA="${GRINS_TEST_DATA_DIR}/axisym_reacting_low_mach_antioch_cea_constant_regression.xdr"

# This problem does *not* like much additive Schwarz overlap...
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_type ilu -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 10"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='u v T p w_N2 w_N' norms='L2 H1' tol='2.0e-8' $PETSC_OPTIONS
