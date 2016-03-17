#!/bin/bash

PROG="${GRINS_TEST_DIR}/reacting_low_mach_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/reacting_low_mach_antioch_cea_constant_prandtl_regression.in ${GRINS_TEST_DATA_DIR}/reacting_low_mach_antioch_cea_constant_prandtl_regression.xdr"

# A MOAB preconditioner
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type ilu -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 10"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
