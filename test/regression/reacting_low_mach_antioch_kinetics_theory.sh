#!/bin/bash

PROG="${GRINS_TEST_DIR}/reacting_low_mach_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/reacting_low_mach_antioch_kinetics_theory_regression.in ${GRINS_TEST_DATA_DIR}/reacting_low_mach_antioch_kinetics_theory_regression.xdr"

# A MOAB preconditioner
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type ilu -sub_pc_factor_levels 10 -sub_pc_factor_shift nonzero"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
