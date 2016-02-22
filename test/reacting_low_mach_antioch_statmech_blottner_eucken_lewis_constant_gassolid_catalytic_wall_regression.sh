#!/bin/bash

PROG="${GRINS_TEST_DIR}/reacting_low_mach_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.in ${GRINS_TEST_DATA_DIR}/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.xdr"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
