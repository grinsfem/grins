#!/bin/bash

PROG="${GRINS_TEST_DIR}/grins_flow_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/locally_refine.in ${GRINS_TEST_DATA_DIR}/locally_refine.xdr 1.0e-8"

# Make sure libMesh is still respecting ksp_rtol by using it to
# override a lousy input file tolerance
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 4 -ksp_rtol 1.0e-8"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
