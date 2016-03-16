#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/dirichlet_nan.in ${GRINS_TEST_DATA_DIR}/dirichlet_nan.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
