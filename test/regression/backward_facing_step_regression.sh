#!/bin/bash

PROG="${GRINS_TEST_DIR}/grins_flow_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/backward_facing_step.in ${GRINS_TEST_DATA_DIR}/backward_facing_step.xdr 1.0e-8"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
