#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/thermally_driven_2d_flow.in ${GRINS_TEST_DATA_DIR}/thermally_driven_2d.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type ilu -sub_pc_factor_levels 10"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
