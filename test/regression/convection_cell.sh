#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/convection_cell_regression.in ${GRINS_TEST_DATA_DIR}/convection_cell_regression.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
