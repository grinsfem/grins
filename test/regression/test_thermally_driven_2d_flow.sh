#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/thermally_driven_2d_flow.in ${GRINS_TEST_DATA_DIR}/thermally_driven_2d.xdr"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
