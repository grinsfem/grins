#!/bin/bash

PROG="${GRINS_TEST_DIR}/elastic_sheet_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_mooney_rivlin_circle_hookean_stiffeners_regression.in ${GRINS_TEST_DATA_DIR}/elastic_mooney_rivlin_circle_hookean_stiffeners_regression.xdr"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
