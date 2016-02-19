#!/bin/bash

PROG="${GRINS_TEST_DIR}/elastic_sheet_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_mooney_rivlin_square_hookean_stiffeners_unifref_regression.in ${GRINS_TEST_DATA_DIR}/elastic_mooney_rivlin_square_hookean_stiffeners_unifref_regression.xdr"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

# FIXME: We currently get a libMesh error when trying to run on
# ParallelMesh on more than 4 processors
#${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
$PROG $INPUT $PETSC_OPTIONS
