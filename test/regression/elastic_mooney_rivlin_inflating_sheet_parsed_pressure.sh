#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_mooney_rivlin_inflating_sheet_parsed_pressure_regression.in"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

# Solution output from GRINS run
SOLNDATA="./elastic_mooney_rivlin_inflating_sheet_parsed_pressure_regression.xda"

# Gold data used for regression comparsion.
# This should be the same as the existing data for the
# elastic_mooney_rivlin_inflating_sheet_regression so we use it
# instead of create a new one.
GOLDDATA="${GRINS_TEST_DATA_DIR}/elastic_mooney_rivlin_inflating_sheet_regression.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${GRINS_TEST_DIR}/regression_testing_app \
                 input=$INPUT \
                 system_name='StretchedElasticSheet' \
                 vars='u v w' \
                 norms='L2 H1' \
                 tol='5.0e-8' \
                 gold-data=$GOLDDATA \
                 soln-data=$SOLNDATA

# Now remove the test turd
rm $SOLNDATA
