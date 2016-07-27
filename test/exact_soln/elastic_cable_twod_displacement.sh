#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_cable_twod_displacement.in"
TESTDATA="./elastic_cable_twod_displacement.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app input=$INPUT vars='Ux Uy' norms='L2' tol='1.0e-10' Ux_L2_error='1.0e-10' Uy_L2_error='1.0e-10' Ux_exact_soln='0.1*x/8' Uy_exact_soln='0.2*x/8' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
