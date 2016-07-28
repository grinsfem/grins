#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_cable_oned_newmark_rayleigh_damping.in"
TESTDATA_NOTUSED="./elastic_cable_oned_newmark_rayleigh_damping.xdr"
TESTDATA="./elastic_cable_oned_newmark_rayleigh_damping.99.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app input=$INPUT vars='Ux' norms='L2' tol='1.0e-9' Ux_L2_error='0.000152487' Ux_exact_soln='0.1*x/8' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA $TESTDATA_NOTUSED
