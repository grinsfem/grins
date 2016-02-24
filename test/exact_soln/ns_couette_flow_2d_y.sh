#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/couette_flow_input_2d_y.in"
TESTDATA="./ns_couette_flow_2d_y.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${GRINS_TEST_DIR}/generic_exact_solution_testing_app input=$INPUT vars='u v' norms='L2' tol='1.0e-10' u_L2_error='1.0e-10' v_L2_error='1.0e-10' u_exact_soln='0.0' v_exact_soln='10.0*x' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
