#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/axi_poiseuille_flow_input.in"
TESTDATA="./axi_ns_poiseuille_flow.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app input=$INPUT vars='u_r u_z p' norms='L2' tol='5.0e-8' u_r_L2_error='1.0e-10' u_z_L2_error='1.0e-10' p_L2_error='5.13691e-09' u_r_exact_soln='0.0' u_z_exact_soln='100.0/4.0*(1-x^2)' p_exact_soln='100*(1-y)' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
