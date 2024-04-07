#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/convection_diffusion_steady_1d.in"
TESTDATA="./convection_diffusion_steady_1d.xda.gz"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
                 --input $INPUT \
                 vars='u' \
                 norms='L2' \
                 tol='1.0e-10' \
                 u_L2_error='4.862013076677351e-04' \
                 u_exact_soln='x-(1-exp(40*x))/(1-exp(40))' \
                 test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
