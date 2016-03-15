#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/convection_diffusion_unsteady_2d.in"
TESTDATA_NOTUSED="./convection_diffusion_unsteady_2d.xdr"
TESTDATA="./convection_diffusion_unsteady_2d.49.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
                 input=$INPUT \
                 vars='u' \
                 norms='L2' \
                 tol='1.0e-10' \
                 u_L2_error='6.5547183600866560e-03' \
                 u_exact_soln='tf:=50*0.025;exp(-((x-0.8*tf-0.2)^2+(y-0.8*tf-0.2)^2)/(0.01*(4.0*tf+1.0)))/(4.0*tf+1.0)' \
                 test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA $TESTDATA_NOTUSED
