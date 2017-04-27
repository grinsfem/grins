#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/poisson_periodic_3d_yz.in"
TESTDATA="./poisson_periodic_3d_yz.xda"

PETSC_OPTIONS="-ksp_type cg -pc_type bjacobi -sub_pc_type icc"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app --input $INPUT vars='u' norms='L2' tol='1.0e-8' u_L2_error='0.0106696' u_exact_soln='sin(pi*x)*cos(2*pi*z)*cos(2*pi*y)' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
