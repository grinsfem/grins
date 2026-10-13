#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/poisson_1d_penalty.in"
TESTDATA="./poisson_1d_penalty.xda"

# First run the case with grins.
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins \
  $INPUT \
  $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
export eps=1e-2
export offset=2.5
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
  --input $INPUT \
  vars='u' norms='L2' tol='1.0e-10' \
  u_L2_error='1e-10' \
  u_exact_soln="x*((1+2*$eps)/(1+$eps)-x)+$offset" \
  test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
