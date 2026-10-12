#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/L2_2d.in"
TESTDATA="./L2_2d.xda"

# First run the case with grins.
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins \
  $INPUT \
  $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
  --input $INPUT \
  vars='u' norms='L2' tol='1.0e-10' \
  u_L2_error='0' \
  u_exact_soln="x*(1-x)*y*(1-y)" \
  test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
