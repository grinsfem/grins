#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/L2_iga.in"
TESTDATA="./L2_iga.xda"

# First run the case with grins.
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins \
  $INPUT \
  $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
  --input $INPUT \
  vars='u' norms='L2' tol='1.0e-10' \
  u_L2_error='0' \
  u_exact_soln="8+4*x-2*y+z" \
  test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
