#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/poisson_periodic_2d_x.in"
TESTDATA="./poisson_periodic_2d_x.xda"

PETSC_OPTIONS="-ksp_type cg -pc_type bjacobi -sub_pc_type icc"

# First run the case with grins.  We'll override part of the input
# file on the command line, too, to test that feature.
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins \
  $INPUT \
  Mesh/Generation/n_elems_x=20 \
  Mesh/Generation/n_elems_y=20 \
  $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
  input=$INPUT \
  Mesh/Generation/n_elems_x=20 \
  Mesh/Generation/n_elems_y=20 \
  vars='u' norms='L2' tol='1.0e-8' \
  u_L2_error='0.00349104' u_exact_soln='sin(pi*y)*cos(2*pi*x)' \
  test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
