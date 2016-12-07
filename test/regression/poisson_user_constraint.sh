#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/poisson_user_constraint.in"
TESTDATA="./poisson_user_constraint.xda.gz"
GOLDDATA="${GRINS_TEST_DATA_DIR}/poisson_user_constraint.xda.gz"

PETSC_OPTIONS="-ksp_type cg -pc_type bjacobi -sub_pc_type icc"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/regression_testing_app input=$INPUT vars='u' norms='L2' tol='1.0e-8' gold-data=$GOLDDATA soln-data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
