#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/elastic_mooney_rivlin_membrane_cantilever_unsteady_newmark_regression.in"
TESTDATA="./elastic_mooney_rivlin_membrane_cantilever_unsteady_newmark_regression.99.xda.gz"
TESTDATA_NOTUSED="./elastic_mooney_rivlin_membrane_cantilever_unsteady_newmark_regression.xda.gz"
GOLDDATA="${GRINS_TEST_DATA_DIR}/elastic_mooney_rivlin_membrane_cantilever_unsteady_newmark_regression.99.xda.gz"

PETSC_OPTIONS="-ksp_type cg -pc_type bjacobi -sub_pc_type icc"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/regression_testing_app input=$INPUT vars='Ux Uy' norms='L2' tol='1.0e-8' gold-data=$GOLDDATA soln-data=$TESTDATA

# Now remove the test turd
rm $TESTDATA $TESTDATA_NOTUSED
