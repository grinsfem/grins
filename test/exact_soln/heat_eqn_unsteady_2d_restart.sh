#!/bin/bash

set -e

INPUT_1="${GRINS_TEST_INPUT_DIR}/heat_eqn_unsteady_2d_restart_pt1.in"
INPUT_2="${GRINS_TEST_INPUT_DIR}/heat_eqn_unsteady_2d_restart_pt2.in"

TESTDATA_NOTUSED="./heat_eqn_unsteady_2d_restart_pt1.xdr ./heat_eqn_unsteady_2d_restart_pt2.xdr ./heat_eqn_unsteady_2d_restart_pt1_mesh.xda ./heat_eqn_unsteady_2d_restart_pt1.24.xdr ./heat_eqn_unsteady_2d_restart_pt1.24_mesh.xda"
TESTDATA="./heat_eqn_unsteady_2d_restart_pt2.24.xdr"

# First run the case to generate the file to restart from
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT_1

# Next run the case restarting from the file dumped out in part 1
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT_2

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
                 input=$INPUT_2 \
                 vars='u' \
                 norms='L2' \
                 tol='1.0e-7' \
                 u_L2_error='0.00205805' \
                 u_exact_soln='tf:=50*0.02;sin(pi*x)*sin(pi*y)*sin(pi*t)' \
                 test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA $TESTDATA_NOTUSED
