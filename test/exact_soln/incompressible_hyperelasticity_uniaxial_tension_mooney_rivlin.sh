#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/incompressible_hyperelasticity_uniaxial_tension_mooney_rivlin.in"
TESTDATA="./incompressible_hyperelasticity_uniaxial_tension_mooney_rivlin.xda.gz"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app --input $INPUT \
 L0='1.0' uz='0.2' lam1='(uz+L0)/L0' lam2='1/(sqrt(${lam1})' \
 ulat='L0:=${L0};uz:=${uz};(sqrt(L0^3/(uz+L0))-L0)' \
 Ux_exact_soln='${ulat}*x' \
 Uy_exact_soln='${ulat}*y' \
 Uz_exact_soln='${uz}*z' \
 vars='Ux Uy Uz' norms='L2' tol='1.0e-10' Ux_L2_error='1.0e-10' Uy_L2_error='1.0e-10' Uz_L2_error='1.0e-10' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
