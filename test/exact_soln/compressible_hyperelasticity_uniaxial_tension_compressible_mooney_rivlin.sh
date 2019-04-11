#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/compressible_hyperelasticity_uniaxial_tension_compressible_mooney_rivlin.in"
TESTDATA="./compressible_hyperelasticity_uniaxial_tension_compressible_mooney_rivlin.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app --input $INPUT \
 L0='1.0' C1='3.1' C2='4.159' lambda2='0.8' \
 lamb1exp='(1+sqrt(1+4*c1/c2*(1-l2*l2)))/(2*l2^2)' \
 uext='L0:=${L0};c1:=${C1};c2:=${C2};l2:=${lambda2};l1:=${lamb1exp};L0*(l1-1)' \
 ulat='L0:=${L0};l2:=${lambda2};L0*(1-l2)' \
 Ux_exact_soln='${ulat}*(-1)*x' \
 Uy_exact_soln='${ulat}*(-1)*y' \
 Uz_exact_soln='${uext}*z' \
 vars='Ux Uy Uz' norms='L2' tol='1.0e-10' Ux_L2_error='1.0e-10' Uy_L2_error='1.0e-10' Uz_L2_error='1.0e-10' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
