#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/plane_strain_incompressible_hyperelasticity_uniaxial_tension_mooney_rivlin.in"
TESTDATA="./plane_strain_incompressible_hyperelasticity_uniaxial_tension_mooney_rivlin.xda.gz"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app --input $INPUT \
 L0='1.0' ux='0.2' \
 ulat='L0:=${L0};ux:=${ux};((L0^2/(ux+L0))-L0)' \
 Ux_exact_soln='${ux}*x' \
 Uy_exact_soln='${ulat}*y' \
 vars='Ux Uy' norms='L2' tol='1.0e-10' Ux_L2_error='1.0e-10' Uy_L2_error='1.0e-10' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
