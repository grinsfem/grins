#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/poiseuille_flow_input.in"
TESTDATA="./ns_poiseuille_flow.xda"

# A more partitioning-robust preconditioner
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_type ilu -sub_pc_factor_mat_ordering_type 1wd -sub_pc_factor_levels 4 -sub_pc_factor_shift_type nonzero"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app --input $INPUT vars='u v p' norms='L2' tol='1.0e-10' u_L2_error='1.0e-10' v_L2_error='1.0e-10' p_L2_error='1.0e-10' u_exact_soln='4*y*(1-y)' v_exact_soln='0.0' p_exact_soln='120.0+(80.0-120.0)/5.0*x' test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
