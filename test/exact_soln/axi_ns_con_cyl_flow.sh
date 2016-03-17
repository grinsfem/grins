#!/bin/bash

set -e

INPUT="${GRINS_TEST_INPUT_DIR}/axi_con_cyl_flow.in"
TESTDATA="./axi_ns_con_cyl_flow.xdr"

# A MOAB preconditioner
PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 12 -sub_pc_type ilu -sub_pc_factor_mat_ordering_type 1wd -sub_pc_factor_levels 4 -sub_pc_factor_shift_type nonzero"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/generic_exact_solution_testing_app \
                 input=$INPUT \
                 vars='u_r u_z' \
                 norms='L2' \
                 tol='1.0e-10' \
                 u_r_L2_error='1.0e-10' \
                 u_z_L2_error='1.0e-10' \
                 u_r_exact_soln='0.0' \
                 u_z_exact_soln='2*log(2/x)/log(2)' \
                 test_data=$TESTDATA

# Now remove the test turd
rm $TESTDATA
