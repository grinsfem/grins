#!/bin/bash

set -e

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/hot_cylinder.in"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_type lu -sub_pc_factor_shift_type nonzero"

# Solution output from GRINS run
SOLNDATA="./hot_cylinder.xda"

# Gold data used for regression comparsion
GOLDDATA="${GRINS_TEST_DATA_DIR}/hot_cylinder_gold.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${GRINS_TEST_DIR}/regression_testing_app \
   input=$INPUT \
   vars='Ux Uy p T' \
   norms='L2 H1' \
   tol='2.0e-7' \
   gold-data=$GOLDDATA \
   soln-data=$SOLNDATA

# Now remove the test turd
rm $SOLNDATA
