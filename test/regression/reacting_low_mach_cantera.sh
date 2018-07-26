#!/bin/bash

set -e

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/reacting_low_mach_cantera_regression.in"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type ilu -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 10"

# Solution output from GRINS run
SOLNDATA="./reacting_low_mach_cantera_regression.xda"

# Gold data used for regression comparsion
GOLDDATA="${GRINS_TEST_DATA_DIR}/reacting_low_mach_cantera_regression.xda.gz"

if [ $GRINS_CANTERA_ENABLED == 1 ]; then
   # First run the case with grins
   ${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

   # Now run the test part to make sure we're getting the correct thing
   ${GRINS_TEST_DIR}/regression_testing_app \
      input=$INPUT \
      vars='u v p T w_N w_N2' \
      norms='L2 H1' \
      tol='2.0e-7' \
      gold-data=$GOLDDATA \
      soln-data=$SOLNDATA

   # Now remove the test turd
   rm $SOLNDATA
else
   exit 77;
fi
