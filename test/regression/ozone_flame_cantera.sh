#!/bin/bash

set -e

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/ozone_flame_cantera_regression.in"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type lu -sub_pc_factor_shift_type nonzero"

# Solution output from GRINS run
SOLNDATA="./ozone_flame_cantera_regression.xda"

# Gold data used for regression comparsion
GOLDDATA="${GRINS_TEST_DATA_DIR}/ozone_flame_cantera_regression.xdr"

if [ $GRINS_CANTERA_ENABLED == 1 ]; then
   # First run the case with grins
   ${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

   # Now run the test part to make sure we're getting the correct thing
   ${GRINS_TEST_DIR}/regression_testing_app \
      input=$INPUT \
      system_name='Ozone' \
      vars='Ux Uy p T Y_O Y_O2 Y_O3' \
      norms='L2 H1' \
      tol='3.0e-6' \
      gold-data=$GOLDDATA \
      soln-data=$SOLNDATA

   # Now remove the test turd
   rm $SOLNDATA
else
   exit 77;
fi
