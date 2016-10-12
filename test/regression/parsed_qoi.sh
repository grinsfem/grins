#!/bin/bash

set -e

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/parsed_qoi.in"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_mat_ordering_type 1wd -sub_pc_factor_shift_type nonzero"

# Solution output from GRINS run
SOLNDATA="./parsed_qoi.xda"

# QoI output from GRINS run
QOIDATA="./parsed_qoi.qoi.dat"

# Gold data used for regression
GOLDDATA="${GRINS_TEST_DATA_DIR}/parsed_qoi_gold.xdr"

# First run the case with grins
${LIBMESH_RUN:-} ${GRINS_BUILDSRC_DIR}/grins $INPUT $PETSC_OPTIONS

# Now run the test part to make sure we're getting the correct thing
${LIBMESH_RUN:-} ${GRINS_TEST_DIR}/regression_testing_app \
                 input=$INPUT \
                 vars='u v' \
                 norms='L2 H1' \
                 tol='1.0e-10' \
                 gold-data=$GOLDDATA \
                 soln-data=$SOLNDATA \
                 qoi-data=$QOIDATA \
                 gold-qoi-names='parsed_interior parsed_boundary weighted_flux' \
                 gold-qoi-values='8.3333333341759397e-01 1.1666666666666670e+00 -4.8958333341067979e+00'

# Now remove the test turd
rm $SOLNDATA $QOIDATA
