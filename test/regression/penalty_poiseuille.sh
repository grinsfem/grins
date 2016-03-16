#!/bin/bash

PROG="${GRINS_TEST_DIR}/grins_flow_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/penalty_poiseuille.in ${GRINS_TEST_DATA_DIR}/penalty_poiseuille.xdr 1e-9"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 2 -sub_pc_type ilu -sub_pc_factor_shift_type nonzero -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
