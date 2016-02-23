#!/bin/bash

PROG="${GRINS_TEST_DIR}/grins_flow_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/penalty_poiseuille_stab.in ${GRINS_TEST_DATA_DIR}/penalty_poiseuille_stab.xdr 1e-9"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
