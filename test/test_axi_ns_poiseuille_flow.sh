#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_axi_ns_poiseuille_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/axi_poiseuille_flow_input.in"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
