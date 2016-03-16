#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/parsed_qoi.in"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_mat_ordering_type 1wd -sub_pc_factor_shift_type nonzero"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS

