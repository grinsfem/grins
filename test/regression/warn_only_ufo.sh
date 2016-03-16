#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/ufo_unit.in"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

GRINS_OPTIONS="--warn-only-unused-var"

${LIBMESH_RUN:-} $PROG $INPUT $GRINS_OPTIONS $PETSC_OPTIONS
