#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/ufo_unit.in"

GRINS_OPTIONS="--warn-only-unused-var"

${LIBMESH_RUN:-} $PROG $INPUT $GRINS_OPTIONS
