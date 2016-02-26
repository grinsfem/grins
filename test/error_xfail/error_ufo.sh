#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/ufo_unit.in"

${LIBMESH_RUN:-} $PROG $INPUT
