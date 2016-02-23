#!/bin/bash

PROG="${GRINS_TEST_DIR}/ufo_unit"

INPUT="${GRINS_TEST_INPUT_DIR}/ufo_unit.in"

${LIBMESH_RUN:-} $PROG $INPUT
