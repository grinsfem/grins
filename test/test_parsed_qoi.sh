#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/parsed_qoi.in"

${LIBMESH_RUN:-} $PROG $INPUT

