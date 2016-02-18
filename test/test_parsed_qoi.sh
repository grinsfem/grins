#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_DIR}/input_files/parsed_qoi.in"

${LIBMESH_RUN:-} $PROG $INPUT

