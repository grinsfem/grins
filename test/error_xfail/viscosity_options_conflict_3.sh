#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/viscosity_options_conflict_3.in"

${LIBMESH_RUN:-} $PROG $INPUT
