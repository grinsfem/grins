#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/ins_invalid_pin_location_unit.in"

${LIBMESH_RUN:-} $PROG $INPUT
