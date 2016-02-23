#!/bin/bash

PROG="${GRINS_TEST_DIR}/ns_couette_flow_2d_y"

INPUT="${GRINS_TEST_INPUT_DIR}/couette_flow_input_2d_y.in"

${LIBMESH_RUN:-} $PROG $INPUT
