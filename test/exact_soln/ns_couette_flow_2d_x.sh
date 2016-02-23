#!/bin/bash

PROG="${GRINS_TEST_DIR}/ns_couette_flow_2d_x"

INPUT="${GRINS_TEST_INPUT_DIR}/couette_flow_input_2d_x.in"

${LIBMESH_RUN:-} $PROG $INPUT
