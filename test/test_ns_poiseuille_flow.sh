#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_ns_poiseuille_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/poiseuille_flow_input.in"

${LIBMESH_RUN:-} $PROG $INPUT 
