#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_stokes_poiseuille_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/stokes_poiseuille_flow_input.in"

${LIBMESH_RUN:-} $PROG $INPUT 
