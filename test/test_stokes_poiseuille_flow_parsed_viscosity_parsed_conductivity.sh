#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_stokes_poiseuille_flow_parsed_viscosity_parsed_conductivity"

INPUT="${GRINS_TEST_INPUT_DIR}/stokes_poiseuille_flow_parsed_viscosity_parsed_conductivity_input.in"

${LIBMESH_RUN:-} $PROG $INPUT 
