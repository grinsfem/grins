#!/bin/bash

PROG="${GRINS_TEST_DIR}/stokes_poiseuille_flow_parsed_viscosity"

INPUT="${GRINS_TEST_INPUT_DIR}/stokes_poiseuille_flow_parsed_viscosity_input.in"

${LIBMESH_RUN:-} $PROG $INPUT 
