#!/bin/bash

PROG="${GRINS_TEST_DIR}/gas_solid_catalytic_wall"

INPUT="antioch ${GRINS_TEST_INPUT_DIR}/gas_surface.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
