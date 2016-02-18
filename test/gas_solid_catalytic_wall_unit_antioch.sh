#!/bin/bash

PROG="${GRINS_TEST_DIR}/gas_solid_catalytic_wall_unit"

INPUT="antioch ${GRINS_TEST_DIR}/input_files/gas_surface.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
