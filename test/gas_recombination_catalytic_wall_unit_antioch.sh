#!/bin/bash

PROG="${GRINS_TEST_DIR}/gas_recombination_catalytic_wall_unit"

INPUT="antioch ${GRINS_TEST_INPUT_DIR}/antioch.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
