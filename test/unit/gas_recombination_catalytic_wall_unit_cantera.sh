#!/bin/bash

PROG="${GRINS_TEST_DIR}/gas_recombination_catalytic_wall_unit"

INPUT="cantera ${GRINS_TEST_INPUT_DIR}/cantera_chem_thermo.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
