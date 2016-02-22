#!/bin/bash

PROG="${GRINS_TEST_DIR}/cantera_mixture_unit"

INPUT="${GRINS_TEST_INPUT_DIR}/cantera_chem_thermo.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
