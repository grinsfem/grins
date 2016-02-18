#!/bin/bash

PROG="${GRINS_TEST_DIR}/antioch_mixture_unit"

INPUT="${GRINS_TEST_DIR}/input_files/antioch.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
