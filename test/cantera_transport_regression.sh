#!/bin/bash

PROG="${GRINS_TEST_DIR}/cantera_transport_regression"

INPUT="${GRINS_TEST_DIR}/input_files/cantera_transport.in"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
