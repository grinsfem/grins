#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/dirichlet_fem.in ${GRINS_TEST_DATA_DIR}/dirichlet_fem.xdr"

${LIBMESH_RUN:-} $PROG $INPUT
