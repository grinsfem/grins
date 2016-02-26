#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

INPUT="${GRINS_TEST_INPUT_DIR}/axi_thermally_driven_flow.in ${GRINS_TEST_DATA_DIR}/axi_thermally_driven.xdr"

${LIBMESH_RUN:-} $PROG $INPUT
