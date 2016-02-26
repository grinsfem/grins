#!/bin/bash

PROG="${GRINS_TEST_DIR}/suspended_cable_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/suspended_cable_test.in ${GRINS_TEST_DATA_DIR}/suspended_cable_test.xdr"

${LIBMESH_RUN:-} $PROG $INPUT
