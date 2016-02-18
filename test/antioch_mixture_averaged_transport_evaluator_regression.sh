#!/bin/bash

PROG="${GRINS_TEST_DIR}/antioch_mixture_averaged_transport_evaluator_regression"

INPUT="${GRINS_TEST_DIR}/input_files/antioch.in"

${LIBMESH_RUN:-} $PROG $INPUT
