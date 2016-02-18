#!/bin/bash

PROG="${GRINS_TEST_DIR}/vorticity_qoi"

INPUT="${GRINS_TEST_DIR}/input_files/vorticity_qoi.in"

${LIBMESH_RUN:-} $PROG $INPUT

