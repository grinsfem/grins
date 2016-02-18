#!/bin/sh

PROG="${GRINS_TEST_DIR}/3d_low_mach_jacobians_xy"

INPUT="${GRINS_TEST_INPUT_DIR}/3d_low_mach_jacobians_xy.in"

${LIBMESH_RUN:-} $PROG $INPUT

