#!/bin/sh

PROG="${GRINS_TEST_DIR}/3d_low_mach_jacobians_xz"

INPUT="${GRINS_TEST_INPUT_DIR}/3d_low_mach_jacobians_xz.in"

${LIBMESH_RUN:-} $PROG $INPUT

