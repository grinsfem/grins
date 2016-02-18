#!/bin/sh

PROG="${GRINS_TEST_DIR}/3d_low_mach_jacobians_yz"

INPUT="${GRINS_TEST_INPUT_DIR}/3d_low_mach_jacobians_yz.in"

${LIBMESH_RUN:-} $PROG $INPUT

