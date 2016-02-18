#!/bin/bash

PROG="${GRINS_BUILDSRC_DIR}/grins"

INPUT="${GRINS_TEST_INPUT_DIR}/simple_ode.in"

# FIXME: In theory we should be able to solve a scalar problem on
# multiple processors, where ranks 1+ just twiddle their thumbs.
# In practice we get libMesh errors.
#${LIBMESH_RUN:-} $PROG $INPUT
$PROG $INPUT
