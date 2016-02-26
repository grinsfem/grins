#!/bin/bash

if [ "${GRINS_LIBMESH_DIM}" -gt 2 ]
then
   PROG="${GRINS_TEST_DIR}/test_thermally_driven_flow"

   INPUT="${GRINS_TEST_INPUT_DIR}/thermally_driven_3d_flow.in ${GRINS_TEST_DATA_DIR}/thermally_driven_3d.xdr"

   LIBMESH_OPTIONS="--n_threads=6"

   ${LIBMESH_RUN:-} $PROG $INPUT
else
   # If LIBMESH_DIM != 3, then we skip this test
   exit 77
fi
