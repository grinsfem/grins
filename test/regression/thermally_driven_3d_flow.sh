#!/bin/bash

if [ "${GRINS_LIBMESH_DIM}" -gt 2 ]
then
   PROG="${GRINS_TEST_DIR}/generic_solution_regression"

   INPUT="${GRINS_TEST_INPUT_DIR}/thermally_driven_3d_flow.in"
   DATA="${GRINS_TEST_DATA_DIR}/thermally_driven_3d.xdr"

   ${LIBMESH_RUN:-} $PROG --input $INPUT soln-data=$DATA vars='u v p T' norms='L2 H1' tol='8.0e-9'
else
   # If LIBMESH_DIM != 3, then we skip this test
   exit 77
fi
