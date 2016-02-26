#!/bin/bash

if [ "${GRINS_LIBMESH_DIM}" -gt 1 ]; then
   PROG="${GRINS_TEST_DIR}/grins_flow_regression"

   INPUT="${GRINS_TEST_INPUT_DIR}/2d_pseudofan.in ${GRINS_TEST_DATA_DIR}/2d_pseudofan.xdr 1.0e-8"

   PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

   $PROG $INPUT $PETSC_OPTIONS 
else
   # If LIBMESH_DIM < 2, we skip this test
   exit 77
fi
