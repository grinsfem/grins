#!/bin/bash

if [ "${GRINS_LIBMESH_DIM}" -gt 2 ]                                                                                                                                                                                                                 
then
   PROG="${GRINS_TEST_DIR}/grins_flow_regression"

   INPUT="${GRINS_TEST_INPUT_DIR}/redistribute.in ${GRINS_TEST_DATA_DIR}/redistributed.xdr 1.0e-8"

   #PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps"
   PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

   ${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
else
   # If LIBMESH_DIM !=3, we skip this test
   exit 77
fi
