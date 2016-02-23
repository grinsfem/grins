#!/bin/bash

PROG="${GRINS_TEST_DIR}/test_turbulent_channel"

INPUT="${GRINS_TEST_INPUT_DIR}/sa_2d_turbulent_channel_regression.in"

MESH_1D="${GRINS_TEST_DATA_DIR}/turbulent_channel_Re944_grid.xda"
DATA_1D="${GRINS_TEST_DATA_DIR}/turbulent_channel_soln.xda"

DATA="${GRINS_TEST_DATA_DIR}/sa_2d_turbulent_channel_regression.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 8 -sub_pc_factor_mat_ordering_type 1wd -sub_pc_type ilu -sub_pc_factor_levels 6"

${LIBMESH_RUN:-} $PROG $INPUT soln-data=$DATA vars='u v p nu' norms='L2 H1' tol='2.0e-8' mesh-1d=$MESH_1D data-1d=$DATA_1D $PETSC_OPTIONS
