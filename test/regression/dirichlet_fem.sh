#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/dirichlet_fem.in"
DATA="${GRINS_TEST_DATA_DIR}/dirichlet_fem.xdr"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='u v p T' norms='L2 H1' tol='8.0e-9'
