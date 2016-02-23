#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/laplace_parsed_source_regression.in"
DATA="${GRINS_TEST_DATA_DIR}/laplace_parsed_source_regression.xdr"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='T' norms='L2 H1' tol='1.0e-10'
