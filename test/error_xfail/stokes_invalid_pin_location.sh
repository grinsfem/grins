#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/stokes_invalid_pin_location_unit.in"

# THIS IS DUMMY, the test is expected to fail before getting to the data comparison
DATA="${GRINS_TEST_DATA_DIR}/laplace_parsed_source_regression.xdr"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='T' norms='L2 H1' tol='1.0e-10'
