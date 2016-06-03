#!/bin/bash

PROG="${GRINS_TEST_DIR}/generic_solution_regression"

# Note the tolerance is high because of the H1 error in the pressure.
# If we had a per-variable tolerance the others could be much lower.
INPUT="${GRINS_TEST_INPUT_DIR}/penalty_poiseuille_stab.in"
DATA="${GRINS_TEST_DATA_DIR}/penalty_poiseuille_stab.xdr"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type lu -sub_pc_factor_shift_type nonzero"

${LIBMESH_RUN:-} $PROG input=$INPUT soln-data=$DATA vars='u v p' norms='L2 H1' tol='5.0e-8' $PETSC_OPTIONS
