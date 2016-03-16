#!/bin/bash

PROG="${GRINS_TEST_DIR}/grins_flow_regression"

# Note the tolerance is high because of the H1 error in the pressure.
# If we had a per-variable tolerance the others could be much lower.
INPUT="${GRINS_TEST_INPUT_DIR}/penalty_poiseuille_stab.in ${GRINS_TEST_DATA_DIR}/penalty_poiseuille_stab.xdr 5.0e-8"

PETSC_OPTIONS="-pc_type asm -pc_asm_overlap 10 -sub_pc_type ilu -sub_pc_factor_levels 10 -sub_pc_factor_shift_type nonzero"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
