#!/bin/bash
set -e

# Multigrid run script for regression test.
PROG="${GRINS_BUILDSRC_DIR}/grins"
INPUT="${GRINS_TEST_INPUT_DIR}/multigrid_stokes.in"

PETSC_OPTIONS="--use_petsc_dm --node-major-dofs \
-snes_view -snes_monitor -snes_converged_reason -snes_rtol 1.0e-4 \
-ksp_type fgmres -ksp_converged_reason -ksp_monitor_true_residual -ksp_rtol 1.0e-5 \
-pc_type fieldsplit -pc_fieldsplit_0_fields 0,1 -pc_fieldsplit_1_fields 2 \
-pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition a11 \
-fieldsplit_0_pc_type mg -fieldsplit_0_pc_mg_galerkin both -fieldsplit_0_pc_mg_type full \
-pc_mg_levels 3 -fieldsplit_0_pc_mg_levels 3 \
\
-fieldsplit_0_mg_levels_pc_type sor -fieldsplit_0_mg_levels_pc_sor_its 5 -fieldsplit_0_mg_levels_ksp_type richardson -fieldsplit_0_mg_levels_ksp_richardson_self_scale \
-fieldsplit_0_mg_levels_ksp_max_it 5 -fieldsplit_0_mg_levels_ksp_monitor_true_residual -fieldsplit_0_mg_levels_ksp_converged_reason \
\
-fieldsplit_0_mg_coarse_pc_type lu -fieldsplit_0_mg_coarse_ksp_type preonly \
\
-fieldsplit_0_ksp_max_it 1 -fieldsplit_0_ksp_type gmres -fieldsplit_0_ksp_monitor_true_residual -fieldsplit_0_ksp_converged_reason \
\
-fieldsplit_p_ksp_max_it 5 -fieldsplit_p_ksp_type gmres -fieldsplit_p_pc_type none -fieldsplit_p_ksp_rtol 1.0e-4 -fieldsplit_p_ksp_converged_reason -fieldsplit_p_ksp_monitor_true_residual \
\
-fieldsplit_p_inner_pc_type lu -fieldsplit_p_inner_ksp_type preonly \
-fieldsplit_p_upper_pc_type lu -fieldsplit_p_upper_ksp_type preonly"


# Gold data used for regression comparsion
GOLDDATA="${GRINS_TEST_DATA_DIR}/multigrid_stokes_gold.xdr"

# We dont pre-run the case with GRINS as we have to do it again from
# inside the app to get iteration counts from the solver.

# Filename to be generated since were running directly through
# regression_testing_app but want to maintain existing error checking
# when running app under standard paradigm
SOLNDATA="multigrid_stokes.xdr"
touch $SOLNDATA

# Now run the test part to make sure we're getting the correct thing
# No point in running through external app since we need to rerun
# internally anyway to collect iteration counts.
${GRINS_TEST_DIR}/regression_testing_app \
   input=$INPUT \
   vars='u v p' \
   norms='L2' \
   tol='1.0e-6' \
   nonlinear_its='1' \
   linear_its='13' \
   gold-data=$GOLDDATA \
   soln-data=$SOLNDATA \
   $PETSC_OPTIONS

# Remove generated test data
rm $SOLNDATA
