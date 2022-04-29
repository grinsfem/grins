#!/bin/bash
set -e

# Multigrid run script for regression test.
PROG="${GRINS_BUILDSRC_DIR}/grins"
INPUT="${GRINS_TEST_INPUT_DIR}/multigrid_poisson.in"

PETSC_OPTIONS="--use_petsc_dm --node_major_dofs \
-snes_atol 1.0e-6 -snes_view -snes_monitor -snes_converged_reason \
-pc_type mg -pc_mg_levels 3 -pc_mg_type full -pc_mg_galerkin both \
-mg_levels_ksp_type richardson -mg_levels_ksp_richardson_self_scale -mg_levels_ksp_monitor -mg_levels_ksp_converged_reason  \
-mg_coarse_ksp_type preonly -mg_coarse_pc_type lu \
-ksp_type richardson -ksp_richardson_self_scale -ksp_rtol 1.0e-8 -ksp_monitor_true_residual"

# Gold data used for regression comparsion
GOLDDATA="${GRINS_TEST_DATA_DIR}/multigrid_poisson_gold.xdr"

# We dont pre-run the case with GRINS as we have to do it again from
# inside the app to get iteration counts from the solver.

# Filename to be generated since were running directly through
# regression_testing_app but want to maintain existing error checking
# when running app under standard paradigm
SOLNDATA="multigrid_poisson.xda.gz"
touch $SOLNDATA

# Now run the test part to make sure we're getting the correct thing
# No point in running through external app since we need to rerun
# internally anyway to collect iteration counts.
${GRINS_TEST_DIR}/regression_testing_app \
   input=$INPUT \
   vars='Temp' \
   norms='L2' \
   tol='1.0e-6' \
   nonlinear_its='1' \
   linear_its='2' \
   gold-data=$GOLDDATA \
   soln-data=$SOLNDATA \
   $PETSC_OPTIONS

# Remove generated test data
rm $SOLNDATA
