
AM_TESTS_ENVIRONMENT  = GRINS_TEST_DIR=$(top_builddir)/test; export GRINS_TEST_DIR;
AM_TESTS_ENVIRONMENT += GRINS_BUILDSRC_DIR=$(top_builddir)/src; export GRINS_BUILDSRC_DIR;
AM_TESTS_ENVIRONMENT += GRINS_TEST_SRCDIR_DIR=$(top_srcdir)/test; export GRINS_TEST_SRCDIR_DIR;
AM_TESTS_ENVIRONMENT += GRINS_TEST_INPUT_DIR=$(top_srcdir)/test/input_files; export GRINS_TEST_INPUT_DIR;
AM_TESTS_ENVIRONMENT += GRINS_TEST_DATA_DIR=$(top_srcdir)/test/test_data; export GRINS_TEST_DATA_DIR;
AM_TESTS_ENVIRONMENT += GRINS_LIBMESH_DIM=@LIBMESH_DIM@; export GRINS_LIBMESH_DIM;

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/src/bc_handling/include
AM_CPPFLAGS += -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/common/include
AM_CPPFLAGS += -I$(top_srcdir)/src/error_estimation/include
AM_CPPFLAGS += -I$(top_srcdir)/src/physics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/properties/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include
AM_CPPFLAGS += -I$(top_srcdir)/src/solver/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_srcdir)/src/variables/include
AM_CPPFLAGS += -I$(top_srcdir)/src/visualization/include

AM_CPPFLAGS += -I$(top_srcdir)/test/interface
AM_CPPFLAGS += -I$(top_builddir)/test/common
AM_CPPFLAGS += -I$(top_srcdir)/test/common
AM_CPPFLAGS += -I$(top_builddir)/test/common

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)

LIBS =
LIBS += $(LIBGRINS_LIBS)

AM_LDFLAGS =

#----------------
# Cantera support
#----------------
if CANTERA_ENABLED
   AM_CPPFLAGS += $(CANTERA_CPPFLAGS)
   AM_LDFLAGS += $(CANTERA_LDFLAGS)
   LIBS += $(CANTERA_LIBS)
endif

#----------------
# Antioch support
#----------------
if ANTIOCH_ENABLED
   AM_CPPFLAGS += $(ANTIOCH_CPPFLAGS)
endif

#--------------
# MASA support
#--------------
if MASA_ENABLED
  AM_CPPFLAGS += $(MASA_CXXFLAGS)
  LIBS += $(MASA_LIBS)
endif

#----------------
# CppUnit support
#----------------
if GRINS_ENABLE_CPPUNIT
  AM_CPPFLAGS += $(CPPUNIT_CPPFLAGS)
  LIBS += $(CPPUNIT_LIBS)
endif

# These are programs that need to be built in order to run
# the various tests. They will be built with `make check`
check_PROGRAMS =

# Unit Tests
check_PROGRAMS += antioch_mixture
check_PROGRAMS += arrhenius_catalycity
check_PROGRAMS += cantera_mixture
check_PROGRAMS += composite_function
check_PROGRAMS += constant_catalycity
check_PROGRAMS += gas_recombination_catalytic_wall
check_PROGRAMS += gas_solid_catalytic_wall
check_PROGRAMS += power_law_catalycity
check_PROGRAMS += split_string

# Interface Tests
check_PROGRAMS += interface_cppunit_test
check_PROGRAMS += antioch_evaluator_regression
check_PROGRAMS += antioch_mixture_averaged_transport_evaluator_regression
check_PROGRAMS += cantera_chem_thermo_test
check_PROGRAMS += cantera_evaluator_regression
check_PROGRAMS += cantera_transport_regression

# Exact solution tests
check_PROGRAMS += generic_exact_solution_testing_app
check_PROGRAMS += vorticity_qoi

# Regression Tests
check_PROGRAMS += test_thermally_driven_flow
check_PROGRAMS += test_turbulent_channel
check_PROGRAMS += gaussian_profiles
check_PROGRAMS += low_mach_cavity_benchmark_regression
check_PROGRAMS += grins_flow_regression
check_PROGRAMS += reacting_low_mach_regression
check_PROGRAMS += elastic_sheet_regression
check_PROGRAMS += 3d_low_mach_jacobians_xy
check_PROGRAMS += 3d_low_mach_jacobians_xz
check_PROGRAMS += 3d_low_mach_jacobians_yz
check_PROGRAMS += suspended_cable_regression
check_PROGRAMS += generic_solution_regression
check_PROGRAMS += axisym_reacting_low_mach_regression


# AMR Tests
check_PROGRAMS += generic_amr_testing_app

# Unit test source files
antioch_mixture_SOURCES = unit/antioch_mixture.C
arrhenius_catalycity_SOURCES = unit/arrhenius_catalycity.C
cantera_mixture_SOURCES = unit/cantera_mixture.C
composite_function_SOURCES = unit/composite_function.C
constant_catalycity_SOURCES = unit/constant_catalycity.C
gas_recombination_catalytic_wall_SOURCES = unit/gas_recombination_catalytic_wall.C
gas_solid_catalytic_wall_SOURCES = unit/gas_solid_catalytic_wall.C
power_law_catalycity_SOURCES = unit/power_law_catalycity.C
split_string_SOURCES = unit/split_string.C

# Interface test source files
interface_cppunit_test_SOURCES = interface/interface_cppunit_test.C \
                                 interface/air_5sp.C

antioch_evaluator_regression_SOURCES = interface/antioch_evaluator_regression.C
antioch_mixture_averaged_transport_evaluator_regression_SOURCES = interface/antioch_mixture_averaged_transport_evaluator_regression.C
cantera_chem_thermo_test_SOURCES = interface/cantera_chem_thermo_test.C
cantera_evaluator_regression_SOURCES = interface/cantera_evaluator_regression.C
cantera_transport_regression_SOURCES = interface/cantera_transport_regression.C


# Exact solution test source files
generic_exact_solution_testing_app_SOURCES = exact_soln/generic_exact_solution_testing_app.C
vorticity_qoi_SOURCES = exact_soln/vorticity_qoi.C


# Regression test source files
test_thermally_driven_flow_SOURCES = regression/test_thermally_driven_flow.C
test_turbulent_channel_SOURCES = regression/test_turbulent_channel.C
gaussian_profiles_SOURCES = regression/gaussian_profiles.C
low_mach_cavity_benchmark_regression_SOURCES = regression/low_mach_cavity_benchmark_regression.C
grins_flow_regression_SOURCES = regression/grins_flow_regression.C
reacting_low_mach_regression_SOURCES = regression/reacting_low_mach_regression.C
elastic_sheet_regression_SOURCES = regression/elastic_sheet_regression.C
3d_low_mach_jacobians_xy_SOURCES = regression/3d_low_mach_jacobians.C
3d_low_mach_jacobians_xz_SOURCES = regression/3d_low_mach_jacobians.C
3d_low_mach_jacobians_yz_SOURCES = regression/3d_low_mach_jacobians.C
suspended_cable_regression_SOURCES = regression/suspended_cable_regression.C
generic_solution_regression_SOURCES = regression/generic_solution_regression.C
axisym_reacting_low_mach_regression_SOURCES = regression/axisym_reacting_low_mach_regression.C

# AMR test source files
generic_amr_testing_app_SOURCES = amr/generic_amr_testing_app.C


#Define tests to actually be run
# The .sh are in the <subdir> directory but the programs
# get built directly in the test directory, not test/<subdir>
# so they don't get a <subdir>/ prefix
TESTS =

# These are tests that are *expected* to fail
XFAIL_TESTS =

# Unit TESTS
TESTS += unit/antioch_mixture.sh
TESTS += arrhenius_catalycity
TESTS += unit/cantera_mixture.sh
TESTS += composite_function
TESTS += constant_catalycity
TESTS += unit/gas_recombination_catalytic_wall_antioch.sh
TESTS += unit/gas_recombination_catalytic_wall_cantera.sh
TESTS += unit/gas_solid_catalytic_wall_antioch.sh
TESTS += unit/gas_solid_catalytic_wall_cantera.sh
TESTS += power_law_catalycity
TESTS += split_string

# Interface TESTS
TESTS += interface_cppunit_test
TESTS += interface/antioch_evaluator_regression.sh
TESTS += interface/antioch_mixture_averaged_transport_evaluator_regression.sh
TESTS += interface/cantera_chem_thermo_test.sh
TESTS += interface/cantera_evaluator_regression.sh
TESTS += interface/cantera_transport_regression.sh


# Exact solution TESTS
TESTS += exact_soln/axi_ns_con_cyl_flow.sh
TESTS += exact_soln/axi_ns_poiseuille_flow.sh
TESTS += exact_soln/convection_diffusion_steady_1d.sh
TESTS += exact_soln/convection_diffusion_unsteady_2d.sh
TESTS += exact_soln/laplace_parsed_source.sh
TESTS += exact_soln/ns_couette_flow_2d_x.sh
TESTS += exact_soln/ns_couette_flow_2d_y.sh
TESTS += exact_soln/ns_poiseuille_flow.sh
TESTS += exact_soln/stokes_poiseuille_flow.sh
TESTS += exact_soln/stokes_poiseuille_flow_parsed_viscosity.sh
TESTS += exact_soln/stokes_poiseuille_flow_parsed_viscosity_parsed_conductivity.sh
TESTS += exact_soln/vorticity_qoi.sh


# Error/xfail TESTS
TESTS += error_xfail/error_ufo.sh
XFAIL_TESTS += error_xfail/error_ufo.sh

TESTS += error_xfail/stokes_invalid_pin_location.sh
XFAIL_TESTS += error_xfail/stokes_invalid_pin_location.sh

TESTS += error_xfail/ins_invalid_pin_location.sh
XFAIL_TESTS += error_xfail/ins_invalid_pin_location.sh

TESTS += error_xfail/lmns_invalid_pin_location.sh
XFAIL_TESTS += error_xfail/lmns_invalid_pin_location.sh

TESTS += error_xfail/viscosity_options_conflict_1.sh
XFAIL_TESTS += error_xfail/viscosity_options_conflict_1.sh

TESTS += error_xfail/viscosity_options_conflict_2.sh
XFAIL_TESTS += error_xfail/viscosity_options_conflict_2.sh

TESTS += error_xfail/viscosity_options_conflict_3.sh
XFAIL_TESTS += error_xfail/viscosity_options_conflict_3.sh

# Regression TESTS
TESTS += gaussian_profiles
TESTS += regression/warn_only_ufo.sh
TESTS += regression/poisson_weighted_flux.sh
TESTS += regression/sa_2d_turbulent_channel.sh
TESTS += regression/thermally_driven_2d_flow.sh
TESTS += regression/axi_thermally_driven_flow.sh
TESTS += regression/thermally_driven_3d_flow.sh
TESTS += regression/convection_cell.sh
TESTS += regression/convection_cell_parsed.sh
TESTS += regression/2d_pseudofan.sh
TESTS += regression/2d_pseudoprop.sh
TESTS += regression/2d_fantrick.sh
TESTS += regression/2d_proptrick.sh
TESTS += regression/2d_vsource.sh
TESTS += regression/dirichlet_fem.sh
TESTS += regression/dirichlet_nan.sh
TESTS += regression/simple_ode.sh
TESTS += regression/parsed_qoi.sh
TESTS += regression/low_mach_cavity_benchmark.sh
TESTS += regression/backward_facing_step.sh
TESTS += regression/locally_refine.sh
TESTS += regression/penalty_poiseuille.sh
TESTS += regression/penalty_poiseuille_stab.sh
TESTS += regression/redistribute.sh
TESTS += regression/coupled_stokes_ns.sh

TESTS += regression/reacting_low_mach_cantera.sh
XFAIL_TESTS += regression/reacting_low_mach_cantera.sh

TESTS += regression/reacting_low_mach_antioch_statmech_constant.sh
TESTS += regression/reacting_low_mach_antioch_statmech_constant_prandtl.sh
TESTS += regression/reacting_low_mach_antioch_cea_constant.sh
TESTS += regression/axisym_reacting_low_mach_antioch_cea_constant.sh
TESTS += regression/reacting_low_mach_antioch_cea_constant_prandtl.sh
TESTS += regression/reacting_low_mach_antioch_statmech_blottner_eucken_lewis.sh
TESTS += regression/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_catalytic_wall.sh
TESTS += regression/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_arrhenius_catalytic_wall.sh
TESTS += regression/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_power_catalytic_wall.sh
TESTS += regression/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall.sh
TESTS += regression/reacting_low_mach_antioch_kinetics_theory.sh
TESTS += regression/3d_low_mach_jacobians_xy.sh
TESTS += regression/3d_low_mach_jacobians_xz.sh
TESTS += regression/3d_low_mach_jacobians_yz.sh
TESTS += regression/elastic_mooney_rivlin_sheet.sh
TESTS += regression/elastic_mooney_rivlin_inflating_sheet.sh
TESTS += regression/suspended_cable_test.sh
TESTS += regression/elastic_mooney_rivlin_circle_hookean_stiffeners.sh
TESTS += regression/elastic_mooney_rivlin_square_hookean_stiffeners.sh
TESTS += regression/elastic_mooney_rivlin_square_hookean_stiffeners_unifref.sh
TESTS += regression/elastic_mooney_rivlin_square_hookean_stiffeners_newmark.sh
TESTS += regression/elastic_mooney_rivlin_square_hookean_stiffeners_rayleigh_damping_newmark.sh
TESTS += regression/extra_quadrature_order_laplace_arefee_amr.sh


#AMR TESTS
TESTS += amr/convection_diffusion_unsteady_2d.sh

CLEANFILES =

if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

# If we have tests fail, libmesh_error() may leave output
CLEANFILES += traceout*.txt temp_print_trace.*

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
