BUILT_SOURCES   = .license.stamp

CLEANFILES = .license.stamp

check_PROGRAMS =
check_PROGRAMS += test_ns_couette_flow_2d_x
check_PROGRAMS += test_ns_couette_flow_2d_y
check_PROGRAMS += test_stokes_poiseuille_flow
check_PROGRAMS += test_ns_poiseuille_flow
check_PROGRAMS += test_axi_ns_poiseuille_flow
check_PROGRAMS += test_axi_ns_con_cyl_flow
check_PROGRAMS += test_thermally_driven_flow
check_PROGRAMS += gaussian_profiles
check_PROGRAMS += vorticity_qoi
check_PROGRAMS += low_mach_cavity_benchmark_regression
check_PROGRAMS += grins_flow_regression
check_PROGRAMS += cantera_mixture_unit
check_PROGRAMS += cantera_chem_thermo_test
check_PROGRAMS += cantera_transport_regression
check_PROGRAMS += cantera_evaluator_regression
check_PROGRAMS += reacting_low_mach_regression
check_PROGRAMS += antioch_mixture_unit
check_PROGRAMS += antioch_kinetics_regression
check_PROGRAMS += antioch_evaluator_regression
check_PROGRAMS += antioch_wilke_evaluator_regression
check_PROGRAMS += composite_function_unit
check_PROGRAMS += gas_recombination_catalytic_wall_unit
check_PROGRAMS += gas_solid_catalytic_wall_unit
check_PROGRAMS += constant_catalycity_unit
check_PROGRAMS += arrhenius_catalycity_unit
check_PROGRAMS += power_law_catalycity_unit

AM_CPPFLAGS = 
AM_CPPFLAGS += -I$(top_srcdir)/src/bc_handling/include
AM_CPPFLAGS += -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/error_estimation/include
AM_CPPFLAGS += -I$(top_srcdir)/src/physics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/properties/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include
AM_CPPFLAGS += -I$(top_srcdir)/src/solver/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_srcdir)/src/visualization/include

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

# Sources for these tests
test_ns_couette_flow_2d_x_SOURCES = $(top_srcdir)/test/test_ns_couette_flow_2d_x.C
test_ns_couette_flow_2d_y_SOURCES = $(top_srcdir)/test/test_ns_couette_flow_2d_y.C
test_stokes_poiseuille_flow_SOURCES = $(top_srcdir)/test/test_stokes_poiseuille_flow.C
test_ns_poiseuille_flow_SOURCES = $(top_srcdir)/test/test_ns_poiseuille_flow.C
test_axi_ns_poiseuille_flow_SOURCES = $(top_srcdir)/test/test_axi_ns_poiseuille_flow.C
test_axi_ns_con_cyl_flow_SOURCES = $(top_srcdir)/test/test_axi_ns_con_cyl_flow.C
test_thermally_driven_flow_SOURCES = $(top_srcdir)/test/test_thermally_driven_flow.C
gaussian_profiles_SOURCES = $(top_srcdir)/test/gaussian_profiles.C
vorticity_qoi_SOURCES = $(top_srcdir)/test/test_vorticity_qoi.C
low_mach_cavity_benchmark_regression_SOURCES = $(top_srcdir)/test/low_mach_cavity_benchmark_regression.C
grins_flow_regression_SOURCES = $(top_srcdir)/test/grins_flow_regression.C
cantera_mixture_unit_SOURCES = $(top_srcdir)/test/cantera_mixture_unit.C
cantera_chem_thermo_test_SOURCES = $(top_srcdir)/test/cantera_chem_thermo_test.C
cantera_transport_regression_SOURCES = $(top_srcdir)/test/cantera_transport_regression.C
cantera_evaluator_regression_SOURCES = $(top_srcdir)/test/cantera_evaluator_regression.C
reacting_low_mach_regression_SOURCES = $(top_srcdir)/test/reacting_low_mach_regression.C
antioch_mixture_unit_SOURCES = $(top_srcdir)/test/antioch_mixture_unit.C
antioch_kinetics_regression_SOURCES = $(top_srcdir)/test/antioch_kinetics_regression.C
antioch_evaluator_regression_SOURCES = $(top_srcdir)/test/antioch_evaluator_regression.C
antioch_wilke_evaluator_regression_SOURCES = $(top_srcdir)/test/antioch_wilke_evaluator_regression.C
composite_function_unit_SOURCES = $(top_srcdir)/test/composite_function_unit.C
gas_recombination_catalytic_wall_unit_SOURCES = $(top_srcdir)/test/gas_recombination_catalytic_wall_unit.C
gas_solid_catalytic_wall_unit_SOURCES = $(top_srcdir)/test/gas_solid_catalytic_wall_unit.C
constant_catalycity_unit_SOURCES = $(top_srcdir)/test/constant_catalycity_unit.C
arrhenius_catalycity_unit_SOURCES = $(top_srcdir)/test/arrhenius_catalycity_unit.C
power_law_catalycity_unit_SOURCES = $(top_srcdir)/test/power_law_catalycity_unit.C

# List of source files for license stamping
STAMPED_FILES = 
STAMPED_FILES += $(test_ns_couette_flow_2d_x_SOURCES)
STAMPED_FILES += $(test_ns_couette_flow_2d_y_SOURCES)
STAMPED_FILES += $(test_stokes_poiseuille_flow_SOURCES)
STAMPED_FILES += $(test_ns_poiseuille_flow_SOURCES)
STAMPED_FILES += $(test_axi_ns_poiseuille_flow_SOURCES)
STAMPED_FILES += $(test_axi_ns_con_cyl_flow_SOURCES)
STAMPED_FILES += $(test_thermally_driven_flow_SOURCES)
STAMPED_FILES += $(gaussian_profiles_SOURCES)
STAMPED_FILES += $(vorticity_qoi_SOURCES)
STAMPED_FILES += $(low_mach_cavity_benchmark_regression_SOURCES)
STAMPED_FILES += $(grins_flow_regression_SOURCES)
STAMPED_FILES += $(cantera_mixture_unit_SOURCES)
STAMPED_FILES += $(cantera_chem_thermo_test_SOURCES)
STAMPED_FILES += $(cantera_transport_regression_SOURCES)
STAMPED_FILES += $(cantera_evaluator_regression_SOURCES)
STAMPED_FILES += $(reacting_low_mach_regression_SOURCES)
STAMPED_FILES += $(grins_kinetics_regression_SOURCES)
STAMPED_FILES += $(catalytic_wall_unit_SOURCES)
STAMPED_FILES += $(antioch_mixture_unit_SOURCES)
STAMPED_FILES += $(antioch_kinetics_regression_SOURCES)
STAMPED_FILES += $(antioch_evaluator_regression_SOURCES)
STAMPED_FILES += $(antioch_wilke_evaluator_regression_SOURCES)
STAMPED_FILES += $(composite_function_unit_SOURCES)
STAMPED_FILES += $(gas_recombination_catalytic_wall_unit_SOURCES)
STAMPED_FILES += $(gas_solid_catalytic_wall_unit_SOURCES)
STAMPED_FILES += $(constant_catalycity_unit_SOURCES)
STAMPED_FILES += $(arrhenius_catalycity_unit_SOURCES)
STAMPED_FILES += $(power_law_catalycity_unit_SOURCES)

#Define tests to actually be run
TESTS =
TESTS += gaussian_profiles
TESTS += cantera_mixture_unit.sh
TESTS += cantera_chem_thermo_test.sh
TESTS += cantera_transport_regression.sh
TESTS += cantera_evaluator_regression.sh
TESTS += antioch_mixture_unit.sh
TESTS += antioch_kinetics_regression.sh
TESTS += antioch_evaluator_regression.sh
TESTS += antioch_wilke_evaluator_regression.sh
TESTS += gas_recombination_catalytic_wall_unit_antioch.sh
TESTS += gas_recombination_catalytic_wall_unit_cantera.sh
TESTS += gas_solid_catalytic_wall_unit_antioch.sh
TESTS += gas_solid_catalytic_wall_unit_cantera.sh
TESTS += constant_catalycity_unit
TESTS += arrhenius_catalycity_unit
TESTS += power_law_catalycity_unit

TESTS += test_ns_couette_flow_2d_x.sh
TESTS += test_ns_couette_flow_2d_y.sh
TESTS += test_stokes_poiseuille_flow.sh
TESTS += test_ns_poiseuille_flow.sh
TESTS += test_axi_ns_poiseuille_flow.sh
TESTS += test_axi_ns_con_cyl_flow.sh
TESTS += test_thermally_driven_2d_flow.sh
TESTS += test_axi_thermally_driven_flow.sh
TESTS += test_thermally_driven_3d_flow.sh
TESTS += test_vorticity_qoi.sh
TESTS += low_mach_cavity_benchmark_regression.sh
TESTS += backward_facing_step_regression.sh
TESTS += locally_refine_regression.sh
TESTS += redistribute_regression.sh
TESTS += coupled_stokes_ns.sh

TESTS += reacting_low_mach_cantera_regression.sh
XFAIL_TESTS = reacting_low_mach_cantera_regression.sh

TESTS += reacting_low_mach_antioch_statmech_constant_regression.sh
TESTS += reacting_low_mach_antioch_statmech_constant_prandtl_regression.sh
TESTS += reacting_low_mach_antioch_cea_constant_regression.sh
TESTS += reacting_low_mach_antioch_cea_constant_prandtl_regression.sh
TESTS += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_regression.sh
TESTS += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_catalytic_wall_regression.sh
TESTS += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_arrhenius_catalytic_wall_regression.sh
TESTS += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_power_catalytic_wall_regression.sh
TESTS += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.sh
TESTS += composite_function_unit

shellfiles_src =
shellfiles_src += test_ns_couette_flow_2d_x.sh
shellfiles_src += test_ns_couette_flow_2d_y.sh
shellfiles_src += test_stokes_poiseuille_flow.sh
shellfiles_src += test_ns_poiseuille_flow.sh
shellfiles_src += test_axi_ns_poiseuille_flow.sh
shellfiles_src += test_axi_ns_con_cyl_flow.sh
shellfiles_src += test_thermally_driven_2d_flow.sh
shellfiles_src += test_axi_thermally_driven_flow.sh
shellfiles_src += test_thermally_driven_3d_flow.sh
shellfiles_src += test_vorticity_qoi.sh
shellfiles_src += low_mach_cavity_benchmark_regression.sh
shellfiles_src += backward_facing_step_regression.sh
shellfiles_src += locally_refine_regression.sh
shellfiles_src += redistribute_regression.sh
shellfiles_src += coupled_stokes_ns.sh
shellfiles_src += reacting_low_mach_cantera_regression.sh
shellfiles_src += cantera_mixture_unit.sh
shellfiles_src += cantera_chem_thermo_test.sh
shellfiles_src += cantera_transport_regression.sh
shellfiles_src += cantera_evaluator_regression.sh
shellfiles_src += antioch_mixture_unit.sh
shellfiles_src += antioch_kinetics_regression.sh
shellfiles_src += antioch_evaluator_regression.sh
shellfiles_src += antioch_wilke_evaluator_regression.sh
shellfiles_src += gas_recombination_catalytic_wall_unit_antioch.sh
shellfiles_src += gas_recombination_catalytic_wall_unit_cantera.sh
shellfiles_src += gas_solid_catalytic_wall_unit_antioch.sh
shellfiles_src += gas_solid_catalytic_wall_unit_cantera.sh
shellfiles_src += reacting_low_mach_antioch_statmech_constant_regression.sh
shellfiles_src += reacting_low_mach_antioch_statmech_constant_prandtl_regression.sh
shellfiles_src += reacting_low_mach_antioch_cea_constant_regression.sh
shellfiles_src += reacting_low_mach_antioch_cea_constant_prandtl_regression.sh
shellfiles_src += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_regression.sh
shellfiles_src += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_catalytic_wall_regression.sh
shellfiles_src += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_arrhenius_catalytic_wall_regression.sh
shellfiles_src += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_power_catalytic_wall_regression.sh
shellfiles_src += reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.sh

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src) input_files test_data grids

if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(STAMPED_FILES)
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
