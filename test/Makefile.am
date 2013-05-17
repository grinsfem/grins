BUILT_SOURCES   = .license.stamp

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
check_PROGRAMS += chem_mixture
check_PROGRAMS += cantera_chemistry_unit
check_PROGRAMS += cea_thermo_test
check_PROGRAMS += cantera_chem_thermo_test
check_PROGRAMS += cantera_transport_regression
check_PROGRAMS += ideal_gas_mixture_cantera_regression
check_PROGRAMS += reacting_low_mach_regression
check_PROGRAMS += arrhenius_rate_test
check_PROGRAMS += grins_kinetics_regression
check_PROGRAMS += catalytic_wall_unit

AM_CPPFLAGS = 
AM_CPPFLAGS += -I$(top_srcdir)/src/bc_handling/include
AM_CPPFLAGS += -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/error_estimation/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/mat_props/include
AM_CPPFLAGS += -I$(top_srcdir)/src/physics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/properties/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include
AM_CPPFLAGS += -I$(top_srcdir)/src/solver/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_srcdir)/src/visualization/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)

AM_LDFLAGS = 
AM_LDFLAGS += $(LIBGRINS_LIBS)

#----------------
# Cantera support
#----------------
if CANTERA_ENABLED
   AM_CPPFLAGS += $(CANTERA_CPPFLAGS)
   AM_LDFLAGS += $(CANTERA_LDFLAGS)
endif

#--------------
# MASA support
#--------------
if MASA_ENABLED
  AM_CPPFLAGS += $(MASA_CXXFLAGS)
  AM_LDFLAGS += $(MASA_LIBS)
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
chem_mixture_SOURCES = $(top_srcdir)/test/chem_mixture.C
cantera_chemistry_unit_SOURCES = $(top_srcdir)/test/cantera_chemistry_unit.C
cea_thermo_test_SOURCES = $(top_srcdir)/test/cea_thermo_test.C
cantera_chem_thermo_test_SOURCES = $(top_srcdir)/test/cantera_chem_thermo_test.C
cantera_transport_regression_SOURCES = $(top_srcdir)/test/cantera_transport_regression.C
ideal_gas_mixture_cantera_regression_SOURCES = $(top_srcdir)/test/ideal_gas_mixture_cantera_regression.C
reacting_low_mach_regression_SOURCES = $(top_srcdir)/test/reacting_low_mach_regression.C
arrhenius_rate_test_SOURCES = $(top_srcdir)/test/arrhenius_rate_test.C
grins_kinetics_regression_SOURCES = $(top_srcdir)/test/grins_kinetics_regression.C
catalytic_wall_unit_SOURCES = $(top_srcdir)/test/catalytic_wall_unit.C

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
STAMPED_FILES += $(chem_mixture_SOURCES)
STAMPED_FILES += $(cantera_chemistry_unit_SOURCES)
STAMPED_FILES += $(cea_thermo_test_SOURCES)
STAMPED_FILES += $(cantera_chem_thermo_test_SOURCES)
STAMPED_FILES += $(cantera_transport_regression_SOURCES)
STAMPED_FILES += $(ideal_gas_mixture_cantera_regression_SOURCES)
STAMPED_FILES += $(reacting_low_mach_regression_SOURCES)
STAMPED_FILES += $(arrhenius_rate_test_SOURCES)
STAMPED_FILES += $(grins_kinetics_regression_SOURCES)
STAMPED_FILES += $(catalytic_wall_unit_SOURCES)

#Define tests to actually be run
TESTS =
TESTS += gaussian_profiles
TESTS += chem_mixture
TESTS += cantera_chemistry_unit.sh
TESTS += cea_thermo_test
TESTS += catalytic_wall_unit
TESTS += cantera_chem_thermo_test.sh
TESTS += cantera_transport_regression.sh
TESTS += ideal_gas_mixture_cantera_regression.sh
TESTS += arrhenius_rate_test
TESTS += grins_kinetics_regression.sh

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
TESTS += coupled_stokes_ns.sh

TESTS += reacting_low_mach_cantera_regression.sh
TESTS += reacting_low_mach_grins_kinetics_regression.sh
TESTS += reacting_low_mach_catalytic_wall_grins_kinetics_regression.sh

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
shellfiles_src += coupled_stokes_ns.sh
shellfiles_src += cantera_chem_thermo_test.sh
shellfiles_src += reacting_low_mach_cantera_regression.sh
shellfiles_src += cantera_transport_regression.sh
shellfiles_src += cantera_chemistry_unit.sh
shellfiles_src += ideal_gas_mixture_cantera_regression.sh

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src) input_files test_data grids

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(STAMPED_FILES)
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
