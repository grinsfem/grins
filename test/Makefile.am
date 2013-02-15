BUILT_SOURCES   = .license.stamp

check_PROGRAMS  = test_ns_couette_flow_2d_x
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
check_PROGRAMS += cea_thermo_test
check_PROGRAMS += cantera_chem_thermo_test
check_PROGRAMS += cantera_transport_test
check_PROGRAMS += ideal_gas_mixture_cantera
check_PROGRAMS += reacting_low_mach_regression
check_PROGRAMS += arrhenius_rate_test
check_PROGRAMS += grins_kinetics_regression
check_PROGRAMS += catalytic_wall_unit

AM_CPPFLAGS = -I$(top_srcdir)/src/solver/include \
              -I$(top_srcdir)/src/physics/include \
              -I$(top_srcdir)/src/boundary_conditions/include \
              -I$(top_srcdir)/src/utilities/include \
              -I$(top_srcdir)/src/visualization/include \
              -I$(top_srcdir)/src/mat_props/include \
              -I$(top_srcdir)/src/bc_handling/include \
              -I$(top_srcdir)/src/qoi/include \
              -I$(top_srcdir)/src/kinetics/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_CPPFLAGS) \
               $(GRVY_CFLAGS)

AM_LDFLAGS = $(LIBGRINS_LIBS)

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
  AM_CPPFLAGS           += $(MASA_CXXFLAGS)
  AM_LDFLAGS            += $(MASA_LIBS)
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
cea_thermo_test_SOURCES = $(top_srcdir)/test/cea_thermo_test.C
cantera_chem_thermo_test_SOURCES = $(top_srcdir)/test/cantera_chem_thermo_test.C
cantera_transport_test_SOURCES = $(top_srcdir)/test/cantera_transport_test.C
ideal_gas_mixture_cantera_SOURCES = $(top_srcdir)/test/ideal_gas_mixture_cantera.C
reacting_low_mach_regression_SOURCES = $(top_srcdir)/test/reacting_low_mach_regression.C
arrhenius_rate_test_SOURCES = $(top_srcdir)/test/arrhenius_rate_test.C
grins_kinetics_regression_SOURCES = $(top_srcdir)/test/grins_kinetics_regression.C
catalytic_wall_unit_SOURCES = $(top_srcdir)/test/catalytic_wall_unit.C

# List of source files for license stamping
testsrc = $(test_ns_couette_flow_2d_x_SOURCES)       \
          $(test_ns_couette_flow_2d_y_SOURCES)       \
          $(test_stokes_poiseuille_flow_SOURCES)     \
          $(test_ns_poiseuille_flow_SOURCES)         \
          $(test_axi_ns_poiseuille_flow_SOURCES)     \
          $(test_axi_ns_con_cyl_flow_SOURCES)        \
          $(test_thermally_driven_flow_SOURCES)      \
          $(gaussian_profiles_SOURCES)               \
          $(vorticity_qoi_SOURCES)                   \
          $(low_mach_cavity_benchmark_regression_SOURCES) \
          $(backward_facing_step_regression_SOURCES) \
          $(reacting_low_mach_cantera_regression_SOURCES) \
          $(arrhenius_rate_test_SOURCES) \
          $(grins_kinetics_regression_SOURCES)

#Define tests to actually be run
TESTS  = gaussian_profiles
TESTS += chem_mixture
TESTS += cea_thermo_test
TESTS += catalytic_wall_unit
TESTS += cantera_chem_thermo_test.sh
TESTS += cantera_transport_test
TESTS += ideal_gas_mixture_cantera
TESTS += arrhenius_rate_test
TESTS += grins_kinetics_regression.sh

TESTS += $(top_builddir)/test/test_ns_couette_flow_2d_x.sh            \
         $(top_builddir)/test/test_ns_couette_flow_2d_y.sh            \
         $(top_builddir)/test/test_stokes_poiseuille_flow.sh          \
         $(top_builddir)/test/test_ns_poiseuille_flow.sh              \
         $(top_builddir)/test/test_axi_ns_poiseuille_flow.sh          \
         $(top_builddir)/test/test_axi_ns_con_cyl_flow.sh             \
         $(top_builddir)/test/test_thermally_driven_2d_flow.sh        \
         $(top_builddir)/test/test_axi_thermally_driven_flow.sh       \
         $(top_builddir)/test/test_thermally_driven_3d_flow.sh        \
         $(top_builddir)/test/test_vorticity_qoi.sh                   \
         $(top_builddir)/test/low_mach_cavity_benchmark_regression.sh \
         $(top_builddir)/test/backward_facing_step_regression.sh      \
         $(top_builddir)/test/coupled_stokes_ns.sh                    \
         $(top_builddir)/test/reacting_low_mach_cantera_regression.sh

TESTS += $(top_builddir)/test/reacting_low_mach_grins_kinetics_regression.sh
TESTS += $(top_builddir)/test/reacting_low_mach_catalytic_wall_grins_kinetics_regression.sh

shellfiles_src = test_ns_couette_flow_2d_x.sh       \
                 test_ns_couette_flow_2d_y.sh       \
                 test_stokes_poiseuille_flow.sh     \
                 test_ns_poiseuille_flow.sh         \
                 test_axi_ns_poiseuille_flow.sh     \
                 test_axi_ns_con_cyl_flow.sh        \
                 test_thermally_driven_2d_flow.sh   \
                 test_axi_thermally_driven_flow.sh  \
                 test_thermally_driven_3d_flow.sh   \
                 test_vorticity_qoi.sh              \
                 low_mach_cavity_benchmark_regression.sh \
                 backward_facing_step_regression.sh \
                 coupled_stokes_ns.sh \
                 cantera_chem_thermo_test.sh \
                 reacting_low_mach_cantera_regression.sh

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src) input_files test_data grids

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(testsrc)
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
