BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = grins grins_version 

if CANTERA_ENABLED
   bin_PROGRAMS += cantera_kinetic_rates
endif

lib_LTLIBRARIES = libgrins.la

#--------------------------------------------------
# List of directories with source and include files
#--------------------------------------------------
solver_srcdir = $(top_srcdir)/src/solver/src
solver_incdir = $(top_srcdir)/src/solver/include

physics_srcdir = $(top_srcdir)/src/physics/src
physics_incdir = $(top_srcdir)/src/physics/include

boundary_conditions_srcdir = $(top_srcdir)/src/boundary_conditions/src
boundary_conditions_incdir = $(top_srcdir)/src/boundary_conditions/include

utilities_srcdir = $(top_srcdir)/src/utilities/src
utilities_incdir = $(top_srcdir)/src/utilities/include

visualization_srcdir = $(top_srcdir)/src/visualization/src
visualization_incdir = $(top_srcdir)/src/visualization/include

matprops_srcdir = $(top_srcdir)/src/mat_props/src
matprops_incdir = $(top_srcdir)/src/mat_props/include

bc_handling_srcdir = $(top_srcdir)/src/bc_handling/src
bc_handling_incdir = $(top_srcdir)/src/bc_handling/include

qoi_srcdir = $(top_srcdir)/src/qoi/src
qoi_incdir = $(top_srcdir)/src/qoi/include

kinetics_srcdir = $(top_srcdir)/src/kinetics/src
kinetics_incdir = $(top_srcdir)/src/kinetics/include

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
# Sources in src/grins
solver_SRC   = $(solver_srcdir)/grins_solver.C \
               $(solver_srcdir)/mesh_builder.C \
               $(solver_srcdir)/simulation.C \
               $(solver_srcdir)/solver_factory.C \
               $(solver_srcdir)/grins_steady_solver.C \
               $(solver_srcdir)/grins_unsteady_solver.C \
               $(solver_srcdir)/simulation_builder.C \
               $(solver_srcdir)/solver_context.C

physics_SRC  = $(physics_srcdir)/multiphysics_sys.C \
               $(physics_srcdir)/physics.C \
               $(physics_srcdir)/stokes.C \
               $(physics_srcdir)/inc_navier_stokes_base.C \
               $(physics_srcdir)/inc_navier_stokes.C \
               $(physics_srcdir)/inc_navier_stokes_stab_base.C \
               $(physics_srcdir)/inc_navier_stokes_adjoint_stab.C \
               $(physics_srcdir)/axisym_inc_navier_stokes.C \
               $(physics_srcdir)/heat_transfer_base.C \
               $(physics_srcdir)/heat_transfer.C \
               $(physics_srcdir)/heat_transfer_source.C \
               $(physics_srcdir)/heat_transfer_stab_base.C \
               $(physics_srcdir)/heat_transfer_adjoint_stab.C \
               $(physics_srcdir)/axisym_heat_transfer.C \
               $(physics_srcdir)/boussinesq_buoyancy.C \
               $(physics_srcdir)/axisym_boussinesq_buoyancy.C \
               $(physics_srcdir)/physics_factory.C \
               $(physics_srcdir)/low_mach_navier_stokes_base.C \
               $(physics_srcdir)/low_mach_navier_stokes.C \
               $(physics_srcdir)/low_mach_navier_stokes_stab_base.C \
               $(physics_srcdir)/low_mach_navier_stokes_spgsm_stab.C \
               $(physics_srcdir)/low_mach_navier_stokes_vms_stab.C \
               $(physics_srcdir)/low_mach_navier_stokes_braack_stab.C \
               $(physics_srcdir)/stab_helper.C \
               $(physics_srcdir)/inc_navier_stokes_stab_helper.C \
               $(physics_srcdir)/heat_transfer_stab_helper.C \
               $(physics_srcdir)/low_mach_navier_stokes_stab_helper.C \
               $(physics_srcdir)/reacting_low_mach_navier_stokes_base.C \
               $(physics_srcdir)/reacting_low_mach_navier_stokes.C \
               $(physics_srcdir)/heat_conduction.C

boundary_conditions_SRC = $(boundary_conditions_srcdir)/boundary_conditions.C \
                          $(boundary_conditions_srcdir)/bc_factory.C \
                          $(boundary_conditions_srcdir)/neumann_func_obj.C \
                          $(boundary_conditions_srcdir)/parabolic_profile.C \
                          $(boundary_conditions_srcdir)/gaussian_profile.C \
                          $(boundary_conditions_srcdir)/gaussian_xy_profile.C \
                          $(boundary_conditions_srcdir)/concentric_cylinder_profile.C \
                          $(boundary_conditions_srcdir)/dbc_container.C \
                          $(boundary_conditions_srcdir)/pressure_pinning.C \
                          $(boundary_conditions_srcdir)/pbc_container.C \
                          $(boundary_conditions_srcdir)/nbc_container.C \
                          $(boundary_conditions_srcdir)/catalytic_wall.C \
                          $(boundary_conditions_srcdir)/catalytic_wall_helper.C

utilities_SRC = $(utilities_srcdir)/grins_version.C \
                $(utilities_srcdir)/input_utils.C \
                $(utilities_srcdir)/cached_values.C \
                $(utilities_srcdir)/tinyxml2.cpp

visualization_SRC = $(visualization_srcdir)/steady_visualization.C \
                    $(visualization_srcdir)/unsteady_visualization.C \
                    $(visualization_srcdir)/visualization.C \
                    $(visualization_srcdir)/visualization_factory.C \
                    $(visualization_srcdir)/postprocessed_quantities.C \
                    $(visualization_srcdir)/postprocessing_factory.C

matprops_SRC = $(matprops_srcdir)/grins_constant_function.C \
               $(matprops_srcdir)/constant_conductivity.C   \
               $(matprops_srcdir)/constant_viscosity.C      \
               $(matprops_srcdir)/constant_specific_heat.C  \
               $(matprops_srcdir)/constant_source_func.C    \
               $(matprops_srcdir)/chemical_species.C        \
               $(matprops_srcdir)/chemical_mixture.C        \
               $(matprops_srcdir)/cea_curve_fit.C           \
               $(matprops_srcdir)/cea_thermo.C              \
               $(matprops_srcdir)/blottner_viscosity.C      \
               $(matprops_srcdir)/cantera_singleton.C       \
               $(matprops_srcdir)/cantera_thermo.C          \
               $(matprops_srcdir)/cantera_transport.C       \
               $(matprops_srcdir)/cantera_kinetics.C        \
               $(matprops_srcdir)/ideal_gas_mixture.C       \
               $(matprops_srcdir)/constant_transport.C

bc_handling_SRC = $(bc_handling_srcdir)/bc_handling_base.C \
                  $(bc_handling_srcdir)/heat_transfer_bc_handling.C \
                  $(bc_handling_srcdir)/axisym_heat_transfer_bc_handling.C \
                  $(bc_handling_srcdir)/inc_navier_stokes_bc_handling.C \
                  $(bc_handling_srcdir)/axisym_inc_navier_stokes_bc_handling.C \
                  $(bc_handling_srcdir)/low_mach_navier_stokes_bc_handling.C \
                  $(bc_handling_srcdir)/reacting_low_mach_navier_stokes_bc_handling.C

qoi_SRC = $(qoi_srcdir)/average_nusselt_number.C \
          $(qoi_srcdir)/qoi_base.C \
          $(qoi_srcdir)/qoi_factory.C \
          $(qoi_srcdir)/vorticity.C

kinetics_SRC = $(kinetics_srcdir)/arrhenius_rate.C \
               $(kinetics_srcdir)/reaction.C \
               $(kinetics_srcdir)/reaction_set.C \
               $(kinetics_srcdir)/grins_kinetics.C

sfilesdir = $(prefix)/src
sfiles_DATA = $(physics_srcdir)/axisym_heat_transfer.C

#----------------------------
#INCLUDES we want distributed
#----------------------------
solver_INC   = $(solver_incdir)/grins/grins_solver.h \
               $(solver_incdir)/grins/mesh_builder.h \
               $(solver_incdir)/grins/simulation.h \
               $(solver_incdir)/grins/solver_factory.h \
               $(solver_incdir)/grins/grins_steady_solver.h \
               $(solver_incdir)/grins/grins_unsteady_solver.h \
               $(solver_incdir)/grins/simulation_builder.h \
               $(solver_incdir)/grins/solver_context.h

physics_INC  = $(physics_incdir)/grins/multiphysics_sys.h \
               $(physics_incdir)/grins/physics.h \
               $(physics_incdir)/grins/variable_name_defaults.h \
               $(physics_incdir)/grins/var_typedefs.h \
               $(physics_incdir)/grins/axisym_inc_navier_stokes.h \
               $(physics_incdir)/grins/stokes.h \
               $(physics_incdir)/grins/inc_navier_stokes_base.h \
               $(physics_incdir)/grins/inc_navier_stokes.h \
               $(physics_incdir)/grins/inc_navier_stokes_stab_base.h \
               $(physics_incdir)/grins/inc_navier_stokes_adjoint_stab.h \
               $(physics_incdir)/grins/heat_transfer_base.h \
               $(physics_incdir)/grins/heat_transfer.h \
               $(physics_incdir)/grins/heat_transfer_source.h \
               $(physics_incdir)/grins/heat_transfer_stab_base.h \
               $(physics_incdir)/grins/heat_transfer_adjoint_stab.h \
               $(physics_incdir)/grins/axisym_heat_transfer.h \
               $(physics_incdir)/grins/boussinesq_buoyancy.h \
               $(physics_incdir)/grins/axisym_boussinesq_buoyancy.h \
               $(physics_incdir)/grins/physics_factory.h \
               $(physics_incdir)/grins/grins_physics_names.h \
               $(physics_incdir)/grins/low_mach_navier_stokes_base.h \
               $(physics_incdir)/grins/low_mach_navier_stokes.h \
               $(physics_incdir)/grins/low_mach_navier_stokes_stab_base.h \
               $(physics_incdir)/grins/low_mach_navier_stokes_spgsm_stab.h \
               $(physics_incdir)/grins/low_mach_navier_stokes_vms_stab.h \
               $(physics_incdir)/grins/low_mach_navier_stokes_braack_stab.h \
               $(physics_incdir)/grins/stab_helper.h \
               $(physics_incdir)/grins/inc_navier_stokes_stab_helper.h \
               $(physics_incdir)/grins/heat_transfer_stab_helper.h \
               $(physics_incdir)/grins/low_mach_navier_stokes_stab_helper.h \
               $(physics_incdir)/grins/reacting_low_mach_navier_stokes_base.h \
               $(physics_incdir)/grins/reacting_low_mach_navier_stokes.h \
               $(physics_incdir)/grins/heat_conduction.h

boundary_conditions_INC = $(boundary_conditions_incdir)/grins/bc_types.h \
                          $(boundary_conditions_incdir)/grins/boundary_conditions.h \
                          $(boundary_conditions_incdir)/grins/bc_factory.h \
                          $(boundary_conditions_incdir)/grins/neumann_func_obj.h \
                          $(boundary_conditions_incdir)/grins/parabolic_profile.h \
                          $(boundary_conditions_incdir)/grins/gaussian_profile.h \
                          $(boundary_conditions_incdir)/grins/gaussian_xy_profile.h \
                          $(boundary_conditions_incdir)/grins/concentric_cylinder_profile.h \
                          $(boundary_conditions_incdir)/grins/dbc_container.h \
                          $(boundary_conditions_incdir)/grins/pressure_pinning.h \
                          $(boundary_conditions_incdir)/grins/pbc_container.h \
                          $(boundary_conditions_incdir)/grins/nbc_container.h \
                          $(boundary_conditions_incdir)/grins/catalytic_wall.h \
                          $(boundary_conditions_incdir)/grins/catalytic_wall_helper.h

utilities_INC = $(top_builddir)/src/utilities/include/grins/grins_version.h \
                $(utilities_incdir)/grins/physical_constants.h \
                $(utilities_incdir)/grins/input_utils.h \
                $(utilities_incdir)/grins/math_constants.h \
                $(utilities_incdir)/grins/cached_values.h \
                $(utilities_incdir)/grins/cached_quantities_enum.h \
                $(utilities_incdir)/grins/tinyxml2.h \
                $(utilities_incdir)/grins/string_utils.h

visualization_INC = $(visualization_incdir)/grins/steady_visualization.h \
                    $(visualization_incdir)/grins/unsteady_visualization.h \
                    $(visualization_incdir)/grins/visualization.h \
                    $(visualization_incdir)/grins/visualization_factory.h \
                    $(visualization_incdir)/grins/postprocessed_quantities.h \
                    $(visualization_incdir)/grins/postprocessing_factory.h

matprops_INC = $(matprops_incdir)/grins/grins_constant_function.h \
               $(matprops_incdir)/grins/constant_conductivity.h   \
               $(matprops_incdir)/grins/constant_viscosity.h      \
               $(matprops_incdir)/grins/constant_specific_heat.h  \
               $(matprops_incdir)/grins/constant_source_func.h    \
               $(matprops_incdir)/grins/species_enum.h            \
               $(matprops_incdir)/grins/chemical_species.h        \
               $(matprops_incdir)/grins/chemical_mixture.h        \
               $(matprops_incdir)/grins/cea_curve_fit.h           \
               $(matprops_incdir)/grins/cea_thermo.h              \
               $(matprops_incdir)/grins/blottner_viscosity.h      \
               $(matprops_incdir)/grins/cantera_singleton.h       \
               $(matprops_incdir)/grins/cantera_thermo.h          \
               $(matprops_incdir)/grins/cantera_transport.h       \
               $(matprops_incdir)/grins/cantera_kinetics.h        \
               $(matprops_incdir)/grins/ideal_gas_mixture.h       \
               $(matprops_incdir)/grins/constant_transport.h

bc_handling_INC = $(bc_handling_incdir)/grins/bc_handling_base.h                      \
                  $(bc_handling_incdir)/grins/heat_transfer_bc_handling.h             \
                  $(bc_handling_incdir)/grins/axisym_heat_transfer_bc_handling.h      \
                  $(bc_handling_incdir)/grins/inc_navier_stokes_bc_handling.h         \
                  $(bc_handling_incdir)/grins/axisym_inc_navier_stokes_bc_handling.h  \
                  $(bc_handling_incdir)/grins/low_mach_navier_stokes_bc_handling.h    \
                  $(bc_handling_incdir)/grins/reacting_low_mach_navier_stokes_bc_handling.h

qoi_INC = $(qoi_incdir)/grins/average_nusselt_number.h \
          $(qoi_incdir)/grins/qoi_base.h \
          $(qoi_incdir)/grins/qoi_factory.h \
          $(qoi_incdir)/grins/qoi_names.h \
          $(qoi_incdir)/grins/vorticity.h

kinetics_INC = $(kinetics_incdir)/grins/arrhenius_rate.h \
               $(kinetics_incdir)/grins/reaction.h \
               $(kinetics_incdir)/grins/reaction_set.h \
               $(kinetics_incdir)/grins/grins_kinetics.h \
               $(kinetics_incdir)/grins/reaction_enum.h

includedir = $(prefix)/include/grins
include_HEADERS = $(solver_INC) $(physics_INC) $(boundary_conditions_INC) \
                  $(utilities_INC) $(visualization_INC) $(matprops_INC) \
                  $(bc_handling_INC) $(qoi_INC) $(kinetics_INC)

libgrins_la_SOURCES = $(solver_SRC) \
                      $(physics_SRC) \
                      $(boundary_conditions_SRC) \
                      $(utilities_SRC) \
                      $(visualization_SRC) \
                      $(matprops_SRC) \
                      $(bc_handling_SRC) \
                      $(qoi_SRC) \
                      $(kinetics_SRC)

if LIBMESH_LIBTOOL
   libgrins_la_LIBADD = $(LIBMESH_LIBDIR)/libmesh_$(LIBMESH_METHOD).la
else
   libgrins_la_LIBADD = $(LIBMESH_LDFLAGS)
endif

if GRVY_ENABLED
   libgrins_la_LIBADD += $(GRVY_PREFIX)/lib/libgrvy.la
endif

#------------------------------------
# Source for stand alone GRINS driver
#------------------------------------

grins_SOURCES = $(solver_srcdir)/grins.C
grins_LDADD = libgrins.la
if !LIBMESH_LIBTOOL
grins_LDADD += $(LIBMESH_LDFLAGS)
endif

grins_version_SOURCES = $(utilities_srcdir)/version.C
grins_version_LDADD = libgrins.la
if !LIBMESH_LIBTOOL
grins_version_LDADD += $(LIBMESH_LDFLAGS)
endif

if CANTERA_ENABLED
   cantera_kinetic_rates_SOURCES = $(top_srcdir)/src/apps/cantera_kinetic_rates.C
   cantera_kinetic_rates_LDADD = libgrins.la
endif

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(solver_incdir) \
              -I$(physics_incdir) \
              -I$(boundary_conditions_incdir) \
              -I$(utilities_incdir) \
              -I$(visualization_incdir) \
              -I$(matprops_incdir) \
              -I$(bc_handling_incdir) \
              -I$(qoi_incdir) \
              -I$(kinetics_incdir) \
              -I$(top_builddir)/src/utilities/include #For grins_version.h

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_CPPFLAGS) \
               $(GRVY_CFLAGS)

#----------------
# Cantera support
#----------------
if CANTERA_ENABLED
   AM_CPPFLAGS += $(CANTERA_CPPFLAGS)
   libgrins_la_LIBADD  += $(CANTERA_LDFLAGS)
endif

#-------------
# MASA support
#-------------
if MASA_ENABLED
  AM_CPPFLAGS += $(MASA_CXXFLAGS)
  AM_LDFLAGS  = $(MASA_LIBS)
endif

#--------------
# Code coverage
#--------------
if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(solver_SRC) $(solver_INC) \
	$(generic_driver_SOURCES) $(physics_SRC) $(physics_INC) \
	$(boundary_conditions_INC) $(boundary_conditions_SRC) \
        $(utilities_INC) $(utilities_SRC) \
	$(visualization_INC) $(visualization_SRC) \
        $(matprops_INC) $(matprops_SRC) \
        $(bc_handling_INC) $(bc_handling_SRC) \
        $(solver_SOURCES) $(grins_version_SOURCES) \
        $(qoi_SRC) $(kinetics_SRC)
	echo 'updated source license headers' >$@

EXTRA_DIST = common

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
