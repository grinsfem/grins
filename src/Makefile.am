BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = grins grins_version
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
               $(solver_srcdir)/simulation_builder.C

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
               $(physics_srcdir)/low_mach_navier_stokes_stab_helper.C

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
                          $(boundary_conditions_srcdir)/nbc_container.C

utilities_SRC = $(utilities_srcdir)/grins_version.C

visualization_SRC = $(visualization_srcdir)/steady_visualization.C \
                    $(visualization_srcdir)/unsteady_visualization.C \
                    $(visualization_srcdir)/visualization.C \
                    $(visualization_srcdir)/visualization_factory.C

matprops_SRC = $(matprops_srcdir)/grins_constant_function.C \
               $(matprops_srcdir)/constant_conductivity.C   \
               $(matprops_srcdir)/constant_viscosity.C      \
               $(matprops_srcdir)/constant_specific_heat.C \
               $(matprops_srcdir)/constant_source_func.C

bc_handling_SRC = $(bc_handling_srcdir)/bc_handling_base.C \
                  $(bc_handling_srcdir)/heat_transfer_bc_handling.C \
                  $(bc_handling_srcdir)/axisym_heat_transfer_bc_handling.C \
                  $(bc_handling_srcdir)/inc_navier_stokes_bc_handling.C \
                  $(bc_handling_srcdir)/axisym_inc_navier_stokes_bc_handling.C \
                  $(bc_handling_srcdir)/low_mach_navier_stokes_bc_handling.C

qoi_SRC = $(qoi_srcdir)/average_nusselt_number.C \
          $(qoi_srcdir)/qoi_base.C \
          $(qoi_srcdir)/qoi_factory.C \
          $(qoi_srcdir)/vorticity.C

sfilesdir = $(prefix)/src
sfiles_DATA = $(physics_srcdir)/axisym_heat_transfer.C

#----------------------------
#INCLUDES we want distributed
#----------------------------
solver_INC   = $(solver_incdir)/grins_solver.h \
               $(solver_incdir)/mesh_builder.h \
               $(solver_incdir)/simulation.h \
               $(solver_incdir)/solver_factory.h \
               $(solver_incdir)/grins_steady_solver.h \
               $(solver_incdir)/grins_unsteady_solver.h \
               $(solver_incdir)/simulation_builder.h

physics_INC  = $(physics_incdir)/multiphysics_sys.h \
               $(physics_incdir)/physics.h \
               $(physics_incdir)/variable_name_defaults.h \
               $(physics_incdir)/var_typedefs.h \
               $(physics_incdir)/axisym_inc_navier_stokes.h \
               $(physics_incdir)/stokes.h \
               $(physics_incdir)/inc_navier_stokes_base.h \
               $(physics_incdir)/inc_navier_stokes.h \
               $(physics_incdir)/inc_navier_stokes_stab_base.h \
               $(physics_incdir)/inc_navier_stokes_adjoint_stab.h \
               $(physics_incdir)/heat_transfer_base.h \
               $(physics_incdir)/heat_transfer.h \
               $(physics_incdir)/heat_transfer_source.h \
               $(physics_incdir)/heat_transfer_stab_base.h \
               $(physics_incdir)/heat_transfer_adjoint_stab.h \
               $(physics_incdir)/axisym_heat_transfer.h \
               $(physics_incdir)/boussinesq_buoyancy.h \
               $(physics_incdir)/axisym_boussinesq_buoyancy.h \
               $(physics_incdir)/physics_factory.h \
               $(physics_incdir)/grins_physics_names.h \
               $(physics_incdir)/low_mach_navier_stokes_base.h \
               $(physics_incdir)/low_mach_navier_stokes.h \
               $(physics_incdir)/low_mach_navier_stokes_stab_base.h \
               $(physics_incdir)/low_mach_navier_stokes_spgsm_stab.h \
               $(physics_incdir)/low_mach_navier_stokes_vms_stab.h \
               $(physics_incdir)/low_mach_navier_stokes_braack_stab.h \
               $(physics_incdir)/stab_helper.h \
               $(physics_incdir)/inc_navier_stokes_stab_helper.h \
               $(physics_incdir)/heat_transfer_stab_helper.h \
               $(physics_incdir)/low_mach_navier_stokes_stab_helper.h

boundary_conditions_INC = $(boundary_conditions_incdir)/bc_types.h \
                          $(boundary_conditions_incdir)/boundary_conditions.h \
                          $(boundary_conditions_incdir)/bc_factory.h \
                          $(boundary_conditions_incdir)/neumann_func_obj.h \
                          $(boundary_conditions_incdir)/parabolic_profile.h \
                          $(boundary_conditions_incdir)/gaussian_profile.h \
                          $(boundary_conditions_incdir)/gaussian_xy_profile.h \
                          $(boundary_conditions_incdir)/concentric_cylinder_profile.h \
                          $(boundary_conditions_incdir)/dbc_container.h \
                          $(boundary_conditions_incdir)/pressure_pinning.h \
                          $(boundary_conditions_incdir)/pbc_container.h \
                          $(boundary_conditions_incdir)/nbc_container.h

utilities_INC = $(utilities_incdir)/grins_version.h

visualization_INC = $(visualization_incdir)/steady_visualization.h \
                    $(visualization_incdir)/unsteady_visualization.h \
                    $(visualization_incdir)/visualization.h \
                    $(visualization_incdir)/visualization_factory.h

matprops_INC = $(matprops_incdir)/grins_constant_function.h \
               $(matprops_incdir)/constant_conductivity.h   \
               $(matprops_incdir)/constant_viscosity.h      \
               $(matprops_incdir)/constant_specific_heat.h  \
               $(matprops_incdir)/constant_source_func.h

bc_handling_INC = $(bc_handling_incdir)/bc_handling_base.h \
                  $(bc_handling_incdir)/heat_transfer_bc_handling.h \
                  $(bc_handling_incdir)/axisym_heat_transfer_bc_handling.h \
                  $(bc_handling_incdir)/inc_navier_stokes_bc_handling.h \
                  $(bc_handling_incdir)/axisym_inc_navier_stokes_bc_handling.h \
                  $(bc_handling_incdir)/low_mach_navier_stokes_bc_handling.h

qoi_INC = $(qoi_incdir)/average_nusselt_number.h \
          $(qoi_incdir)/qoi_base.h \
          $(qoi_incdir)/qoi_factory.h \
          $(qoi_incdir)/qoi_names.h \
          $(qoi_incdir)/vorticity.h

include_HEADERS = $(solver_INC) $(physics_INC) $(boundary_conditions_INC) \
                  $(utilities_INC) $(visualization_INC) $(matprops_INC) \
                  $(bc_handling_INC) $(qoi_INC)

libgrins_la_SOURCES = $(solver_SRC) \
                      $(physics_SRC) \
                      $(boundary_conditions_SRC) \
                      $(utilities_SRC) \
                      $(visualization_SRC) \
                      $(matprops_SRC) \
                      $(bc_handling_SRC) \
                      $(qoi_SRC)

if LIBMESH_LIBTOOL
   libgrins_la_LIBADD = $(LIBMESH_LIBDIR)/libmesh.la
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

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(solver_incdir) -I$(physics_incdir) -I$(boundary_conditions_incdir) -I$(utilities_incdir) -I$(visualization_incdir) -I$(matprops_incdir) -I$(bc_handling_incdir) -I$(qoi_incdir)

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_INCLUDE) \
               $(GRVY_CFLAGS)

#-------------
# MASA support
#-------------

if MASA_ENABLED
  AM_CPPFLAGS           += $(MASA_CXXFLAGS)
  LIBS                  += $(MASA_LIBS)
endif

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
        $(qoi_SRC)
	echo 'updated source license headers' >$@

EXTRA_DIST = common

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
