BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = grins grins_version
lib_LTLIBRARIES = libgrins.la

#--------------------------------------------------
# List of directories with source and include files
#--------------------------------------------------
grins_srcdir = $(top_srcdir)/src/grins_solver/src
grins_incdir = $(top_srcdir)/src/grins_solver/include

physics_srcdir = $(top_srcdir)/src/physics/src
physics_incdir = $(top_srcdir)/src/physics/include

boundary_conditions_srcdir = $(top_srcdir)/src/boundary_conditions/src
boundary_conditions_incdir = $(top_srcdir)/src/boundary_conditions/include

utilities_srcdir = $(top_srcdir)/src/utilities/src
utilities_incdir = $(top_srcdir)/src/utilities/include

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
# Sources in src/grins
grins_SRC    = $(grins_srcdir)/grins_solver.C \
               $(grins_srcdir)/grins_mesh_manager.C

physics_SRC  = $(physics_srcdir)/multiphysics_sys.C \
               $(physics_srcdir)/physics.C \
               $(physics_srcdir)/inc_navier_stokes.C \
               $(physics_srcdir)/axisym_inc_navier_stokes.C \
               $(physics_srcdir)/heat_transfer.C \
               $(physics_srcdir)/axisym_heat_transfer.C \
               $(physics_srcdir)/boussinesq_buoyancy.C \
               $(physics_srcdir)/axisym_boussinesq_buoyancy.C \
               $(physics_srcdir)/axisym_mushy_zone_solidification.C

boundary_conditions_SRC = $(boundary_conditions_srcdir)/boundary_conditions.C \
                          $(boundary_conditions_srcdir)/dirichlet_func_obj.C \
                          $(boundary_conditions_srcdir)/neumann_func_obj.C \
                          $(boundary_conditions_srcdir)/parabolic_profile.C \
                          $(boundary_conditions_srcdir)/concentric_cylinder_profile.C

utilities_SRC = $(utilities_srcdir)/grins_version.C

#----------------------------
#INCLUDES we want distributed
#----------------------------
grins_INC    = $(grins_incdir)/grins_solver.h \
               $(grins_incdir)/grins_mesh_manager.h

physics_INC  = $(physics_incdir)/multiphysics_sys.h \
               $(physics_incdir)/physics.h \
               $(physics_incdir)/variable_name_defaults.h \
               $(physics_incdir)/var_typedefs.h \
               $(physics_incdir)/axisym_inc_navier_stokes.h \
               $(physics_incdir)/inc_navier_stokes.h \
               $(physics_incdir)/heat_transfer.h \
               $(physics_incdir)/axisym_heat_transfer.h \
               $(physics_incdir)/boussinesq_buoyancy.h \
               $(physics_incdir)/axisym_boussinesq_buoyancy.h \
               $(physics_incdir)/axisym_mushy_zone_solidification.h

boundary_conditions_INC = $(boundary_conditions_incdir)/bc_types.h \
                          $(boundary_conditions_incdir)/boundary_conditions.h \
                          $(boundary_conditions_incdir)/dirichlet_func_obj.h \
                          $(boundary_conditions_incdir)/neumann_func_obj.h \
                          $(boundary_conditions_incdir)/parabolic_profile.h \
                          $(boundary_conditions_incdir)/concentric_cylinder_profile.h

utilities_INC = $(utilities_incdir)/grins_version.h

include_HEADERS = $(grins_INC) $(physics_INC) $(boundary_conditions_INC) \
                  $(utilities_INC)

libgrins_la_SOURCES = $(grins_SRC) \
                      $(physics_SRC) \
                      $(boundary_conditions_SRC) \
                      $(utilities_SRC)

#------------------------------------
# Source for stand alone GRINS driver
#------------------------------------

grins_SOURCES = $(grins_srcdir)/grins.C
grins_LDADD = libgrins.la

grins_version_SOURCES = $(utilities_srcdir)/version.C
grins_version_LDADD = libgrins.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(grins_incdir) -I$(physics_incdir) -I$(boundary_conditions_incdir) -I$(utilities_incdir)

AM_CPPFLAGS += $(LIBMESH_INCLUDE) \
               $(GRVY_CFLAGS)

AM_LDFLAGS = $(LIBMESH_LIBS) \
             $(GRVY_LIBS)

#-------------
# MASA support
#-------------

if MASA_ENABLED
  AM_CPPFLAGS           += $(MASA_CXXFLAGS)
  AM_LDFLAGS            += $(MASA_LIBS)
endif

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(grins_SRC) $(grins_INC) \
	$(generic_driver_SOURCES) $(physics_SRC) $(physics_INC) \
	$(boundary_conditions_INC) $(boundary_conditions_SRC) \
        $(utilities_INC) $(utilities_SRC) \
        $(grins_SOURCES) $(grins_version_SOURCES)
	echo 'updated source license headers' >$@

EXTRA_DIST = common

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
