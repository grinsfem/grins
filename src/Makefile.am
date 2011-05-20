BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = grins
lib_LTLIBRARIES = libgrins.la

#--------------------------------------------------
# List of directories with source and include files
#--------------------------------------------------
grins_srcdir = $(top_srcdir)/src/grins_solver/src
grins_incdir = $(top_srcdir)/src/grins_solver/include

physics_srcdir = $(top_srcdir)/src/physics/src
physics_incdir = $(top_srcdir)/src/physics/include

# utilities may have more sub-dirs.
utils_dir    = $(top_srcdir)/src/utilities

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
# Sources in src/grins
grins_SRC    = $(grins_srcdir)/grins_solver.C \
               $(grins_srcdir)/grins_mesh_manager.C \
               $(grins_srcdir)/multiphysics_sys.C

physics_SRC  = $(physics_srcdir)/physics.C \
               $(physics_srcdir)/inc_navier_stokes.C

#Sources in utilities/src
#utils_SRC = $(utils_dir)/<NULL>/src/<NULL.C>

#----------------------------
#INCLUDES we want distributed
#----------------------------
grins_INC    = $(grins_incdir)/grins_solver.h \
               $(grins_incdir)/grins_mesh_manager.h \
               $(grins_incdir)/multiphysics_sys.h

physics_INC  = $(physics_incdir)/physics.h \
               $(physics_incdir)/inc_navier_stokes.h

#utils_INC = $(utils_dir)/<NULL>/include/<NULL.h>

include_HEADERS = $(grins_INC) $(utils_INC) $(physics_INC)

libgrins_la_SOURCES = $(grins_SRC) \
                      $(utils_SRC) \
                      $(physics_SRC)

#------------------------------------
# Source for stand alone GRINS driver
#------------------------------------

grins_SOURCES = $(grins_srcdir)/grins.C

grins_LDADD = libgrins.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(grins_incdir) -I$(physics_incdir)
#              \ -I$(utils_dir)/<NULL>/inc

#-------------------------------------
# External library include directories
#-------------------------------------
###AM_CPPFLAGS += $(HDF5_CFLAGS) \
###               $(QUESO_CPPFLAGS) \
###               $(BOOST_CPPFLAGS) \
###               $(GSL_CFLAGS) \
###               $(GLPK_CFLAGS) \
###               $(LIBMESH_CXXFLAGS) \
###               $(EPETRA_INCLUDES)

AM_CPPFLAGS += $(LIBMESH_INCLUDE) \
               $(GRVY_CFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
###AM_LDFLAGS = $(QUESO_LIBS) \
###             $(GSL_LIBS) \
###             $(GLPK_LIBS) \
###             $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
###             $(HDF5_LIBS)  \
###             $(EPETRA_LIBS) \
###             $(LIBMESH_LIBS)

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
	$(generic_driver_SOURCES) $(physics_SRC) $(physics_INC)
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
