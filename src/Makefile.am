#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS =  grins_driver

lib_LTLIBRARIES = libgrins.la

#--------------------------------------------------
# List of directories with source and include files
#--------------------------------------------------
grins_srcdir = $(top_srcdir)/src/grins_solver/src
grins_incdir = $(top_srcdir)/src/grins_solver/include

# utilities may have more sub-dirs.
utils_dir = $(top_srcdir)/src/utilities

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
# Sources in src/grins
grins_SRC = $(grins_srcdir)/GRINS_Solver.C

#Sources in utilities/src
#utils_SRC = $(utils_dir)/<NULL>/src/<NULL.C>

#----------------------------
#INCLUDES we want distributed
#----------------------------
grins_INC = $(grins_incdir)/GRINS_Solver.h

#utils_INC = $(utils_dir)/<NULL>/include/<NULL.h>

include_HEADERS = $(grins_INC) $(utils_INC)

libgrins_la_SOURCES = $(grins_SRC) \
                      $(utils_SRC)

#------------------------------------
# Source for stand alone GRINS driver
#------------------------------------
grins_driver_SOURCES = $(grins_srcdir)/GRINS_Driver.C

grins_driver_LDADD = libgrins.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(grins_incdir)
#              \ -I$(utils_dir)/<NULL>/inc

#-------------------------------------
# External library include directories
#-------------------------------------
###AM_CPPFLAGS += $(HDF5_CFLAGS) \
###               $(QUESO_CPPFLAGS) \
###               $(BOOST_CPPFLAGS) \
###               $(GRVY_CFLAGS) \
###               $(GSL_CFLAGS) \
###               $(GLPK_CFLAGS) \
###               $(LIBMESH_CXXFLAGS) \
###               $(EPETRA_INCLUDES)

AM_CPPFLAGS += $(LIBMESH_CXXFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
###AM_LDFLAGS = $(QUESO_LIBS) \
###             $(GSL_LIBS) \
###             $(GLPK_LIBS) \
###             $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
###             $(GRVY_LIBS) \
###             $(HDF5_LIBS)  \
###             $(EPETRA_LIBS) \
###             $(LIBMESH_LIBS)

AM_LDFLAGS = $(LIBMESH_LIBS)

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

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
