check_PROGRAMS  = # Append below
check_PROGRAMS += test_libmesh_init

AM_CPPFLAGS = -I$(top_srcdir)/src/grins_solver/include
#	      \ -I$(top_srcdir)/

AM_CPPFLAGS +=  $(QUESO_CPPFLAGS) \
                $(BOOST_CPPFLAGS) \
                $(GRVY_CFLAGS) \
                $(GSL_CFLAGS) \
                $(GLPK_CFLAGS) \
                $(LIBMESH_CXXFLAGS) \
                $(HDF5_CFLAGS)  \
                $(EPETRA_INCLUDES)

AM_LDFLAGS = $(QUESO_LIBS) \
             $(GSL_LIBS) \
             $(GLPK_LIBS) \
             $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
             $(GRVY_LIBS) \
             $(HDF5_LIBS) \
             $(EPETRA_LIBS) \
             $(LIBMESH_LIBS) \
             -L$(top_builddir)/src -lgrins

#--------------
# MASA support
#--------------

if MASA_ENABLED
  AM_CPPFLAGS           += $(MASA_CXXFLAGS)
  AM_LDFLAGS            += $(MASA_LIBS)
endif

# Sources for these tests
test_libmesh_init_SOURCES = test_libmesh_init.C

#Define tests to actually be run
TESTS = test_libmesh_init

#shellfiles_src = <NULL.sh>

# Want these put with the distro so we can run make check
#EXTRA_DIST = $(shellfiles_src)

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
