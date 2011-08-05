BUILT_SOURCES   = .license.stamp

check_PROGRAMS  = test_ns_couette_flow
check_PROGRAMS += test_ns_poiseuille_flow

AM_CPPFLAGS = -I$(top_srcdir)/src/grins_solver/include \
              -I$(top_srcdir)/src/physics/include \
              -I$(top_srcdir)/src/boundary_conditions/include \
              -I$(top_srcdir)/src/funcs/include

#	      \ -I$(top_srcdir)/src/utilities/<NULL>/include

###AM_CPPFLAGS +=  $(QUESO_CPPFLAGS) \
###                $(BOOST_CPPFLAGS) \
###                $(GRVY_CFLAGS) \
###                $(GSL_CFLAGS) \
###                $(GLPK_CFLAGS) \
###                $(LIBMESH_CXXFLAGS) \
###                $(HDF5_CFLAGS)  \
###                $(EPETRA_INCLUDES)

AM_CPPFLAGS += $(LIBMESH_CXXFLAGS) \
                $(GRVY_CFLAGS)

###AM_LDFLAGS = $(QUESO_LIBS) \
###             $(GSL_LIBS) \
###             $(GLPK_LIBS) \
###             $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
###             $(GRVY_LIBS) \
###             $(HDF5_LIBS) \
###             $(EPETRA_LIBS) \
###             $(LIBMESH_LIBS) \
###             -L$(top_builddir)/src -lgrins

AM_LDFLAGS = $(LIBMESH_LIBS) \
             $(GRVY_LIBS) \
             -L$(top_builddir)/src -lgrins

#--------------
# MASA support
#--------------

if MASA_ENABLED
  AM_CPPFLAGS           += $(MASA_CXXFLAGS)
  AM_LDFLAGS            += $(MASA_LIBS)
endif

# Sources for these tests
test_ns_couette_flow_SOURCES = test_ns_couette_flow.C
test_ns_poiseuille_flow_SOURCES = test_ns_poiseuille_flow.C

# List of source files for license stamping
testsrc = $(test_ns_couette_flow_SOURCES)            \
          $(test_ns_poiseuille_flow_SOURCES)

#Define tests to actually be run
TESTS = navier_stokes_tests.sh             \
        test_ns_couette_flow.sh            \
        test_ns_poiseuille_flow.sh \
        finalize.sh

#shellfiles_src = <NULL.sh>

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src)

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(testsrc)
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
