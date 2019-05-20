 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(grins, 0.9.0, grins-users@googlegroups.com)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(grins_config.h.tmp)
dnl Prefixes everything in the temporary file with FINS_
AX_PREFIX_CONFIG_H(grins_config.h,GRINS,grins_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-pax 1.11 no-define color-tests parallel-tests subdir-objects])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

LT_INIT

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=grins

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])



dnl------------------------------
dnl Checks for compilers and MPI
dnl------------------------------

AC_PROG_CC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_PROG_CXX
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

dnl Added for AX_CODE_COVERAGE macro to work correctly, even though there's no Fortran here.
AC_PROG_FC

dnl----------------
dnl Libmesh Config
dnl----------------
AX_PATH_LIBMESH_NEW(1.0.0,yes)
AX_LIBMESH_DIM()

AC_ARG_ENABLE([libmesh-flags],
               AC_HELP_STRING([--disable-libmesh-flags],[disable use of libmesh cxxflags for CXXFLAGS]),
               [case "${enableval}" in
                   yes) USE_LIBMESH_CXX_FLAGS=1 ;;
                    no) USE_LIBMESH_CXX_FLAGS=0 ;;
 		     *) AC_MSG_ERROR(bad value ${enableval} for libmesh-flags) ;;
                esac],
               [USE_LIBMESH_CXX_FLAGS=1])

if test "$USE_LIBMESH_CXX_FLAGS" = "1"; then
   AC_DEFINE(USE_LIBMESH_CXX_FLAGS,1,[Define if libmesh cxxflags enabled])
   CXXFLAGS=`$LIBMESH_CONFIG --cppflags --cxxflags`
else
 dnl TODO: WE NEED TO DISTINGUISH BETWEEN COMPILERS IN A BETTER WAY
   AC_DEFINE(USE_LIBMESH_CXX_FLAGS,0,[Define if libmesh cxxflags enabled])
   CXXFLAGS="-O0 -g -fno-unsafe-math-optimizations"
fi

dnl -Wall warnings, -Wall the time.
AX_CXXFLAGS_WARN_ALL

dnl--------------------------
dnl Even if no dependency needed it,
dnl C++11 makes our asserts nicer,
dnl so let's try and enable it.
dnl--------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--enable-cxx11],
                             [build with C++11 support]),
              [case "${enableval}" in
                yes)  enablecxx11=yes ;;
                 no)  enablecxx11=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cxx11) ;;
               esac],
               [enablecxx11=optional])

if (test x$enablecxx11 = xyes); then
  AX_CXX_COMPILE_STDCXX_11(noext, mandatory)
elif (test x$enablecxx11 = xoptional); then
  AX_CXX_COMPILE_STDCXX_11(noext, optional)
fi

dnl---------------------------------------------------------
dnl Add libMesh flags manually if it's not a libtool build
dnl---------------------------------------------------------
if test "x$LIBMESH_HAVE_LIBTOOL" == "x1"; then
  LIBGRINS_LIBS='$(top_builddir)/src/libgrins.la'
else
  LIBGRINS_LIBS='$(top_builddir)/src/libgrins.la $(LIBMESH_LDFLAGS)'
fi
AC_SUBST(LIBGRINS_LIBS)

dnl---------------------------------------------------------
dnl Checks for  library prerequisites for other libraries...
dnl---------------------------------------------------------

dnl Optional check for MASA
AX_PATH_MASA(0.20,no)

dnl QUESO dependencies

dnl----------------
dnl Check for GSL
dnl----------------
dnl### AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))
dnl AC_CACHE_SAVE


dnl------------------------------------------
dnl Check if libMesh detected std::shared_ptr
dnl------------------------------------------
AC_MSG_CHECKING([for C++11 std::shared_ptr support from libMesh])
grins_save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$LIBMESH_CPPFLAGS $CPPFLAGS"
libmesh_have_cxx11_shared_ptr=no
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
                                      @%:@include <libmesh/libmesh_config.h>
                                   ]],
                                   [[
                                      #if LIBMESH_HAVE_CXX11_SHARED_PTR
                                      #else
                                      #error "libMesh shared_ptr support not found!"
                                      #endif
                                   ]])
                  ],
                                   [
                libmesh_have_cxx11_shared_ptr=yes
                AC_MSG_RESULT(yes)
            ],[
                libmesh_have_cxx11_shared_ptr=no
                AC_MSG_RESULT(no)
])
CPPFLAGS=$grins_save_CPPFLAGS


dnl-----------------------
dnl Check for Cantera 2.0+
dnl-----------------------
AC_ARG_VAR([CANTERA_DIR],[Location of Cantera installation])
AC_ARG_ENABLE(cantera,
  [  --enable-cantera        Compile with Cantera chemistry support],
       enable_cantera=$enableval,
       enable_cantera=yes)

HAVE_CANTERA=0

if test "$enable_cantera" != no; then
  # User-specific include path
  AC_ARG_WITH(cantera,
               AC_HELP_STRING([--with-cantera=PATH],[Specify the path for Cantera]),
               with_cantera=$withval,
               with_cantera=$CANTERA_DIR)


   ac_cantera_save_CPPFLAGS="$CPPFLAGS"
   ac_cantera_save_LDFLAGS="$LDFLAGS"
   ac_cantera_save_LIBS="$LIBS"

   CANTERA_PREFIX=$with_cantera

   CPPFLAGS="-I$CANTERA_PREFIX/include"
   # v2.0.x with scons
   LDFLAGS="-L$CANTERA_PREFIX/lib"
   LIBS="-lcantera_shared"

   if (test -f $CANTERA_PREFIX/lib/libcvode_shared*); then
      LIBS="$LIBS -lcvode_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libctmath_shared*); then
      LIBS="$LIBS -lctmath_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libctf2c_shared*); then
      LIBS="$LIBS -lctf2c_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libexecstream_shared*); then
      LIBS="$LIBS -lexecstream_shared"
   fi

   CANTERA_CPPFLAGS="${CPPFLAGS}"
   CANTERA_LDFLAGS="${LDFLAGS}"
   CANTERA_LIBS="${LIBS}"

   dnl We don't want these in the final link line since we'll already have BLAS
   dnl from PETSC etc., but these are needed for the test to pass.
   if (test -f $CANTERA_PREFIX/lib/libctblas_shared*); then
      LIBS="$LIBS -lctblas_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libctlapack_shared*); then
      LIBS="$LIBS -lctlapack_shared"
   fi

   #--------------------------------------------------------------
   # Check for Cantera library linking
   #--------------------------------------------------------------
   AC_MSG_CHECKING([for Cantera linkage])

   AC_LANG_PUSH([C++])
   AC_LINK_IFELSE( [AC_LANG_PROGRAM([#include "cantera/IdealGasMix.h"],
                                    [Cantera::IdealGasMix()])],
                                    [AC_MSG_RESULT(yes)
                                     found_cantera_library=yes],
                                    [AC_MSG_RESULT(no)
                                     found_cantera_library=no] )
   AC_LANG_POP([C++])

   CPPFLAGS="$ac_cantera_save_CPPFLAGS"
   LDFLAGS="$ac_cantera_save_LDFLAGS"
   LIBS="$ac_cantera_save_LIBS"

   if test "x${found_cantera_library}" = "xyes" ; then
      HAVE_CANTERA=1
      AC_DEFINE(HAVE_CANTERA, 1, [Flag indicating support for Cantera chemistry])
      AC_SUBST(CANTERA_CPPFLAGS)
      AC_SUBST(CANTERA_LDFLAGS)
      AC_SUBST(CANTERA_LIBS)
      AC_SUBST(CANTERA_PREFIX)
   else
      AC_MSG_NOTICE([Disabling optional Cantera chemistry support])
      CANTERA_CPPFLAGS=""
      CANTERA_LDFLAGS=""
      CANTERA_LIBS=""
   fi
fi

AC_SUBST(HAVE_CANTERA)
AM_CONDITIONAL(CANTERA_ENABLED,test x$HAVE_CANTERA = x1)


AX_PATH_ANTIOCH(0.4.0,no)

# -------------------------------------------------------------
# cppunit C++ unit testing -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(cppunit,
              AS_HELP_STRING([--disable-cppunit],
                             [Build without cppunit C++ unit testing support]))
AS_IF([test "x$enable_cppunit" != "xno"], [
   AM_PATH_CPPUNIT([1.10.0],[enablecppunit=yes],[enablecppunit=no])
])

AM_CONDITIONAL(GRINS_ENABLE_CPPUNIT, test x$enablecppunit = xyes)

dnl-------------------
dnl Check for TRILINOS
dnl-------------------
dnl###AX_TRILINOS_BASE(8.0.7,,AC_MSG_ERROR([Could not find required Trilinos version.]))
dnl###AX_TRILINOS_EPETRA(AC_DEFINE(HAVE_LIBEPETRA,1,[Define if Epetra is available]),AC_MSG_ERROR([Could not find Epetra.]))
dnl AC_CACHE_SAVE


dnl---------------------------------------------------------------
dnl Checks for libraries
dnl---------------------------------------------------------------

dnl QUESO REQUIRES HDF5
dnl###AX_PATH_HDF5


dnl----------------
dnl Check for QUESO
dnl----------------
dnl###AX_PATH_QUESO
dnl AC_CACHE_SAVE


dnl--------------------------
dnl Checks for code coverage
dnl--------------------------
AX_CODE_COVERAGE


dnl---------------------------------
dnl Query configuration environment
dnl Must be done *before* Doxygen to
dnl get some happy variables defined.
dnl---------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------------------
dnl            Doxygen support
dnl------------------------------------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(GRINS, doxygen/grins.dox, docs/doxygen)

AC_CONFIG_FILES(doxygen/txt_common/about_vpath.page)

dnl-------------------------------------------------------------------------
dnl Generate symlink to helper script for running tests through sequence of
dnl LIBMESH_RUN settings.
dnl-------------------------------------------------------------------------
AC_CONFIG_LINKS([test/run_tests_parallel_loop.sh:test/run_tests_parallel_loop.sh],[chmod +x test/run_tests_parallel_loop.sh])

dnl-------------------------------------------------------------------------
dnl Generate symlinks to allow tests to read grids without having to require
dnl AC_CONFIG_FILES to generate the input files.
dnl NOTE: AC_CONFIG_LINKS will generate the rules for `make dist`
dnl-------------------------------------------------------------------------
AC_CONFIG_LINKS([test/grids/backward_facing_step.e:test/grids/backward_facing_step.e])
AC_CONFIG_LINKS([test/grids/circular_sheet.e:test/grids/circular_sheet.e])
AC_CONFIG_LINKS([test/grids/coupled_stokes_ns.e:test/grids/coupled_stokes_ns.e])
AC_CONFIG_LINKS([test/grids/mixed_dim_circle.exo:test/grids/mixed_dim_circle.exo])
AC_CONFIG_LINKS([test/grids/mixed_dim_symmetric_square.exo:test/grids/mixed_dim_symmetric_square.exo])
AC_CONFIG_LINKS([test/grids/poiseuille_verification.e:test/grids/poiseuille_verification.e])
AC_CONFIG_LINKS([test/grids/mixed_quad_tri_square_mesh.xda:test/grids/mixed_quad_tri_square_mesh.xda])
AC_CONFIG_LINKS([test/grids/cylinder.exo:test/grids/cylinder.exo])
AC_CONFIG_LINKS([test/grids/one_elem_overlap.exo:test/grids/one_elem_overlap.exo])

dnl-------------------------------------------------------------------------
dnl Generate symlinks to allow tests to read restarts without having to
dnl require AC_CONFIG_FILES to generate the .sh files
dnl NOTE: AC_CONFIG_LINKS will generate the rules for `make dist`
dnl-------------------------------------------------------------------------
AC_CONFIG_LINKS([test/test_data/elastic_mooney_rivlin_square_hookean_stiffeners_regression.xdr:test/test_data/elastic_mooney_rivlin_square_hookean_stiffeners_regression.xdr])
AC_CONFIG_LINKS([test/test_data/CO2_data.dat:test/test_data/CO2_data.dat])
AC_CONFIG_LINKS([test/test_data/CO2_partition_function.dat:test/test_data/CO2_partition_function.dat])
AC_CONFIG_LINKS([test/test_data/CO2_N2.xml:test/test_data/CO2_N2.xml])

dnl------------------------------------------------------------------------------
dnl Generate symlinks to allow tests to read these files from the test input file
dnl NOTE: AC_CONFIG_LINKS will generate the rules for `make dist`
dnl------------------------------------------------------------------------------
AC_CONFIG_LINKS([test/input_files/air_2sp.xml:test/input_files/air_2sp.xml])
AC_CONFIG_LINKS([test/input_files/air_4sp.xml:test/input_files/air_4sp.xml])
AC_CONFIG_LINKS([test/input_files/air_5sp_test.xml:test/input_files/air_5sp_test.xml])
AC_CONFIG_LINKS([test/input_files/air.xml:test/input_files/air.xml])
AC_CONFIG_LINKS([test/input_files/elements.xml:test/input_files/elements.xml])
AC_CONFIG_LINKS([test/input_files/ozone_flame_antioch_constant_restart.xdr:test/input_files/ozone_flame_antioch_constant_restart.xdr])
AC_CONFIG_LINKS([test/input_files/ozone_flame_cantera_restart.xdr:test/input_files/ozone_flame_cantera_restart.xdr])
AC_CONFIG_LINKS([test/input_files/ozone.xml:test/input_files/ozone.xml])
AC_CONFIG_LINKS([test/input_files/ozone_cea_data.dat:test/input_files/ozone_cea_data.dat])
AC_CONFIG_LINKS([test/input_files/ozone_species_data.dat:test/input_files/ozone_species_data.dat])
AC_CONFIG_LINKS([test/input_files/ozone_transport_data.dat:test/input_files/ozone_transport_data.dat])

dnl----------------------------------------------
dnl Generate header with paths for use in testing
dnl----------------------------------------------
AC_CONFIG_FILES(test/common/grins_test_paths.h)

dnl-----------------------------------------------
dnl Generate run scripts for examples
dnl-----------------------------------------------
AC_CONFIG_FILES(examples/cavity_benchmark/run.sh,	     [chmod +x examples/cavity_benchmark/run.sh])
AC_CONFIG_FILES(examples/backward_facing_step/run.sh,	     [chmod +x examples/backward_facing_step/run.sh])
AC_CONFIG_FILES(examples/lid_driven_cavity/run.sh,	         [chmod +x examples/lid_driven_cavity/run.sh])
AC_CONFIG_FILES(examples/velocity_penalty/run.sh,	         [chmod +x examples/velocity_penalty/run.sh])
AC_CONFIG_FILES(examples/convection_cell/run.sh,	         [chmod +x examples/convection_cell/run.sh])
AC_CONFIG_FILES(examples/rayleigh_taylor/run.sh,	         [chmod +x examples/rayleigh_taylor/run.sh])
AC_CONFIG_FILES(examples/simple_fan/run.sh,	                 [chmod +x examples/simple_fan/run.sh])
AC_CONFIG_FILES(examples/simple_prop/run.sh,	             [chmod +x examples/simple_prop/run.sh])
AC_CONFIG_FILES(examples/coupled_stokes_navierstokes/run.sh, [chmod +x examples/coupled_stokes_navierstokes/run.sh])
AC_CONFIG_FILES(examples/elastic_sheet/run.sh,               [chmod +x examples/elastic_sheet/run.sh])
AC_CONFIG_FILES(examples/cable_displacement/run.sh,          [chmod +x examples/cable_displacement/run.sh])
AC_CONFIG_FILES(examples/inflating_sheet/run.sh,             [chmod +x examples/inflating_sheet/run.sh])
AC_CONFIG_FILES(examples/sa_airfoil/run.sh,         [chmod +x examples/sa_airfoil/run.sh])
AC_CONFIG_FILES(examples/suspended_cable/run.sh,             [chmod +x examples/suspended_cable/run.sh])
AC_CONFIG_FILES(examples/mixed_dim_inflating_sheet/run.sh,   [chmod +x examples/mixed_dim_inflating_sheet/run.sh])
AC_CONFIG_FILES(examples/transient_amr/run.sh,   [chmod +x examples/transient_amr/run.sh])
AC_CONFIG_FILES(examples/ozone_flame/run.sh,                 [chmod +x examples/ozone_flame/run.sh])
AC_CONFIG_FILES(examples/ozone_flame/crun.sh,                 [chmod +x examples/ozone_flame/crun.sh])
AC_CONFIG_FILES(examples/spectroscopy/run.sh,                 [chmod +x examples/spectroscopy/run.sh])
AC_CONFIG_FILES(examples/multigrid/run_poisson.sh,            [chmod +x examples/multigrid/run_poisson.sh])
AC_CONFIG_FILES(examples/multigrid/run_stokes.sh,             [chmod +x examples/multigrid/run_stokes.sh])
AC_CONFIG_FILES(examples/multigrid/run_thermo_coupled_flow.sh, [chmod +x examples/multigrid/run_thermo_coupled_flow.sh])
AC_CONFIG_FILES(examples/torsion_hyperelasticity/run.sh,     [chmod +x examples/torsion_hyperelasticity/run.sh])

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/grins/grins_version.h)

dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile
                 doxygen/Makefile
                 doxygen/grins.dox
                 examples/Makefile])


dnl-----------------------------------------------
dnl Collect files for "header stamping" here
dnl-----------------------------------------------

# PB: Although the documentation says that ac_top_srcdir should be defined,
#     it doesn't get activated until _AC_SRCDIRS gets called
#     (used to be _AC_SRCPATHS), according to this thread:
#     http://lists.gnu.org/archive/html/autoconf-patches/2003-02/msg00003.html
#     My own hacking revealed that we could also do the following:
#     TOP_SEARCH_DIR=$ac_pwd/$srcdir
#     I'm not sure which will be more "future proof", but this is the alternative
#     to using $(shell) calls in the Makefile.am which 1. requires
#     GNU Make and 2. Causes automake to issue a warning
_AC_SRCDIRS(.)
ABS_TOP_SRC_DIR=$ac_abs_top_srcdir
TOP_SRC_DIR=$ac_srcdir

# Since we don't distribute the lic_utils, check and make sure it's there.
# This way, we won't run the license stamping perl script on distributed
# tarballs, only on the repository clones
AM_CONDITIONAL(GRINS_LICENSESTAMPEXISTS, [test -f $ABS_TOP_SRC_DIR/src/common/lic_utils/update_license.pl])

# We have to do this by subdirectory because otherwise distcheck
# breaks as we start picking up files in the directories
# that it uses.
GRINS_STAMPED_FILES=$(find $ABS_TOP_SRC_DIR/{src,examples,test} -name "*.h" -or -name "*.C" | tr "\n" " ")
AC_SUBST(GRINS_STAMPED_FILES)

# PB: This *must not* be an absolute path since we'll use this variable
#     in EXTRA_DIST. `make dist` always checks './'.
GRINS_TEST_SH_FILES=$(find $TOP_SRC_DIR/test -name "*.sh" | tr "\n" " ")
AC_SUBST(GRINS_TEST_SH_FILES)

# PB: Since the .xml input files are put into `make dist` by AC_CONFIG_LINKS,
#     we need to do the others separately because if add the directory as a whole
#     to EXTRA_DIST in a Makefile.am, somehow the AC_CONFIG_LINKS parts don't
#     end up in the dist tarball...
GRINS_TEST_INPUT_FILES=$(find $TOP_SRC_DIR/test/input_files -name "*.in" | tr "\n" " ")
AC_SUBST(GRINS_TEST_INPUT_FILES)

# PB: Since there are some .xdr files put into `make dist` by AC_CONFIG_LINKS,
#     we need to do the others separately because if add the directory as a whole
#     to EXTRA_DIST in a Makefile.am, somehow the AC_CONFIG_LINKS parts don't
#     end up in the dist tarball...
GRINS_TEST_DATA_FILES=$(find $TOP_SRC_DIR/test/test_data -name "*.xdr" -or -name "*.xda" -or -name "*.gz" | grep -v "elastic_mooney_rivlin_square_hookean_stiffeners_regression.xdr" | tr "\n" " ")
AC_SUBST(GRINS_TEST_DATA_FILES)

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl Hackish fix for Ubuntu/Debian libtool
perl -pi -e 's/link_all_deplibs=no/link_all_deplibs=yes/' libtool

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
