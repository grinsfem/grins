dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(grins, 0.3.0, pbauman@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(grins_config.h.tmp)
dnl Prefixes everything in the temporary file with FINS_
AX_PREFIX_CONFIG_H(grins_config.h,GRINS,grins_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.11 no-define color-tests parallel-tests])
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
AX_PATH_LIBMESH_NEW(0.9.0,yes)

AC_ARG_ENABLE([libmesh-flags],
               AC_HELP_STRING([--enable-libmesh-flags],[enable use of libmesh cxxflags for CXXFLAGS]),
               USE_LIBMESH_CXX_FLAGS=1
               AC_DEFINE(USE_LIBMESH_CXX_FLAGS,1,[Define if libmesh cxxflags enabled]),[])


if test "$USE_LIBMESH_CXX_FLAGS" = "1"; then
   CXXFLAGS=`$LIBMESH_CONFIG --cxxflags`
else
 dnl TODO: WE NEED TO DISTINGUISH BETWEEN COMPILERS IN A BETTER WAY
   CXXFLAGS="-O0 -g -Wall -fno-unsafe-math-optimizations"
fi

dnl---------------------------------------------------------
dnl Add libMesh flags manually if it's not a libtool build
dnl---------------------------------------------------------
if test "x$LIBMESH_LIBTOOL" == "x1"; then
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

dnl----------------
dnl Check for boost
dnl----------------
BOOST_REQUIRE([1.39.0])


dnl-------------------------------------
dnl Check for GNU Linear Programming kit
dnl-------------------------------------
dnl###AX_PATH_GLPK(4.35,,AC_MSG_ERROR([Could not find required GLPK error]))
dnl AC_CACHE_SAVE

dnl-------------------------------
dnl Check for GRVY
dnl
dnl libGRVY Config 
dnl (min_version,required=yes/no)
dnl-------------------------------
AX_PATH_GRVY_NEW([0.30],[no])

dnl Check and see if the user wants to use GRVY timers
AC_ARG_ENABLE([grvy-timers],
              [AS_HELP_STRING([--enable-grvy-timers], [use GRVY timers (default is no)])],
	      [grvy_timers=${enableval}],
	      [grvy_timers=no])

dnl We'll rely on GRVY to make sure it's built correctly, but want to warn the user if
dnl they asked for timers, but GRVY wasn't found.
if test "x${grvy_timers}" = xyes; then
   if test "x$HAVE_GRVY" = x0; then
      AC_MSG_WARN([You asked for GRVY timers, but GRVY install not detected. Disabling timers. ])
   fi

   if test "x$HAVE_GRVY" = x1; then
      AC_DEFINE([USE_GRVY_TIMERS], 1, [GRVY timers enabled])
      USE_GRVY_TIMERS=1
      AC_SUBST(USE_GRVY_TIMERS)
   fi
fi


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
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(GRINS, doxygen/grins.dox, docs/doxygen)

AC_CONFIG_FILES(doxygen/txt_common/about.page)

dnl-----------------------------------------------
dnl Generate files for unit and regression testing
dnl-----------------------------------------------
AC_CONFIG_FILES(test/test_axi_ns_poiseuille_flow.sh,	      [chmod +x test/test_axi_ns_poiseuille_flow.sh])
AC_CONFIG_FILES(test/test_ns_couette_flow_2d_x.sh,            [chmod +x test/test_ns_couette_flow_2d_x.sh])
AC_CONFIG_FILES(test/test_ns_couette_flow_2d_y.sh,            [chmod +x test/test_ns_couette_flow_2d_y.sh])
AC_CONFIG_FILES(test/test_ns_poiseuille_flow.sh,              [chmod +x test/test_ns_poiseuille_flow.sh])
AC_CONFIG_FILES(test/test_stokes_poiseuille_flow.sh,          [chmod +x test/test_stokes_poiseuille_flow.sh])
AC_CONFIG_FILES(test/test_thermally_driven_2d_flow.sh,        [chmod +x test/test_thermally_driven_2d_flow.sh])
AC_CONFIG_FILES(test/test_thermally_driven_3d_flow.sh,        [chmod +x test/test_thermally_driven_3d_flow.sh])
AC_CONFIG_FILES(test/test_axi_thermally_driven_flow.sh,       [chmod +x test/test_axi_thermally_driven_flow.sh])
AC_CONFIG_FILES(test/test_axi_ns_con_cyl_flow.sh,             [chmod +x test/test_axi_ns_con_cyl_flow.sh])
AC_CONFIG_FILES(test/test_vorticity_qoi.sh,                   [chmod +x test/test_vorticity_qoi.sh])
AC_CONFIG_FILES(test/input_files/vorticity_qoi.in)
AC_CONFIG_FILES(test/low_mach_cavity_benchmark_regression.sh, [chmod +x test/low_mach_cavity_benchmark_regression.sh])
AC_CONFIG_FILES(test/backward_facing_step_regression.sh,      [chmod +x test/backward_facing_step_regression.sh])
AC_CONFIG_FILES(test/input_files/backward_facing_step.in)
AC_CONFIG_FILES(test/coupled_stokes_ns.sh,                    [chmod +x test/coupled_stokes_ns.sh])
AC_CONFIG_FILES(test/input_files/coupled_stokes_ns.in)
AC_CONFIG_FILES(test/magnetostatics_regression.sh,            [chmod +x test/magnetostatics_regression.sh])

dnl-----------------------------------------------
dnl Generate run scripts for examples
dnl-----------------------------------------------
AC_CONFIG_FILES(examples/cavity_benchmark/run.sh,	     [chmod +x examples/cavity_benchmark/run.sh])
AC_CONFIG_FILES(examples/backward_facing_step/run.sh,	     [chmod +x examples/backward_facing_step/run.sh])
AC_CONFIG_FILES(examples/lid_driven_cavity/run.sh,	     [chmod +x examples/lid_driven_cavity/run.sh])
AC_CONFIG_FILES(examples/mass_injection/run.sh,	             [chmod +x examples/mass_injection/run.sh])
AC_CONFIG_FILES(examples/rayleigh_taylor/run.sh,	     [chmod +x examples/rayleigh_taylor/run.sh])
AC_CONFIG_FILES(examples/coupled_stokes_navierstokes/run.sh, [chmod +x examples/coupled_stokes_navierstokes/run.sh])

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/grins/grins_version.h)
            
dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile doxygen/Makefile examples/Makefile])

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
