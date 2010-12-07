dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(grins, 0.0.0, pbauman@ices.utexas.edu)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE(grins,0.0.0)

dnl Release versioning

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=0
GENERIC_MICRO_VERSION=0
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

AC_PROG_LIBTOOL

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
AX_PATH_LIBMESH
AC_CHECK_FILE($LIBMESH_PREFIX/contrib/bin/libmesh-config,
              LIBMESH_CONFIG=$LIBMESH_PREFIX/contrib/bin/libmesh-config,
              AC_MSG_ERROR([Cannot find libmesh-config!]))

AC_ARG_ENABLE([libmesh-flags],
               AC_HELP_STRING([--enable-libmesh-flags],[enable use of libmesh cxxflags for CXXFLAGS]),
               USE_LIBMESH_CXX_FLAGS=1
               AC_DEFINE(USE_LIBMESH_CXX_FLAGS,1,[Define if libmesh cxxflags enabled]),[])


if test "$USE_LIBMESH_CXX_FLAGS" = "1"; then
   CXXFLAGS=`$LIBMESH_CONFIG --cxxflags`
else
   CXXFLAGS="-O0 -g -fno-unsafe-math-optimizations"
fi

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
dnl###BOOST_REQUIRE([1.36])
dnl###BOOST_MATH()
dnl###BOOST_PROGRAM_OPTIONS()
dnl AC_CACHE_SAVE

dnl-------------------------------------
dnl Check for GNU Linear Programming kit
dnl-------------------------------------
dnl###AX_PATH_GLPK(4.35,,AC_MSG_ERROR([Could not find required GLPK error]))
dnl AC_CACHE_SAVE

dnl---------------
dnl Check for HPCT
dnl---------------
dnl###AX_PATH_GRVY
dnl AC_CACHE_SAVE

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

dnl---------------------------------------------------------------
dnl Checks for typedefs, structures, and compiler characteristics.
dnl---------------------------------------------------------------


dnl------------------------------
dnl Checks for library functions.
dnl------------------------------


dnl--------------------------
dnl Checks for code coverage
dnl--------------------------
AX_CODE_COVERAGE


dnl-------------------
dnl Generate makefiles
dnl-------------------
AC_OUTPUT( Makefile \
           src/Makefile \
           test/Makefile )

dnl--------------------------------
dnl Query configuration environment
dnl------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
