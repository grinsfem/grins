#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(grins, 0.0.0, pbauman@ices.utexas.edu)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE(grins,0.0.0)

# Release versioning

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=0
GENERIC_MICRO_VERSION=0
GENERIC_LIBRARY_NAME=grins

#shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

# AC_SUBST(GENERIC_LIBRARY_VERSION)
# AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])

AC_PROG_LIBTOOL

#------------------------------
# Checks for compilers and MPI
#------------------------------

AC_PROG_CC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_PROG_CXX
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

# Added for AX_CODE_COVERAGE macro to work correctly, even though there's no Fortran here.
AC_PROG_FC

#----------------
# Libmesh Config
#----------------
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

#---------------------------------------------------------
# Checks for  library prerequisites for other libraries...
#---------------------------------------------------------

# Optional check for MASA
AX_PATH_MASA(0.20,no)

# QUESO dependencies

#----------------
# Check for GSL
#----------------
AX_PATH_GSL(1.10,AM_CONDITIONAL([UQBT_GSL], [test 'TRUE']),AC_MSG_ERROR([Could not find required GSL version.]))
#AC_CACHE_SAVE

#----------------
# Check for boost
#----------------
BOOST_REQUIRE([1.36])
BOOST_MATH()
BOOST_PROGRAM_OPTIONS()
#AC_CACHE_SAVE

#-------------------------------------
# Check for GNU Linear Programming kit
#-------------------------------------
AX_PATH_GLPK(4.35,,AC_MSG_ERROR([Could not find required GLPK error]))
#AC_CACHE_SAVE

#---------------
# Check for HPCT
#---------------
AX_PATH_GRVY
#AC_CACHE_SAVE

#-------------------
# Check for TRILINOS
#-------------------
AX_TRILINOS_BASE(8.0.7,,AC_MSG_ERROR([Could not find required Trilinos version.]))
AX_TRILINOS_EPETRA(AC_DEFINE(HAVE_LIBEPETRA,1,[Define if Epetra is available]),AC_MSG_ERROR([Could not find Epetra.]))
#AC_CACHE_SAVE


#---------------------------------------------------------------
# Checks for libraries
#---------------------------------------------------------------

# QUESO REQUIRES HDF5
AX_PATH_HDF5


#----------------
# Check for QUESO
#----------------
AX_PATH_QUESO
#AC_CACHE_SAVE

#---------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#---------------------------------------------------------------


#------------------------------
# Checks for library functions.
#------------------------------


#--------------------------
# Checks for code coverage
#--------------------------
AX_CODE_COVERAGE


#-------------------
# Generate makefiles
#-------------------
AC_OUTPUT( Makefile \
           src/Makefile \
           test/Makefile )

#--------------------------------
# Query configuration environment
#------------------------------
AX_SUMMARIZE_ENV

#------------------------------
# Final summary
#------------------------------
AX_SUMMARIZE_CONFIG
