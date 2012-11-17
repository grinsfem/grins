 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(grins, 0.3.0, pbauman@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.10 no-define])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)]) # use silent rules where available - automake 1.11

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
AX_PATH_LIBMESH_NEW(0.8.0,yes)

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

   CANTERA_PREFIX=$with_cantera

   CPPFLAGS="-I$CANTERA_PREFIX/include"
   # v2.0.x with scons
   LDFLAGS="-L$CANTERA_PREFIX/lib -lcantera_shared"
   
   if (test -f $CANTERA_PREFIX/lib/libcvode_shared*); then
      LDFLAGS="$LDFLAGS -lcvode_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libctmath_shared*); then
      LDFLAGS="$LDFLAGS -lctmath_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libctf2c_shared*); then
      LDFLAGS="$LDFLAGS -lctf2c_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libexecstream_shared*); then
      LDFLAGS="$LDFLAGS -lexecstream_shared"
   fi

   CANTERA_CPPFLAGS="${CPPFLAGS}"
   CANTERA_LDFLAGS="${LDFLAGS}"

   dnl We don't want these in the final link line since we'll already have BLAS
   dnl from PETSC etc., but these are needed for the test to pass.
   if (test -f $CANTERA_PREFIX/lib/libctblas_shared*); then
      LDFLAGS="$LDFLAGS -lctblas_shared"
   fi
   if (test -f $CANTERA_PREFIX/lib/libctlapack_shared*); then
      LDFLAGS="$LDFLAGS -lctlapack_shared"
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

   if test "x${found_cantera_library}" = "xyes" ; then
      HAVE_CANTERA=1
      AC_DEFINE(HAVE_CANTERA, 1, [Flag indicating support for Cantera chemistry])
      AC_SUBST(CANTERA_CPPFLAGS)
      AC_SUBST(CANTERA_LDFLAGS)
      AC_SUBST(CANTERA_PREFIX)
   else
      AC_MSG_NOTICE([Disabling optional Cantera chemistry support])
      CANTERA_CPPFLAGS=""
      CANTERA_LDFLAGS=""
   fi
fi

AC_SUBST(HAVE_CANTERA)
AM_CONDITIONAL(CANTERA_ENABLED,test x$HAVE_CANTERA = x1)

 
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

DX_INIT_DOXYGEN(GRINS, doxygen/grins.dox, docs/doxygen )


dnl-------------------
dnl Generate makefiles
dnl-------------------
AC_OUTPUT( Makefile \
           src/Makefile src/utilities/include/grins_version.h \
           test/Makefile \
           doxygen/Makefile \
           doxygen/txt_common/about.page \
           examples/Makefile \
           examples/cavity_benchmark/Makefile \
	   examples/rayleigh_taylor/Makefile \
           examples/mass_injection/Makefile \
           examples/vortex/Makefile \
           examples/rlmns/Makefile )



dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
