 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(grins, 0.4.0, pbauman@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(grins_config.h.tmp)
dnl Prefixes everything in the temporary file with FINS_
AX_PREFIX_CONFIG_H(grins_config.h,GRINS,grins_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-pax 1.11 no-define color-tests parallel-tests])
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


AX_PATH_ANTIOCH(0.0.8,no)

 
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

dnl-----------------------------------------------
dnl Generate files for unit and regression testing
dnl-----------------------------------------------
AC_CONFIG_FILES(test/test_axi_ns_poiseuille_flow.sh,	                  [chmod +x test/test_axi_ns_poiseuille_flow.sh])
AC_CONFIG_FILES(test/test_ns_couette_flow_2d_x.sh,                        [chmod +x test/test_ns_couette_flow_2d_x.sh])
AC_CONFIG_FILES(test/test_ns_couette_flow_2d_y.sh,                        [chmod +x test/test_ns_couette_flow_2d_y.sh])
AC_CONFIG_FILES(test/test_ns_poiseuille_flow.sh,                          [chmod +x test/test_ns_poiseuille_flow.sh])
AC_CONFIG_FILES(test/test_stokes_poiseuille_flow.sh,                      [chmod +x test/test_stokes_poiseuille_flow.sh])
AC_CONFIG_FILES(test/test_thermally_driven_2d_flow.sh,                    [chmod +x test/test_thermally_driven_2d_flow.sh])
AC_CONFIG_FILES(test/test_thermally_driven_3d_flow.sh,                    [chmod +x test/test_thermally_driven_3d_flow.sh])
AC_CONFIG_FILES(test/test_dirichlet_nan.sh,                               [chmod +x test/test_dirichlet_nan.sh])
AC_CONFIG_FILES(test/test_axi_thermally_driven_flow.sh,                   [chmod +x test/test_axi_thermally_driven_flow.sh])
AC_CONFIG_FILES(test/test_axi_ns_con_cyl_flow.sh,                         [chmod +x test/test_axi_ns_con_cyl_flow.sh])
AC_CONFIG_FILES(test/test_vorticity_qoi.sh,                               [chmod +x test/test_vorticity_qoi.sh])
AC_CONFIG_FILES(test/input_files/vorticity_qoi.in)
AC_CONFIG_FILES(test/low_mach_cavity_benchmark_regression.sh,             [chmod +x test/low_mach_cavity_benchmark_regression.sh])
AC_CONFIG_FILES(test/backward_facing_step_regression.sh,                  [chmod +x test/backward_facing_step_regression.sh])
AC_CONFIG_FILES(test/locally_refine_regression.sh,                        [chmod +x test/locally_refine_regression.sh])
AC_CONFIG_FILES(test/redistribute_regression.sh,                          [chmod +x test/redistribute_regression.sh])
AC_CONFIG_FILES(test/input_files/backward_facing_step.in)
AC_CONFIG_FILES(test/input_files/locally_refine.in)
AC_CONFIG_FILES(test/input_files/redistribute.in)
AC_CONFIG_FILES(test/coupled_stokes_ns.sh,                                [chmod +x test/coupled_stokes_ns.sh])
AC_CONFIG_FILES(test/input_files/coupled_stokes_ns.in)
AC_CONFIG_FILES(test/cantera_chem_thermo_test.sh,	                  [chmod +x test/cantera_chem_thermo_test.sh])
AC_CONFIG_FILES(test/input_files/cantera_chem_thermo.in)
AC_CONFIG_FILES(test/reacting_low_mach_cantera_regression.sh,             [chmod +x test/reacting_low_mach_cantera_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_cantera_regression.in)
AC_CONFIG_FILES(test/input_files/air_2sp.xml)
AC_CONFIG_FILES(test/input_files/air_4sp.xml)
AC_CONFIG_FILES(test/cantera_mixture_unit.sh,                             [chmod +x test/cantera_mixture_unit.sh])
AC_CONFIG_FILES(test/cantera_transport_regression.sh,                     [chmod +x test/cantera_transport_regression.sh])
AC_CONFIG_FILES(test/input_files/cantera_transport.in)
AC_CONFIG_FILES(test/cantera_evaluator_regression.sh,                     [chmod +x test/cantera_evaluator_regression.sh])
AC_CONFIG_FILES(test/gas_recombination_catalytic_wall_unit_cantera.sh,    [chmod +x test/gas_recombination_catalytic_wall_unit_cantera.sh])
AC_CONFIG_FILES(test/gas_recombination_catalytic_wall_unit_antioch.sh,    [chmod +x test/gas_recombination_catalytic_wall_unit_antioch.sh])
AC_CONFIG_FILES(test/gas_solid_catalytic_wall_unit_cantera.sh,            [chmod +x test/gas_solid_catalytic_wall_unit_cantera.sh])
AC_CONFIG_FILES(test/gas_solid_catalytic_wall_unit_antioch.sh,            [chmod +x test/gas_solid_catalytic_wall_unit_antioch.sh])
AC_CONFIG_FILES(test/antioch_mixture_unit.sh,                             [chmod +x test/antioch_mixture_unit.sh])
AC_CONFIG_FILES(test/antioch_kinetics_regression.sh,                      [chmod +x test/antioch_kinetics_regression.sh])
AC_CONFIG_FILES(test/antioch_evaluator_regression.sh,                     [chmod +x test/antioch_evaluator_regression.sh])
AC_CONFIG_FILES(test/antioch_wilke_evaluator_regression.sh,               [chmod +x test/antioch_wilke_evaluator_regression.sh])
AC_CONFIG_FILES(test/input_files/antioch.in)
AC_CONFIG_FILES(test/input_files/gas_surface.in)


AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_constant_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_constant_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_constant_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_cea_constant_regression.sh, [chmod +x test/reacting_low_mach_antioch_cea_constant_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_cea_constant_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_constant_prandtl_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_constant_prandtl_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_constant_prandtl_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_cea_constant_prandtl_regression.sh, [chmod +x test/reacting_low_mach_antioch_cea_constant_prandtl_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_cea_constant_prandtl_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_catalytic_wall_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_catalytic_wall_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_catalytic_wall_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_arrhenius_catalytic_wall_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_arrhenius_catalytic_wall_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_arrhenius_catalytic_wall_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_power_catalytic_wall_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_power_catalytic_wall_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_power_catalytic_wall_regression.in)

AC_CONFIG_FILES(test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.sh, [chmod +x test/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.sh])
AC_CONFIG_FILES(test/input_files/reacting_low_mach_antioch_statmech_blottner_eucken_lewis_constant_gassolid_catalytic_wall_regression.in)

dnl-----------------------------------------------
dnl Generate run scripts for examples
dnl-----------------------------------------------
AC_CONFIG_FILES(examples/cavity_benchmark/run.sh,	     [chmod +x examples/cavity_benchmark/run.sh])
AC_CONFIG_FILES(examples/backward_facing_step/run.sh,	     [chmod +x examples/backward_facing_step/run.sh])
AC_CONFIG_FILES(examples/lid_driven_cavity/run.sh,	     [chmod +x examples/lid_driven_cavity/run.sh])
AC_CONFIG_FILES(examples/velocity_penalty/run.sh,	     [chmod +x examples/velocity_penalty/run.sh])
AC_CONFIG_FILES(examples/convection_cell/run.sh,	     [chmod +x examples/convection_cell/run.sh])
AC_CONFIG_FILES(examples/mass_injection/run.sh,	             [chmod +x examples/mass_injection/run.sh])
AC_CONFIG_FILES(examples/rayleigh_taylor/run.sh,	     [chmod +x examples/rayleigh_taylor/run.sh])
AC_CONFIG_FILES(examples/coupled_stokes_navierstokes/run.sh, [chmod +x examples/coupled_stokes_navierstokes/run.sh])
AC_CONFIG_FILES(examples/thermal_vortex/run.sh,              [chmod +x examples/thermal_vortex/run.sh])

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/grins/grins_version.h)
            
dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile doxygen/Makefile
		 doxygen/grins.dox examples/Makefile])

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl Hackish fix for Ubuntu/Debian libtool
perl -pi -e 's/link_all_deplibs=no/link_all_deplibs=yes/' libtool

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
