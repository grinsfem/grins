# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2010-03-24
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo Revision id................... : $BUILD_VERSION
echo
echo Library Dependencies:
echo libMesh....................... : $LIBMESH_PREFIX
echo libMesh CXXFLAGS.............. : $LIBMESH_CXXFLAGS
echo libMesh INCLUDE............... : $LIBMESH_INCLUDE
echo libMesh LDFLAGS............... : $LIBMESH_LDFLAGS
echo libMesh LIBS.................. : $LIBMESH_LIBS
#echo QUESO......................... : $QUESO_PREFIX
#echo Trilinos...................... : $TRILINOS_PREFIX
#echo HDF5.......................... : $HDF5_PREFIX
#echo GSL........................... : $GSL_PREFIX
#echo GLPK.......................... : $GLPK_PREFIX
echo Boost......................... : $BOOST_ROOT

# masa optional check:
echo
echo Optional Features:
if test "$HAVE_ANTIOCH" = "0"; then
  echo '   'Link with Antioch............. : no
else
  echo '   'Link with Antioch............. : $ANTIOCH_PREFIX
fi
if test "$HAVE_CANTERA" = "0"; then
  echo '   'Link with Cantera............. : no
else
  echo '   'Link with Cantera............. : $CANTERA_PREFIX
fi
if test "$HAVE_GRVY" = "0"; then
  echo '   'Link with GRVY................ : no
else
  echo '   'Link with GRVY................ :$GRVY_PREFIX
fi
if test "x$USE_GRVY_TIMERS" = "x1"; then
  echo '   'Use GRVY timers............... : yes
else
  echo '   'Use GRVY timers............... : no
fi
if test "$HAVE_MASA" = "0"; then
  echo '   'Link with MASA................ : no
else
  echo '   'Link with MASA................ : $MASA_PREFIX
fi
if test "$HAVE_GCOV_TOOLS" = "0"; then
  echo '   'Enable gcov code coverage..... : no
else     
  echo '   'Enable gcov code coverage..... : yes
fi


echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
