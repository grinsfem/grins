# SYNOPSIS
#
#   Test for LIBMESH_DIM value
#
#   AX_LIBMESH_DIM()
#
# DESCRIPTION
#
#   Provides AC_SUBST(LIBMESH_DIM) depending on the value from the libMesh
#   installation.
#
# COPYLEFT
#
#   Copyright (c) 2014 Paul T. Bauman <pbauman@buffalo.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_LIBMESH_DIM],
[

   AC_REQUIRE([AX_PATH_LIBMESH_NEW])

   ac_libmeshdim_save_CPPFLAGS="$CPPFLAGS"
   CPPFLAGS="${LIBMESH_CPPFLAGS} ${CPPFLAGS}"

   AC_MSG_CHECKING(for libMesh dimension)

   AC_LANG_PUSH([C++])

   # First test if LIBMESH_DIM == 1
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
               @%:@include "libmesh/libmesh_config.h"
               ]], [[
               #if LIBMESH_DIM == 1
               /* Huzzah */
               #else
               #  error libmesh_dim != 1
               #endif
               ]])],[
                     AC_MSG_RESULT(1)
                     dim_one_succeeded=yes
                    ])

   
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
               @%:@include "libmesh/libmesh_config.h"
               ]], [[
               #if LIBMESH_DIM == 2
               /* Huzzah */
               #else
               #  error libmesh_dim != 2
               #endif
               ]])],[
                     AC_MSG_RESULT(2)
                     dim_two_succeeded=yes
                    ])

   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
               @%:@include "libmesh/libmesh_config.h"
               ]], [[
               #if LIBMESH_DIM == 3
               /* Huzzah */
               #else
               #  error libmesh_dim != 3
               #endif
               ]])],[
                     AC_MSG_RESULT(3)
                     dim_three_succeeded=yes
                    ])

   AC_LANG_POP([C++])

   CPPFLAGS="$ac_libmeshdim_save_CPPFLAGS"

   if test "x${dim_one_succeeded}" = "xyes"; then
      AC_SUBST(LIBMESH_DIM,1)
   fi

   if test "x${dim_two_succeeded}" = "xyes"; then
      AC_SUBST(LIBMESH_DIM,2)
   fi

   if test "x${dim_three_succeeded}" = "xyes"; then
      AC_SUBST(LIBMESH_DIM,3)
   fi

])
