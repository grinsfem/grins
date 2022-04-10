dnl
dnl AM_PATH_CPPUNIT([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
AC_DEFUN([AM_PATH_CPPUNIT],
[
  dnl Check the user's environment for CPPUNIT_DIR and if it's set, use it,
  dnl otherwise set some defaults that maybe will work and then later check
  dnl if the user passed some options to configure
  AS_IF([test "x$CPPUNIT_DIR" = "x"],
        [dnl Defaults that might work if cppunit headers are in /usr/include
         dnl and libraries are in /usr/lib, i.e. standard installation locations.
         CPPUNIT_CPPFLAGS=
         CPPUNIT_LIBS=-lcppunit],
        [CPPUNIT_CPPFLAGS="-I$CPPUNIT_DIR/include"
         CPPUNIT_LIBS="-L$CPPUNIT_DIR/lib -lcppunit"])

  dnl User can override with explicit values for:
  dnl --with-cppunit-prefix
  dnl or
  dnl --with-cppunit-include
  dnl --with-cppunit-lib
  AC_ARG_WITH([cppunit-prefix],
              [AS_HELP_STRING([--with-cppunit-prefix=PATH],
                              [Specify a path for cppunit installation])],
              [CPPUNIT_PREFIX=$withval
               CPPUNIT_CPPFLAGS="-I$CPPUNIT_PREFIX/include"
               CPPUNIT_LIBS="-L$CPPUNIT_PREFIX/lib -lcppunit"])

  AC_ARG_WITH(cppunit-include,
              AS_HELP_STRING([--with-cppunit-include=PATH],
                             [Specify a path for cppunit header files]),
              CPPUNIT_CPPFLAGS="-I$withval")

  AC_ARG_WITH(cppunit-lib,
              AS_HELP_STRING([--with-cppunit-lib=PATH],
                             [Specify a path for cppunit libs]),
              CPPUNIT_LIBS="-L$withval -lcppunit")


  AC_MSG_CHECKING(whether we can build a trivial CppUnit program)
  AC_LANG_PUSH([C++])
  saveCPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$saveCPPFLAGS $CPPUNIT_CPPFLAGS"
  saveLIBS="$LIBS"
  LIBS="$CPPUNIT_LIBS $saveLIBS"

  AC_LINK_IFELSE([AC_LANG_SOURCE([[
  @%:@include <cppunit/ui/text/TestRunner.h>
  int main(int argc, char **argv)
  {
    CppUnit::TextUi::TestRunner runner;

    if (runner.run())
      return 0;

    return 1;
  }
  ]])],[
    AC_MSG_RESULT(yes)
    ifelse([$1], , :, [$1])
    AC_SUBST(HAVE_CPPUNIT,[1])
    AC_DEFINE([HAVE_CPPUNIT], [1], [Enable CPPUnit Tests])
  ],[
    AC_MSG_RESULT(no)
    CPPUNIT_CPPFLAGS=
    CPPUNIT_LIBS=
    ifelse([$2], , :, [$2])
  ])

  AC_LANG_POP
  LIBS="$saveLIBS"
  CPPFLAGS="$saveCPPFLAGS"

  AC_SUBST(CPPUNIT_CPPFLAGS)
  AC_SUBST(CPPUNIT_LIBS)
])
