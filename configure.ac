#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(appcontour, 0.7.1, paolini@dmf.unicatt.it)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([contour.c])
AC_CONFIG_HEADER([config.h])

AC_ARG_WITH([glut],
  AC_HELP_STRING([--with-glut], [glut engine for showcontour]),
  glut="yes"; nograph="no",
  glut="no")

AC_ARG_WITH([none],
  AC_HELP_STRING([--with-none], [no graphic engine for showcontour]),
  nograph="yes"; glut="no",
  nograph="no")

if (test "$glut" != "yes" -a "$nograph" != "yes"); then
  AC_MSG_ERROR([You must select an interface for showcontour to build.
     Possible options are:
       --with-glut    GLUT engine
       --with-none    no graphic engine]);
fi

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

if (test "$glut" = "yes"); then
  dnl check for GLUT headers.
  AC_CHECK_HEADERS(GL/glut.h, ,
    [AC_MSG_ERROR([No glut.h found, use --with-none.])])
  dnl check for GLUT.
  AC_CHECK_LIB(glut, glutMainLoop, ,
    [AC_MSG_ERROR([No glut library found, use --with-none.])])
fi

AM_CONDITIONAL(ENABLE_GLUT, test "$glut" = "yes")

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile show/Makefile])
AC_OUTPUT
