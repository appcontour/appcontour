#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(appcontour, 0.7.1, paolini@dmf.unicatt.it)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([contour.c])
AC_CONFIG_HEADER([config.h])

AC_ARG_ENABLE(glut,
[  --enable-glut    Use glut graphics for show],
[case "${enableval}" in
  yes) glut=true ;;
  no)  glut=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-glut) ;;
esac],[glut=false])
AM_CONDITIONAL(GLUT, test x$glut = xtrue)

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile show/Makefile])
AC_OUTPUT
