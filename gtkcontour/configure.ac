#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(gtkmorse, 0.1.1, pasqui@dmf.unicatt.it)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([gtkcontour.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
#CONTOUR="no"
AC_CHECK_PROG(CONTOUR,contour,trovato,no,$PATH)
if (test "$CONTOUR" = "trovato"); then
  dnl ho trovato contour.
  AC_PATH_PROG(PATH_CONTOUR,contour,no)
  AC_DEFINE(HAVE_CONTOUR, 1,[Ho trovato contour])
  AC_DEFINE_UNQUOTED(PATH_CONTOUR,"${PATH_CONTOUR} -q iscontour",[path di contour])
fi

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

LIBS=""
  dnl check for GTK gtk-config.
  AC_PATH_PROG(GTK_CONFIG, pkg-config, no, ${PATH})
  if (test "${GTK_CONFIG}" = "no"); then
     AC_MSG_ERROR([No pkg-config program found.])
  fi

  pkg-config --exists gtk+-2.0 
  if (test $? != "0"); then
     AC_MSG_ERROR([No gtk library found, gtk-devel installed?.])
  fi
  CFLAGS="${CFLAGS} `${GTK_CONFIG} gtk+-2.0 --cflags`"
  LIBS="${LIBS} `${GTK_CONFIG} gtk+-2.0 --libs`"
  AC_DEFINE(HAVE_GTK, 1, [You have a usable gtk library.])
AC_SUBST(LIBS)

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
