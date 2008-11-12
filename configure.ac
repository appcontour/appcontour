#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(appcontour, 1.0.0, paolini@dmf.unicatt.it)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([contour.c])
AC_CONFIG_HEADER([config.h])

#glut="yes"
gtk="yes"

AC_ARG_WITH([glut],
  AC_HELP_STRING([--without-glut], [no glut engine for showcontour]),
  glut="no",
  glut="yes")

AC_ARG_WITH([gtk],
  AC_HELP_STRING([--without-gtk], [no gtk engine for showcontour]),
  gtk="no",
  gtk="yes")

#AC_ARG_WITH([none],
#  AC_HELP_STRING([--with-none], [no graphic engine for showcontour]),
#  nograph="yes"; gtk="no"; glut="no",
#  nograph="no")

if (test "$glut" != "yes" -a "$nograph" != "yes" -a "$gtk" != "yes"); then
  AC_MSG_ERROR([You must select an interface for showcontour to build.
     Possible options are:
       --with-gtk     GTK engine
       --with-glut    GLUT engine
       --with-none    no graphic engine]);
fi

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

GRAPHIC_LIB=""
if (test "$gtk" = "yes"); then
  dnl check for GTK gtk-config.
  AC_PATH_PROG(GTK_CONFIG, gtk-config, no, ${PATH})
  if (test "${GTK_CONFIG}" = "no"); then
     AC_MSG_ERROR([No gtk library found, use --with-none or --with-glut.])
  fi
  CFLAGS="${CFLAGS} `${GTK_CONFIG} --cflags`"
  GRAPHIC_LIB="${GRAPHIC_LIB} `${GTK_CONFIG} --libs`"
  AC_DEFINE(HAVE_GTK, 1, [You have a usable gtk library.])
fi

if (test "$glut" = "yes"); then
  dnl check for GLUT headers.
  AC_CHECK_HEADERS(GL/glut.h, ,
    [AC_MSG_ERROR([No glut.h found, use --with-none or --with-gtk.])])
  dnl check for GLUT.
  AC_CHECK_LIB(glut, glutMainLoop, [GRAPHIC_LIB="${GRAPHIC_LIB} -lglut"],
    [AC_MSG_ERROR([No glut library found, use --with-none or --with-gtk.])])
  AC_DEFINE(HAVE_GLUT, 1, [You have a usable glut library.])
fi

if (test "$GRAPHIC_LIB" = ""); then
  GRAPHIC_LIB="-lm"
fi
AC_SUBST(GRAPHIC_LIB)
AM_CONDITIONAL(ENABLE_GLUT, test "$glut" = "yes")
AM_CONDITIONAL(ENABLE_GTK, test "$gtk" = "yes")

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile show/Makefile utils/Makefile])
AC_OUTPUT
