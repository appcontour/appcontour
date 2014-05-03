#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([appcontour],[2.2.1],[paolini@dmf.unicatt.it])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([contour.c])
AC_CONFIG_HEADER([config.h])

#glut="yes"
#gtk="yes"

AC_ARG_WITH([glut],
  AS_HELP_STRING([--without-glut],[no glut engine for showcontour]),
  [glut=="no"],
  [with_glut=check])

AC_ARG_WITH([gtk],
  AS_HELP_STRING([--without-gtk],[no gtk engine for showcontour]),
  [gtk="no"],
  [with_gtk=check])

#AC_ARG_WITH([none],
#  AS_HELP_STRING([--with-none],[no graphic engine for showcontour]),
#  nograph="yes"; gtk="no"; glut="no",
#  nograph="no")

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

GTKLIBS=
AS_IF([test "x$with_gtk" != xno],
      [AC_CHECK_LIB([gtk], [main],
              [AC_PATH_PROG(GTK_CONFIG, gtk-config, no, ${PATH})
               GTKLIBS="`${GTK_CONFIG} --libs`"
               GTK_CFLAGS="`${GTK_CONFIG} --cflags`"
               AC_DEFINE([HAVE_GTK], [1],
                         [Define if you have gtk])
               gtk="yes"
              ],
              [if test "x$with_gtk" != xcheck; then
                 AC_MSG_FAILURE(
                   [--with-gtk was given, but test for gtk failed])
               else
                 AC_MSG_WARN([gtk libraries not found, gtk graphic disabled])
               fi
              ], "")])

GLUTLIBS=
AS_IF([test "x$with_glut" != xno],
      [AC_CHECK_HEADERS(GL/glut.h,
          [AC_CHECK_LIB([glut], [main],
              [GLUTLIBS="-lglut -lGL"
               GLUT_CFLAGS="-I/usr/include/GL"
               AC_DEFINE([HAVE_GLUT], [1],[Define if you have glut])
               glut="yes"
              ],
              [if test "x$with_glut" != xcheck; then
                 AC_MSG_FAILURE(
                   [--with-glut was given, but test for glut failed])
               else
                 AC_MSG_WARN([glut libraries not found, glut graphic disabled])
               fi
              ], "")
          ],
          [if test "x$with_glut" != xcheck; then
                 AC_MSG_FAILURE(
                   [--with-glut was given, but test for glut failed])
           else
                 AC_MSG_WARN([glut libraries not found, glut graphic disabled])
           fi
          ])
       ])  

GRAPHIC_LIB="${GTKLIBS} ${GLUTLIBS}"
if (test "$GRAPHIC_LIB" = " "); then

  GRAPHIC_LIB="-lm"
fi

AC_SUBST([GRAPHIC_LIB])
AC_SUBST([CFLAGS],["${GTK_CFLAGS} ${GLUT_CFLAGS}"]) 
#AC_SUBST(GRAPHIC_LIB)
AM_CONDITIONAL(ENABLE_GLUT, test "$glut" = "yes")
AM_CONDITIONAL(ENABLE_GTK, test "$gtk" = "yes")

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile show/Makefile utils/Makefile examples/Makefile])
AC_OUTPUT
