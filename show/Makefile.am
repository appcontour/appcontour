bin_PROGRAMS = showcontour

if ENABLE_GLUT
GLUTGRAPHSOURCE = glutcontour.c
else
GLUTGRAPHSOURCE = 
endif

if ENABLE_GTK
GTKGRAPHSOURCE = gtkcontour.c
else
GTKGRAPHSOURCE = 
endif

showcontour_SOURCES = showcontour.c grcommon.c grcommon.h doptimize.c doptimize.h morseevent.c morseevent.h energy.c xfigexport.c xfigexport.h showcontour.h energy.h $(GLUTGRAPHSOURCE) $(GTKGRAPHSOURCE)
showcontour_LDADD = ../parser.o @GRAPHIC_LIB@

AM_CFLAGS = -Wall
