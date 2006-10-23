bin_PROGRAMS = showcontour

if ENABLE_GLUT
GRAPHSOURCE = glutcontour.c
GRAPHLIB = -lglut
else
if ENABLE_GTK
GRAPHSOURCE = gtkcontour.c
else
GRAPHSOURCE = nullcontour.c
GRAPHLIB = -lm
endif
endif

showcontour_SOURCES = showcontour.c grcommon.c grcommon.h doptimize.c doptimize.h morseevent.c morseevent.h energy.c showcontour.h energy.h $(GRAPHSOURCE)
showcontour_LDADD = ../parser.o $(GRAPHLIB)

EXTRA_DIST = glutcontour.c nullcontour.c

AM_CFLAGS = -Wall
