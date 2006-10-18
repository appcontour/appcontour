bin_PROGRAMS = showcontour

if ENABLE_GLUT
GRAPHSOURCE = glutcontour.c
GRAPHLIB = -lglut
else
GRAPHSOURCE = nullcontour.c
GRAPHLIB = -lm
endif

showcontour_SOURCES = showcontour.c morseevent.c morseevent.h energy.c showcontour.h energy.h $(GRAPHSOURCE)
showcontour_LDADD = ../parser.o $(GRAPHLIB)

EXTRA_DIST = glutcontour.c nullcontour.c

AM_CFLAGS = -Wall
