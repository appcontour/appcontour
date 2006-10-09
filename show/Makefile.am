bin_PROGRAMS = showcontour

if ENABLE_GLUT
showcontour_SOURCES = showcontour.c glutcontour.c
showcontour_LDADD = ../parser.o -lglut
else
showcontour_SOURCES = showcontour.c nullcontour.c
showcontour_LDADD = ../parser.o -lm
endif

AM_CFLAGS = -Wall
