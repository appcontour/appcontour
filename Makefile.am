SUBDIRS = show utils examples

bin_PROGRAMS = contour
contour_SOURCES = contour.c readsketch.c readmorse.c rules.c sketch.c parser.c appcontour.c printmorse.c knot2morse.c mendes.c invrules.c fundamental.c alexander.c giovecanonify.c contour.h parser.h mendes.h fundamental.h alexander.h giovecanonify.h
EXTRA_DIST = tests/test.* BUGS utils/*.sh TODO
AM_CFLAGS = -Wall -DEXAMPLES_DIR=\"$(pkgdatadir)/examples\"
AM_LDFLAGS = -ftrapv
# man1_MANS = contour.1
dist_man_MANS = contour.1

TESTS = test
