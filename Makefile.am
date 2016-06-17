SUBDIRS = show utils examples data

bin_PROGRAMS = contour
contour_SOURCES = contour.c readsketch.c readmorse.c readdtcode.c rules.c sketch.c parser.c appcontour.c printmorse.c knot2morse.c mendes.c invrules.c fundamental.c fox.c alexander.c giovecanonify.c laurent.c contour.h parser.h readdtcode.h mendes.h fundamental.h fox.h alexander.h giovecanonify.h laurent.h rolfsen_to_dt.h
EXTRA_DIST = test tests/test.* BUGS utils/*.sh TODO
AM_CFLAGS = -Wall -DEXAMPLES_DIR=\"$(pkgdatadir)/examples\" -DPKGDATA_DIR=\"$(pkgdatadir)\"
AM_LDFLAGS = -ftrapv
# man1_MANS = contour.1
dist_man_MANS = contour.1

TESTS = test
