#include "grcommon.h"

#define XFIG_RED 4
#define XFIG_BLUE 7

void xfig_export0 (struct polyline *contour, char *filen, struct grflags *f);
void xfig_export (struct polyline *contour, FILE *file, struct grflags *f);
