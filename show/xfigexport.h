#include "grcommon.h"

#define XFIG_BLACK 0
#define XFIG_BLUE 1
#define XFIG_GREEN 2
#define XFIG_CYAN 3
#define XFIG_RED 4
#define XFIG_MAGENTA 5
#define XFIG_YELLOW 6
#define XFIG_WHITE 7

void xfig_export0 (struct polyline *contour, char *filen, struct grflags *f);
void xfig_export (struct polyline *contour, FILE *file, struct grflags *f);
