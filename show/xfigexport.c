#include <stdio.h>
#include "showcontour.h"

void
xfig_export (struct polyline *contour, FILE *file)
{
  struct line *line;

  fprintf (file, "#FIG 3.2\n");
  fprintf (file, "Landscape\nCenter\nInches\nLetter\n100.00\nSingle\n");
  fprintf (file, "-2\n1200 2\n");

  for (line = contour->line; line; line = line->next)
  {
    fprintf (file, "2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n");
  }
  return;
}
