#include <stdio.h>
#include <assert.h>
#include "showcontour.h"
#include "grcommon.h"
#include "xfigexport.h"

#define XFIGMAXX 11000
#define XFIGMAXY  8000
#define XOFF 200
#define YOFF 200

#define conv(a) (int)((a->x - minx)*zoom) + XOFF,(int)((a->y - miny)*zoom) + YOFF

void
xfig_export (struct polyline *contour, FILE *file, struct grflags *grflags)
{
  struct line *line, *l, *last, *markline;
  struct vertex *v, *a, *b;
  struct rarc *arc;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;
  int count, markcount, w;

  fprintf (file, "#FIG 3.2\n");
  fprintf (file, "Landscape\nCenter\nInches\nLetter\n100.00\nSingle\n");
  fprintf (file, "-2\n1200 2\n");

  init_rarc (contour);
  maxx = -1000.0;
  maxy = -1000.0;
  minx =  1000.0;
  miny =  1000.0;
  for (v = contour->vertex; v; v = v->next)
  {
    if (maxx < v->x) maxx = v->x;
    if (maxy < v->y) maxy = v->y;
    if (minx > v->x) minx = v->x;
    if (miny > v->y) miny = v->y;
  }
  xmed = (maxx + minx)/2.0;
  ymed = (maxy + miny)/2.0;
  zoomx = XFIGMAXX/(maxx - minx);
  zoomy = XFIGMAXY/(maxy - miny);
  zoom = zoomx;
  if (zoom > zoomy) zoom = zoomy;

  for (line = contour->line; line; line = line->next)
  {
    arc = line->rarc;
    if (arc->first != line && arc->loop != line) continue;
    /* let's walk along the arc */
//printf ("walking (%d)...\n", arc->numsegments);
    last = arc->last;
    if (last == 0)
    {
//printf ("isloop\n");
      last = arc->loop;
      assert (last);
      last = last->a->line[0];
    }
//printf ("last defined\n");
    w = (arc->d == 0) ? 2 : 1;
    fprintf (file, "2 1 0 %d 0 7 50 -1 -1 0.000 0 0 -1 0 0 %d\n", 
                    w, arc->numsegments + 1);
    a = line->a;
    fprintf (file, "  %d %d\n", conv(a));
    l = line;
    count = 0;
    markcount = count/3;
    while (1)
    {
//printf ("   another segment\n");
      b = l->b;
      fprintf (file, "  %d %d\n", conv(b));
      count++;
      if (count == markcount) markline = l;
      if (l == last) break;
      l = b->line[1];
    }
    assert (count == arc->numsegments);
    /* print depth value */
    if (arc->d > 0) printf ("depth = %d\n", arc->d);
  }
  return;
}
