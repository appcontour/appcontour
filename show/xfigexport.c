#include <stdio.h>
#include <assert.h>
#include "showcontour.h"
#include "xfigexport.h"

#define XFIGMAXX 11000
#define XFIGMAXY  8000
#define XOFF 200
#define YOFF 200

//  #define conv(a) (int)((a->x - minx)*zoom) + XOFF,(int)((a->y - miny)*zoom) + YOFF
#define conv(a) (int)((a->x - minx)*zoom) + XOFF,(int)((maxy - a->y)*zoom) + YOFF
#define convx(ax) ((int)((ax - minx)*zoom) + XOFF)
#define convy(ay) ((int)((ay - miny)*zoom) + YOFF)
#define BUFLEN 200

void
xfig_export0 (struct polyline *contour, char *problem, struct grflags *grflags)
{
  static int xfig_count = 0;
  char defaultproblem[] = "contour";
  char buffer[BUFLEN];
  FILE *file;

  if (problem == 0)
  {
    problem = defaultproblem;
  }
  if (xfig_count > 0)
  {
    snprintf (buffer, BUFLEN, "%s_%05d.fig", problem, xfig_count);
  } else {
    snprintf (buffer, BUFLEN, "%s.fig", problem);
  }
  xfig_count++;
  file = fopen (buffer, "w");
  xfig_export (contour, file, grflags);
  fclose (file);
}

void
xfig_export (struct polyline *contour, FILE *file, struct grflags *grflags)
{
  struct line *line, *l, *last, *markline;
  struct vertex *v, *a, *b, *arrowa, *arrowb;
  struct rarc *arc;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;
  int count, markcount, w, st, cap_style, plotarrow;
  double arrowwidth, arrowheight;
  double stv, arrowthick;

  fprintf (file, "#FIG 3.2  produced by showcontour version %s\n", 
                 PACKAGE_VERSION);
  fprintf (file, "Landscape\nCenter\nInches\nLetter\n100.00\nSingle\n");
  fprintf (file, "-2\n1200 2\n");

  if (contour->isempty) {
    fprintf (file, "2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5\n");
    fprintf (file, "   150 150 8250 150 8250 8250 150 8250 150 150\n");
    fprintf (file, "2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5\n");
    fprintf (file, "   3450 3975 4725 3975 4725 4650 3450 4650 3450 3975\n");
    fprintf (file, "4 0 0 50 -1 0 24 0.0000 2 360 1110 3525 4425 Empty\\001\n");
  }
  //init_rarc (contour);
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
    assert (line->rarc);
    arc = line->rarc;
    if (arc->first != line && arc->loop != line) continue;
    /* let's walk along the arc */
    last = arc->last;
    if (last == 0)
    {
      last = arc->loop;
      assert (last);
      last = last->a->line[0];
    }
    w = (arc->d == 0) ? grflags->visiblewidth : grflags->invisiblewidth;
    stv = 0.0;
    st = arc->d;
    if (arc->d > 2) st = 2;
    if (st == 1) stv = grflags->dashlength;
    if (st == 2) stv = grflags->dotspacing;
    a = arrowa = line->a;
    arrowb = line->b;
    plotarrow = 1;
    if (arc->loop)
    {
      if ((a->type & V_TYPES) == V_CUSP) plotarrow = 0;
      fprintf (file, "2 3 %d %d 0 7 50 -1 -1 %5.3f 0 1 -1 0 0 %d\n", 
                    st, w, stv, arc->numsegments + 1);
    } else {
      plotarrow = 0;
      fprintf (file, "2 1 %d %d 0 7 50 -1 -1 %5.3f 0 1 -1 0 0 %d\n", 
                    st, w, stv, arc->numsegments + 1);
    }
    fprintf (file, "  %d %d\n", conv(a));
    l = markline = line;
    count = 0;
    markcount = arc->numsegments/3 + 1;
    while (1)
    {
//printf ("   another segment\n");
      b = l->b;
      if ((b->type & V_TYPES) == V_CUSP) plotarrow = 0;
      fprintf (file, "  %d %d\n", conv(b));
      count++;
      if (count == markcount) markline = l;
      if (l == last) break;
      l = b->line[1];
    }
    assert (count == arc->numsegments);
    if (plotarrow)  /* draw an arrow (a segment with 2 points) */
    {
      fprintf (file, "2 1 %d %d 0 7 49 -1 -1 %5.3f 0 1 -1 1 0 2\n", 
                    st, w, stv);
      arrowthick = w;
      arrowheight = 600.0;
      arrowwidth = 450.0;
      fprintf (file, "    0 0 %5.2f %8.2f %8.2f\n",
               arrowthick, arrowwidth, arrowheight);
      fprintf (file, "  %d %d\n", conv(arrowa));
      fprintf (file, "  %d %d\n", conv(arrowb));
    }
    /* print depth value */
    if (arc->d > st)
    {
      fprintf (file, "4 0 0 50 -1 0 24 0.0000 %d 135 450 ", 
        grflags->xfigspecial?2:0);
      fprintf (file, "%d %d ", conv(markline->a));
      if (grflags->xfigspecial) fprintf (file, "$");
      fprintf (file, "%d", arc->d);
      if (grflags->xfigspecial) fprintf (file, "$");
      fprintf (file, "\\001\n");
      //printf ("depth = %d\n", arc->d);
    }
  }
  for (v = contour->vertex; v; v = v->next)
  {
    switch (v->type & V_TYPES)
    {
      case V_CROSS:
      case V_CUSP:
        st = 0;
        w = grflags->pointsize;
        stv = 0.0;
        cap_style = 1;  /* round points */
        //cap_style = 2;  /* square points */
        fprintf (file, "2 1 %d %d %d 0 30 -1 -1 %5.3f 0 %d -1 0 0 1\n", 
                    st, w-4, XFIG_WHITE, stv, cap_style);
        fprintf (file, "  %d %d\n", conv(v));
        fprintf (file, "2 1 %d %d %d 0 31 -1 -1 %5.3f 0 %d -1 0 0 1\n", 
                    st, w, XFIG_BLACK, stv, cap_style);
        fprintf (file, "  %d %d\n", conv(v));
        break;
    }
  }
  return;
}
