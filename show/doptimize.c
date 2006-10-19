#include <stdlib.h>
#include <stdio.h>
#include "showcontour.h"
#include "doptimize.h"

void
doptimize (struct polyline *contour)
{
  struct line *line;
  int goon = 1;

  while (goon)
  {
    goon = 0;
    for (line = contour->line; line; line = line->next)
    {
      goon += check_u_turn (contour, line);
    }
  }
}

int
check_u_turn (struct polyline *contour, struct line *line)
{
  struct line *l2, *l3, *l4;
  struct vertex *a, *b, *c, *d, *e;
  int linex, liney, l2x, l2y, l3x, l3y, l4x, l4y;
  double vec1, vec2, vec3;

  a = line->a;
  if (a->type == V_CROSS) return (0);
  b = line->b;
  if (b->type == V_CROSS) return (0);
  linex = b->x - a->x;
  liney = b->y - a->y;

  l2 = b->line[1];
  c = l2->b;
  if (c->type == V_CROSS) return (0);
  l2x = c->x - b->x;
  l2y = c->y - b->y;

  vec1 = - linex*l2y + liney*l2x;
  if (abs (vec1) != 1.0) return (0);
  l3 = c->line[1];
  d = l3->b;
  if (d->type == V_CROSS) return (0);
  l3x = d->x - c->x;
  l3y = d->y - c->y;

  vec2 = - l2x*l3y + l2y*l3x;

  if (vec1 != vec2) return (0);
  l4 = d->line[1];
  e = l4->b;
  l4x = e->x - d->x;
  l4y = e->y - d->y;

  vec3 = - l3x*l4y + l3y*l4x;
  if (vec3 == vec2) return (0);

  //printf ("found u turn...\n");
  /* modify contour, remove b, c; l2, l3 */
  line->b = d;
  d->line[0] = line;
  removenode (contour, b);
  removenode (contour, c);
  removeline (contour, l2);
  removeline (contour, l3);

  return (1);
}
