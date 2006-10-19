#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "showcontour.h"
#include "doptimize.h"

#define DOPT_UTURN 1
#define DOPT_PLATEAU 1
#define DOPT_PLATEAU_BACK 1

void
doptimize (struct polyline *contour)
{
  struct line *line;
  int goon = 1;

  while (goon)
  {
    goon = 0;
#ifdef DOPT_UTURN
    for (line = contour->line; line; line = line->next)
    {
      goon += check_u_turn (contour, line);
    }
#endif
#ifdef DOPT_PLATEAU
    for (line = contour->line; line; line = line->next)
    {
      goon += check_plateau (contour, line);
    }
#endif
#ifdef DOPT_PLATEAU_BACK
    for (line = contour->line; line; line = line->next)
    {
      goon += check_plateau_back (contour, line);
    }
#endif
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

int
check_plateau (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b, *c, *p, *q;
  struct line *l2, *ll;
  int linex, liney, l2x, l2y, llx, lly, i;
  int vec1, vecgen, dimplateau, dx, dy, px, py;

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
  p = c;
  i = 0;
  while (1)
  {
    i++;
    ll = p->line[1];
    q = ll->b;
    llx = q->x - p->x;
    lly = q->y - p->y;
    vecgen = - l2x*lly + l2y*llx;
    p = q;
    if (vecgen == - vec1) return (0);
    if (q->type == V_CROSS) break;
    if (vecgen == vec1) break;
  }
  dimplateau = i;
  if (dimplateau <= 1) return (0);

  dx = l2x - linex;
  dy = l2y - liney;

  printf ("trovato plateau: %d (%lf, %lf)\n", dimplateau, b->x, b->y);

  p = b;
  i = 0;
  while (1)
  {
    ll = p->line[1];
    q = ll->b;
    px = p->x;
    py = p->y;
    if (site_occupied (contour, px + dx, py + dy))
    {
      /* site is not free, cannot move node p */
      return (i);
    }
    i++;
    if (i >= dimplateau) return (i - 1);
    /* can move node p */
    p->x += dx;
    p->y += dy;
    p = q;
  }
  assert (0);
  return (0);
}

int
check_plateau_back (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b, *c, *p, *q;
  struct line *l2, *ll;
  int linex, liney, l2x, l2y, llx, lly, i;
  int vec1, vecgen, dimplateau, dx, dy, px, py;

  a = line->b;
  if (a->type == V_CROSS) return (0);
  b = line->a;
  if (b->type == V_CROSS) return (0);
  linex = b->x - a->x;
  liney = b->y - a->y;

  l2 = b->line[0];
  c = l2->a;
  if (c->type == V_CROSS) return (0);
  l2x = c->x - b->x;
  l2y = c->y - b->y;
  vec1 = - linex*l2y + liney*l2x;
  if (abs (vec1) != 1.0) return (0);
  p = c;
  i = 0;
  while (1)
  {
    i++;
    ll = p->line[0];
    q = ll->a;
    llx = q->x - p->x;
    lly = q->y - p->y;
    vecgen = - l2x*lly + l2y*llx;
    p = q;
    if (vecgen == - vec1) return (0);
    if (q->type == V_CROSS) break;
    if (vecgen == vec1) break;
  }
  dimplateau = i;
  if (dimplateau <= 1) return (0);

  dx = l2x - linex;
  dy = l2y - liney;

  printf ("trovato backw plateau: %d (%lf, %lf)\n", dimplateau, b->x, b->y);

  p = b;
  i = 0;
  while (1)
  {
    ll = p->line[0];
    q = ll->a;
    px = p->x;
    py = p->y;
    if (site_occupied (contour, px + dx, py + dy))
    {
      /* site is not free, cannot move node p */
      return (i);
    }
    i++;
    if (i >= dimplateau) return (i - 1);
    /* can move node p */
    p->x += dx;
    p->y += dy;
    p = q;
  }
  assert (0);
  return (0);
}

struct vertex *
site_occupied (struct polyline *contour, int px, int py)
{
  struct vertex *v;

  for (v = contour->vertex; v; v = v->next)
  {
    if (px == (int) v->x && py == (int) v->y) return (v);
  }
  return (0);
}
