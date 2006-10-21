#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "showcontour.h"
#include "doptimize.h"

#define DOPT_UTURN 1
//#define DOPT_PLATEAU 1
//#define DOPT_PLATEAU_BACK 1  /* don't define both: an infinite loop results */
#define DOPT_FLIP_LEFT 1
#define DOPT_CHECK_CROSS_TURN 1
#define DOPT_CHECK_CROSS_SLIDE 1

void
doptimize (struct polyline *contour)
{
  struct line *line;
  int goon = 1;
  int count = 0;

  while (goon && count++ < 10000)
  {
    goon = 0;
#ifdef DOPT_UTURN
    for (line = contour->line; line; line = line->next)
    {
      goon += check_u_turn (contour, line);
    }
#endif
#ifdef DOPT_FLIP_LEFT
    for (line = contour->line; line; line = line->next)
    {
      goon += check_flip_left (contour, line);
    }
#endif
#ifdef DOPT_PLATEAU
    for (line = contour->line; line; line = line->next)
    {
      goon += check_plateau (contour, line, 0);
    }
#endif
#ifdef DOPT_PLATEAU_BACK
    for (line = contour->line; line; line = line->next)
    {
      goon += check_plateau (contour, line, 1);
    }
#endif
#ifdef DOPT_CHECK_CROSS_TURN
    for (line = contour->line; line; line = line->next)
    {
      goon += check_cross_turn (contour, line);
    }
#endif
#ifdef DOPT_CHECK_CROSS_SLIDE
    for (line = contour->line; line; line = line->next)
    {
      goon += check_cross_slide (contour, line);
      //if (check_cross_slide (contour, line)) return;
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
check_flip_left (struct polyline *contour, struct line *lab)
{
  struct vertex *a, *b, *c;
  struct line *lbc;
  int ax, ay, abx, aby, bcx, bcy, vec;

  a = lab->a;
  b = lab->b;
  if (a->type == V_CROSS) return (0);
  if (b->type == V_CROSS) return (0);
  lbc = b->line[1];
  c = lbc->b;
  if (c->type == V_CROSS) return (0);
  ax = a->x;
  ay = a->y;
  abx = b->x - ax;
  aby = b->y - ay;
  bcx = c->x - b->x;
  bcy = c->y - b->y;
  vec = aby*bcx - abx*bcy;
  if (vec != 1) return (0);

  if (site_occupied (contour, ax + bcx, ay + bcy)) return (0);

  /* can move node */
  b->x = ax + bcx;
  b->y = ay + bcy;
  return (1);
}

int
check_plateau (struct polyline *contour, struct line *line, int back)
{
  struct vertex *a, *b, *c, *p, *q;
  struct line *l2, *ll;
  int linex, liney, l2x, l2y, llx, lly, i;
  int vec1, vecgen, dimplateau, dx, dy, px, py;

  a = back ? line->b : line->a;
  if (a->type == V_CROSS) return (0);
  b = back ? line->a : line->b;
  if (b->type == V_CROSS) return (0);
  linex = b->x - a->x;
  liney = b->y - a->y;

  l2 = b->line[1 - back];
  c = back ? l2->a : l2->b;
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
    ll = p->line[1 - back];
    q = back ? ll->a : ll->b;
    llx = q->x - p->x;
    lly = q->y - p->y;
    vecgen = - l2x*lly + l2y*llx;
    p = q;
    if (q->type == V_CROSS) break;
    if (vecgen == vec1 || vecgen == - vec1) break;
  }
  dimplateau = i;
  if (dimplateau <= 1) return (0);

  dx = l2x - linex;
  dy = l2y - liney;

  p = b;
  i = 0;
  while (1)
  {
    ll = p->line[1 - back];
    q = back ? ll->a : ll->b;
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
    //printf ("trovato plateau: %d:%d (%d,%d) (%lf, %lf)\n", dimplateau, i, dx, dy, b->x, b->y);
    p = q;
  }
  assert (0);
  return (0);
}

/*
 * make a cross turn around a corner
 */

static int i1i2[] = {1,3,0,2};

int
check_cross_turn (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b1, *b2, *b3, *b4, *b5, *c1, *c2, *d, *e;
  struct line *lab1, *lac1, *lb1b2, *lb2b3, *lb3b4, *lb4b5;
  struct line *lc1c2, *lc2c3, *lad, *lae;
  int i1, i2, i3, i4;
  int lac1x, lac1y, lc1c2x, lc1c2y, vec1;
  int backward1 = 0;
  int backward2 = 0;
  int backward3 = 0;
  int backward4 = 0;
  double temp;

  lae = line;
  a = lae->a;
  if (a->type != V_CROSS)
  {
    backward4 = 1;
    a = lae->b;
  }
  if (a->type != V_CROSS) return (0);
  e = backward4 ? lae->a : lae->b;
  for (i4 = 0; i4 < 4; i4++)
  {
    if (a->line[i4] == lae) break;
  }
  assert (i4 < 4);
  i2 = 3 - i4;
  lac1 = a->line[i2];
  c1 = lac1->b;
  if (c1 == a) {backward2 = 1; c1 = lac1->a;}
  if (c1->type == V_CROSS) return (0);

  lc1c2 = c1->line[1-backward2];

  c2 = backward2 ? lc1c2->a : lc1c2->b;
  if (c2->type == V_CROSS) return (0);

  lac1x = c1->x - a->x;
  lac1y = c1->y - a->y;
  lc1c2x = c2->x - c1->x;
  lc1c2y = c2->y - c1->y;
  vec1 = lac1x*lc1c2y - lac1y*lc1c2x;
  if (vec1 == 0) return (0);

  i1 = (vec1 == 1) ? i1i2[i4] : i1i2[i2];
  i3 = 3 - i1;

  lab1 = a->line[i1];
  b1 = lab1->b;
  if (b1 == a) {backward1 = 1; b1 = lab1->a;}
  if (b1->type == V_CROSS) return (0);

  lb1b2 = b1->line[1-backward1];
  b2 = backward1 ? lb1b2->a : lb1b2->b;
  if (b2->type == V_CROSS) return (0);

  lad = a->line[i3];
  d = lad->b;
  if (d == a) {backward3 = 1; d = lad->a;}
  if (d->type == V_CROSS) return (0);

  if (c2->x - c1->x != c1->x - b2->x) return (0);
//printf ("quo c1 (%lf,%lf), c2 (%lf,%lf), b2 (%lf,%lf)\n", c1->x, c1->y, c2->x, c2->y, b2->x, b2->y);
  if (c2->y - c1->y != c1->y - b2->y) return (0);

  lb2b3 = b2->line[1-backward1];
  b3 = backward1 ? lb2b3->a : lb2b3->b;
  if (b3->type == V_CROSS) return (0);
  if (b3->x - b2->x != b2->x - b1->x) return (0);
  if (b3->y - b2->y != b2->y - b1->y) return (0);
  lb3b4 = b3->line[1-backward1];
  b4 = backward1 ? lb3b4->a : lb3b4->b;
  if (b4->type == V_CROSS) return (0);
  if (b4->x - c1->x != c1->x - a->x) return (0);
  if (b4->y - c1->y != c1->y - a->y) return (0);
  lb4b5 = b4->line[1-backward1];
  b5 = backward1 ? lb4b5->a : lb4b5->b;
  if (b5->type == V_CROSS) return (0);
  if (b5->x - b4->x != b4->x - b3->x) return (0);
  if (b5->y - b4->y != b4->y - b3->y) return (0);

  /* all conditions are satisfied! */

  //printf ("found a corner turn (%lf,%lf)\n", a->x, a->y);

  temp = a->x;
  a->x = c2->x;
  c2->x = temp;
  temp = a->y;
  a->y = c2->y;
  c2->y = temp;
  a->line[i4] = lac1;
  a->line[i1] = lb4b5;
  lc2c3 = c2->line[1-backward2];
  a->line[i2] = lc2c3;
  // a->line[i3] does not change
  if (backward2) lc2c3->b = a; else lc2c3->a = a;
  if (backward1) lb4b5->b = a; else lb4b5->a = a;
  if (backward4)
  {
    lac1->b = a;
    lac1->a = c1;
    lc1c2->b = c1;
    lc1c2->a = c2;
    lae->b = c2;
    c1->line[1] = lac1;
    c1->line[0] = lc1c2;
    c2->line[1] = lc1c2;
    c2->line[0] = lae;
  } else {
    lac1->a = a;
    lac1->b = c1;
    lc1c2->a = c1;
    lc1c2->b = c2;
    lae->a = c2;
    c1->line[0] = lac1;
    c1->line[1] = lc1c2;
    c2->line[0] = lc1c2;
    c2->line[1] = lae;
  }
  lc1c2->arc->refcount--;
  lac1->arc->refcount--;
  assert (lc1c2->arc->refcount > 0);
  assert (lac1->arc->refcount > 0);
  lc1c2->arc = lae->arc;
  lac1->arc = lae->arc;
  lae->arc->refcount += 2;

  removenode (contour, b1);
  removenode (contour, b2);
  removenode (contour, b3);
  removenode (contour, b4);
  removeline (contour, lab1);
  removeline (contour, lb1b2);
  removeline (contour, lb2b3);
  removeline (contour, lb3b4);

  return (1);
}

int
check_cross_slide (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b1, *b2, *c1, *c2, *d1, *d2;
  struct line *lac1, *lab1, *lad1, *lc1c2, *lb1b2, *ld1d2;
  int i1, i2, i3, i4;
  int backward1 = 0;
  int backward2 = 0;
  int backward3 = 0;
  int backward4 = 0;
  double temp;

  a = line->a;
  if (a->type != V_CROSS)
  {
    backward1 = 1;
    a = line->b;
  }
  if (a->type != V_CROSS) return (0);
  for (i1 = 0; i1 < 4; i1++)
  {
    if (a->line[i1] == line) break;
  }
  i3 = 3 - i1;
  i2 = i1i2[i1];
  i4 = 3 - i2;

  lac1 = a->line[i2];
  lab1 = a->line[i3];
  lad1 = a->line[i4];
  c1 = lac1->b;
  if (c1 == a) {backward2 = 1; c1 = lac1->a;}
  b1 = lab1->b;
  if (b1 == a) {backward3 = 1; b1 = lab1->a;}
  d1 = lad1->b;
  if (d1 == a) {backward4 = 1; d1 = lad1->a;}
  if (b1->type == V_CROSS) return (0);
  if (c1->type == V_CROSS) return (0);
  if (d1->type == V_CROSS) return (0);

  lc1c2 = c1->line[1-backward2];
  lb1b2 = b1->line[1-backward3];
  ld1d2 = d1->line[1-backward4];

  c2 = backward2 ? lc1c2->a : lc1c2->b;
  b2 = backward3 ? lb1b2->a : lb1b2->b;
  d2 = backward4 ? ld1d2->a : ld1d2->b;

  if (c2->x - c1->x != b1->x - a->x) return (0);
  if (c2->y - c1->y != b1->y - a->y) return (0);
  if (b2->x - b1->x != b1->x - a->x) return (0);
  if (b2->y - b1->y != b1->y - a->y) return (0);
  if (d2->x - d1->x != b1->x - a->x) return (0);
  if (d2->y - d1->y != b1->y - a->y) return (0);

  //printf ("found a cross slide (%lf,%lf)\n", a->x, a->y);

  temp = a->x;
  a->x = b1->x;
  b1->x = temp;
  temp = a->y;
  a->y = b1->y;
  b1->y = temp;
  a->line[i1] = lab1;
  a->line[i2] = lc1c2;
  a->line[i3] = lb1b2;
  a->line[i4] = ld1d2;
  lab1->arc->refcount--;
  assert (lab1->arc->refcount > 0);
  lab1->arc = line->arc;
  line->arc->refcount++;
  if (backward2) lc1c2->b = a;
    else lc1c2->a = a;
  if (backward4) ld1d2->b = a;
    else ld1d2->a = a;
  if (backward3) lb1b2->b = a;
    else lb1b2->a = a;
  if (backward1)
  {
    lab1->b = a;
    lab1->a = b1;
    b1->line[0] = line;
    b1->line[1] = lab1;
    line->b = b1;
  } else {
    lab1->a = a;
    lab1->b = b1;
    b1->line[0] = lab1;
    b1->line[1] = line;
    line->a = b1;
  }

  removenode (contour, c1);
  removenode (contour, d1);
  removeline (contour, lac1);
  removeline (contour, lad1);
  return (1);
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
