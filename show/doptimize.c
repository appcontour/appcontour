#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "showcontour.h"
#include "doptimize.h"

#define DOPT_UTURN
//#define DOPT_PLATEAU 1
//#define DOPT_PLATEAU_BACK 1  /* don't define both: an infinite loop results */
#define DOPT_FLIP_LEFT
#define DOPT_CHECK_CROSS_TURN
//#define DOPT_CHECK_CROSS_SLIDE  /* this is contained in lower_plateau */
#define DOPT_CHECK_CROSS_STAIRS
#define DOPT_LOWER_PLATEAU

void
doptimize (struct polyline *contour)
{
  struct line *line;
  int goon = 1;
  int count = 0;
  int goonl;

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
    }
#endif
#ifdef DOPT_CHECK_CROSS_STAIRS
    for (line = contour->line; line; line = line->next)
    {
      goon += check_cross_stairs (contour, line);
    }
#endif
#ifdef DOPT_LOWER_PLATEAU
    for (line = contour->line; line; line = line->next)
    {
      goonl = check_lower_plateau (contour, line);
      if (goonl < 0) return;
      goon += goonl;
      //goon += check_lower_plateau (contour, line);
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

int
check_cross_stairs (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b1, *b2, *c1, *c2, *c3, *c4, *d1, *d2, *d3;
  struct line *lab1, *lac1, *lad1, *lb1b2, *lc1c2, *ld1d2;
  struct line *lb2b3, *lc2c3, *ld2d3, *lc3c4;
  int i1, i2, i3, i4, temp;
  int backward1 = 0;
  int backward2 = 0;
  int backward3 = 0;
  int backward4 = 0;
  int lab1x, lab1y, lb1b2x, lb1b2y, vec;

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
  lab1 = a->line[i3];
  b1 = lab1->b;
  if (b1 == a) {backward3 = 1; b1 = lab1->a;}
  if (b1->type == V_CROSS) return (0);
  lb1b2 = b1->line[1-backward3];
  b2 = backward3 ? lb1b2->a : lb1b2->b;
  if (b2->type == V_CROSS) return (0);

  lab1x = b1->x - a->x;
  lab1y = b1->y - a->y;
  lb1b2x = b2->x - b1->x;
  lb1b2y = b2->y - b1->y;
  vec = lab1x*lb1b2y - lab1y*lb1b2x; 
  /* vec = 1 for descending stairs (ab1=(1,0), b1b2=(0,1)) */

  if (vec == 0) return (0);
  assert (abs(vec) == 1);
  if (vec > 0) i2 = i1i2[i1]; else i2 = i1i2[i3];

  i4 = 3 - i2;

  lac1 = a->line[i2];
  lad1 = a->line[i4];
  c1 = lac1->b;
  if (c1 == a) {backward2 = 1; c1 = lac1->a;}
  d1 = lad1->b;
  if (d1 == a) {backward4 = 1; d1 = lad1->a;}
  if (c1->type == V_CROSS) return (0);
  if (d1->type == V_CROSS) return (0);

  lc1c2 = c1->line[1-backward2];
  ld1d2 = d1->line[1-backward4];

  c2 = backward2 ? lc1c2->a : lc1c2->b;
  d2 = backward4 ? ld1d2->a : ld1d2->b;
  if (c2->type == V_CROSS) return (0);
  if (d2->type == V_CROSS) return (0);
  if (c2->x - b1->x != b1->x - b2->x) return (0);
  if (c2->y - b1->y != b1->y - b2->y) return (0);
  if (d2->x - d1->x != d1->x - a->x) return (0);
  if (d2->y - d1->y != d1->y - a->y) return (0);

  lc2c3 = c2->line[1-backward2];
  lb2b3 = b2->line[1-backward3];
  ld2d3 = d2->line[1-backward4];

  c3 = backward2 ? lc2c3->a : lc2c3->b;
  d3 = backward4 ? ld2d3->a : ld2d3->b;
  if (c3->type == V_CROSS) return (0);
  if (d3->type == V_CROSS) return (0);
  if (c3->x - c2->x != c2->x - c1->x) return (0);
  if (c3->y - c2->y != c2->y - c1->y) return (0);
  if (d3->x - b2->x != b2->x - b1->x) return (0);
  if (d3->y - b2->y != b2->y - b1->y) return (0);

  lc3c4 = c3->line[1-backward2];
  c4 = backward2 ? lc3c4->a : lc3c4->b;
  if (c4->type == V_CROSS) return (0);
  if (c4->x - b1->x != b1->x - a->x) return (0);
  if (c4->y - b1->y != b1->y - a->y) return (0);

  //printf ("stairs found (%lf,%lf)\n", a->x, a->y);

  /* changing coordinates */

  temp = a->x;
  a->x = b2->x;
  b2->x = temp;
  temp = a->y;
  a->y = b2->y;
  b2->y = temp;

  c3->x = b1->x;
  c3->y = b1->y;

  b1->x = d1->x;
  b1->y = d1->y;

  /* changing topological links */

  a->line[i1] = lab1;
  assert (a->line[i2] == lac1);
  a->line[i3] = lb2b3;
  a->line[i4] = ld2d3;

  if (backward2)
  {
    c3->line[1] = lac1;
    lac1->a  = c3;
    assert (lac1->b == a);
  } else {
    c3->line[0] = lac1;
    lac1->b  = c3;
    assert (lac1->a == a);
  }

  if (backward3) lb2b3->b = a; else lb2b3->a = a;
  if (backward4) ld2d3->b = a; else ld2d3->a = a;

  if (backward1)
  {
    lab1->b = a;
    lab1->a = b1;
    lb1b2->b = b1;
    lb1b2->a = b2;
    line->b = b2;
    b1->line[1] = lab1;
    b1->line[0] = lb1b2;
    b2->line[1] = lb1b2;
    b2->line[0] = line;
  } else {
    lab1->a = a;
    lab1->b = b1;
    lb1b2->a = b1;
    lb1b2->b = b2;
    line->a = b2;
    b1->line[0] = lab1;
    b1->line[1] = lb1b2;
    b2->line[0] = lb1b2;
    b2->line[1] = line;
  }

  /* maintain arc references and refcount */
  lab1->arc->refcount--;
  lb1b2->arc->refcount--;
  assert (lab1->arc->refcount > 0);
  assert (lb1b2->arc->refcount > 0);
  lab1->arc = line->arc;
  lb1b2->arc = line->arc;
  line->arc->refcount += 2;

  /* rimozioni */
  removenode (contour, c1);
  removenode (contour, c2);
  removenode (contour, d1);
  removenode (contour, d2);
  removeline (contour, lc1c2);
  removeline (contour, lc2c3);
  removeline (contour, lad1);
  removeline (contour, ld1d2);

  return (1);
}

int
check_lower_plateau (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b, *aback, *bprev, *bt, *btt;
  struct line *l, *lback, *lend, *l1, *l2, *l3;
  int abx, aby, px, py, qx, qy, bprevx, bprevy;
  int i, i1, i2, i3, i4, vec, vecg, count;
  int backward2, crossingfound = 0;

  a = line->a;
  if (a->type == V_CROSS) return (0);
  b = line->b;
  //if (b->type == V_CROSS) return (0);

  abx = b->x - a->x;
  aby = b->y - a->y;

  lback = a->line[0];
  aback = lback->a;

  px = a->x - aback->x;
  py = a->y - aback->y;

  vec = abx*py - aby*px;
  if (vec == 0) return (0);
  count = 0;

  l = line;
  while (1)
  {
    count++;
    if (b->type == V_CROSS)
    {
      for (i1 = 0; i1 < 4; i1++)
      {
        if (b->line[i1] == l) break;
      }
      assert (i1 < 4);
      i3 = 3 - i1;
      l = b->line[i3];
      assert (l->a == b);
      i2 = i1i2[i1];
      if (vec < 0) i2 = i1i2[i3];
      bt = b->line[i2]->b;
      backward2 = 0;
      if (bt == b) {backward2 = 1; bt = b->line[i2]->a;}
      if (bt->type == V_CROSS) return (0);
      if (backward2)
      {
        btt = bt->line[0]->a;
      } else {
        btt = bt->line[1]->b;
      }
      if (btt->x - bt->x != bt->x - b->x) return (0);
      if (btt->y - bt->y != bt->y - b->y) return (0);
      crossingfound = 1;
      assert (l->b != b);
      b = l->b;
    } else {
      l = b->line[1];
      bprev = b;
      b = l->b;

      qx = b->x - bprev->x;
      qy = b->y - bprev->y;
      vecg = qx*aby - qy*abx;
      if (vecg == -vec) return (0);  /* not a plateau */
      if (vecg == vec && bprev->type == V_CROSS) return (0);
      if (vecg == vec)
      {
        if (b->type == V_CROSS) return (0);
        break;        /* this IS a plateau */
      }

      bprevx = bprev->x;
      bprevy = bprev->y;
      if (site_occupied (contour, bprevx - px, bprevy - py)) return (0);
                 /* there is an occupied cell */
    }
  }

  if (count < 2) return (0);

  //printf ("plateau found at (%lf,%lf) of size %d\n", a->x, a->y, count);

  /* change coordinates */
  b = line->b;
  l = line;
  for (i = 1; i < count; i++)
  {
    b->x -= px;
    b->y -= py;
    if (b->type == V_CROSS)
    {
      for (i1 = 0; i1 < 4; i1++)
      {
        if (b->line[i1] == l) break;
      }
      i3 = 3 - i1;
      l = b->line[i3];
      i2 = i1i2[i1];
      if (vec < 0) i2 = i1i2[i3];
      i4 = 3 - i2;
      l2 = b->line[i2];
      bt = l2->b;
      backward2 = 0;
      if (bt == b) {backward2 = 1; bt = b->line[i2]->a;}
      assert (bt->type != V_CROSS);
      bt->x += px;
      bt->y += py;
      l1 = b->line[i4];
      l3 = bt->line[1-backward2];
      if (backward2)
      {
        assert (l1->a == b);
        l1->a = bt;
        l2->b = bt;
        l2->a = b;
        l3->b = b;
      } else {
        assert (l1->b == b);
        l1->b = bt;
        l2->a = bt;
        l2->b = b;
        l3->a = b;
      }
      bt->line[backward2] = l1;
      bt->line[1-backward2] = l2;
      b->line[i2] = l3;
      b->line[i4] = l2;
      l2->arc->refcount--;
      l2->arc = l1->arc;
      l2->arc->refcount++;
    } else {
      l = b->line[1];
    }
    b = l->b;
  }
  /* change topological info at start */
  line->a = lback->a;
  line->a->line[1] = line;
  removenode (contour, a);
  removeline (contour, lback);
  /* change topological info at end */
  lend = b->line[1];
  l->b = lend->b;
  l->b->line[0] = l;
  removenode (contour, b);
  removeline (contour, lend);

  //if (crossingfound)
  //{
    //printf ("crossing... (%lf,%lf)\n", line->b->x, line->b->y);
    //return (-1);
  //}
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
