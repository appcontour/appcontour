#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "showcontour.h"
#include "doptimize.h"

//#define DOPT_PLATEAU 1
//#define DOPT_PLATEAU_BACK 1  /* don't define both: an infinite loop results */
//#define DOPT_CHECK_CROSS_SLIDE  /* this is contained in lower_plateau */
#define DOPT_UTURN
#define DOPT_FLIP_LEFT
#define DOPT_CHECK_CROSS_TURN
#define DOPT_CHECK_CROSS_STAIRS
#define DOPT_LOWER_PLATEAU
#define DOPT_UNWIND

/* macro definitions */

#define TESTALIGNED(n1,n2,n3) \
  if (n3->x - n2->x != n2->x - n1->x) return (0);  \
  if (n3->y - n2->y != n2->y - n1->y) return (0);

#define SWAPCOORDS(n1,n2) \
  temp = n1->x;   \
  n1->x = n2->x;  \
  n2->x = temp;   \
  temp = n1->y;   \
  n1->y = n2->y;  \
  n2->y = temp;
#define RMN(node) removenode (contour, node)
#define RML(line) removeline (contour, line)
/*
 * the following works *only* if we are sure that
 * there is at least another line with the same earc
 * as ld, which is never the case since 'mergearcsinfo'
 * is not done yet.  This call should be replaced
 * by the following one.
 */
#define MVARC_DANGER(ld,lh) \
  ld->earc->refcount--;\
  assert (ld->earc->refcount > 0);\
  ld->earc = lh->earc;\
  ld->earc->refcount++;
/*
 * the third argument is a line in the same
 * arc as ld, into which inherit earc info
 * in case the refcount reaches 0
 */
#define MVARC(ld,lh,ea) \
  if (ld->earc->refcount <= 1) { \
    assert (ld->earc != ea); \
    ld->earc = mergearcinfo (ea, ld->earc); \
  } \
  ld->earc->refcount--;\
  assert (ld->earc->refcount > 0);\
  ld->earc = lh->earc;\
  ld->earc->refcount++;
#define CHAINFC2(bw,cross,l1)\
  if (bw) \
  {       \
    l1->b = cross;    \
  } else {            \
    l1->a = cross;    \
  }
#define CHAINFC3(bw,cross,l1,n1)\
  if (bw) \
  {       \
    l1->b = cross;    \
    l1->a = n1;       \
    n1->line[1] = l1; \
  } else {            \
    l1->a = cross;    \
    l1->b = n1;       \
    n1->line[0] = l1; \
  }
#define CHAINFC4(bw,cross,l1,n1,l2)\
  if (bw) \
  {       \
    l1->b = cross;    \
    l1->a = n1;       \
    n1->line[1] = l1; \
    n1->line[0] = l2; \
    l2->b = n1;       \
  } else {            \
    l1->a = cross;    \
    l1->b = n1;       \
    n1->line[0] = l1; \
    n1->line[1] = l2; \
    l2->a = n1;       \
  }
#define CHAINFC6(bw,cross,l1,n1,l2,n2,l3)\
  if (bw) \
  {       \
    l1->b = cross;    \
    l1->a = n1;       \
    n1->line[1] = l1; \
    n1->line[0] = l2; \
    l2->b = n1;       \
    l2->a = n2;       \
    n2->line[1] = l2; \
    n2->line[0] = l3; \
    l3->b = n2;       \
  } else {            \
    l1->a = cross;    \
    l1->b = n1;       \
    n1->line[0] = l1; \
    n1->line[1] = l2; \
    l2->a = n1;       \
    l2->b = n2;       \
    n2->line[0] = l2; \
    n2->line[1] = l3; \
    l3->a = n2;       \
  }

void
doptimize (struct polyline *contour)
{
  struct line *line;
  int goon = 1;
  int count = 0;
  //int goonl;

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
      //goonl = check_cross_turn (contour, line);
      //if (goonl < 0) return;
      //goon += goonl;
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
      //goonl = check_lower_plateau (contour, line);
      //if (goonl < 0) return;
      //goon += goonl;
      goon += check_lower_plateau (contour, line);
    }
#endif
#ifdef DOPT_UNWIND
    for (line = contour->line; line; line = line->next)
    {
      //goonl = check_unwind (contour, line);
      //if (goonl < 0) return;
      //goon += goonl;
      goon += check_unwind (contour, line);
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
  if (a->type & V_CROSS) return (0);
  b = line->b;
  if (b->type & V_CROSS) return (0);
  linex = b->x - a->x;
  liney = b->y - a->y;

  l2 = b->line[1];
  c = l2->b;
  if (c->type & V_CROSS) return (0);
  l2x = c->x - b->x;
  l2y = c->y - b->y;

  vec1 = - linex*l2y + liney*l2x;
  if (abs (vec1) != 1.0) return (0);
  l3 = c->line[1];
  d = l3->b;
  if (d->type & V_CROSS) return (0);
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
  RMN (b);
  RMN (c);
  RML (l2);
  RML (l3);

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
  if (a->type & V_CROSS) return (0);
  if (b->type & V_CROSS) return (0);
  lbc = b->line[1];
  c = lbc->b;
  if (c->type & V_CROSS) return (0);
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
  if (a->type & V_CROSS) return (0);
  b = back ? line->a : line->b;
  if (b->type & V_CROSS) return (0);
  linex = b->x - a->x;
  liney = b->y - a->y;

  l2 = b->line[1 - back];
  c = back ? l2->a : l2->b;
  if (c->type & V_CROSS) return (0);
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
    if (q->type & V_CROSS) break;
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
  struct vertex *a, *b1, *b2, *b3, *b4, *b5, *c1, *c2, *c3, *d;
  //struct vertex *e;
  struct line *lab1, *lac1, *lb1b2, *lb2b3, *lb3b4, *lb4b5;
  struct line *lc1c2, *lc2c3, *lad, *lae;
  int i1, i2, i3, i4;
  int lac1x, lac1y, lc1c2x, lc1c2y, vec1;
  int backward1 = 0;
  int backward2 = 0;
  //int backward3 = 0;
  int backward4 = 0;
  double temp;

  lae = line;
  a = lae->a;
  if ((a->type & V_CROSS) == 0)
  {
    backward4 = 1;
    a = lae->b;
  }
  if ((a->type & V_CROSS) == 0) return (0);
  //e = backward4 ? lae->a : lae->b;
  for (i4 = 0; i4 < 4; i4++)
  {
    if (a->line[i4] == lae) break;
  }
  assert (i4 < 4);
  i2 = 3 - i4;
  lac1 = a->line[i2];
  c1 = lac1->b;
  if (c1 == a) {backward2 = 1; c1 = lac1->a;}
  if (c1->type & V_CROSS) return (0);

  lc1c2 = c1->line[1-backward2];

  c2 = backward2 ? lc1c2->a : lc1c2->b;
  if (c2->type & V_CROSS) return (0);
  lc2c3 = c2->line[1-backward2];
  c3 = backward2 ? lc2c3->a : lc2c3->b;
  TESTALIGNED (c1, c2, c3);

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
  if (b1->type & V_CROSS) return (0);

  lb1b2 = b1->line[1-backward1];
  b2 = backward1 ? lb1b2->a : lb1b2->b;
  if (b2->type & V_CROSS) return (0);

  lad = a->line[i3];
  d = lad->b;
  if (d == a)
  {
    //backward3 = 1;
    d = lad->a;
  }
  if (d->type & V_CROSS) return (0);

  TESTALIGNED (b2, c1, c2);

  lb2b3 = b2->line[1-backward1];
  b3 = backward1 ? lb2b3->a : lb2b3->b;
  if (b3->type & V_CROSS) return (0);
  TESTALIGNED (b1, b2, b3);
  lb3b4 = b3->line[1-backward1];
  b4 = backward1 ? lb3b4->a : lb3b4->b;
  if (b4->type & V_CROSS) return (0);
  TESTALIGNED (a, c1, b4);
  lb4b5 = b4->line[1-backward1];
  b5 = backward1 ? lb4b5->a : lb4b5->b;
  if (b5->type & V_CROSS) return (0);
  TESTALIGNED (b3, b4, b5);

  /* all conditions are satisfied! */

  //printf ("found a corner turn (%lf,%lf)\n", a->x, a->y);

  SWAPCOORDS (a, c2);
  a->line[i4] = lac1;
  a->line[i1] = lb4b5;
  //lc2c3 = c2->line[1-backward2];
  a->line[i2] = lc2c3;
  // a->line[i3] does not change
  CHAINFC2 (backward1, a, lb4b5);
  CHAINFC2 (backward2, a, lc2c3);
  CHAINFC6 (backward4, a, lac1, c1, lc1c2, c2, lae);
  MVARC_DANGER (lc1c2, lae);
  MVARC_DANGER (lac1, lae);
  RMN (b1);
  RMN (b2);
  RMN (b3);
  RMN (b4);
  RML (lab1);
  RML (lb1b2);
  RML (lb2b3);
  RML (lb3b4);

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
  if ((a->type & V_CROSS) == 0)
  {
    backward1 = 1;
    a = line->b;
  }
  if ((a->type & V_CROSS) == 0) return (0);
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
  if (b1->type & V_CROSS) return (0);
  if (c1->type & V_CROSS) return (0);
  if (d1->type & V_CROSS) return (0);

  lc1c2 = c1->line[1-backward2];
  lb1b2 = b1->line[1-backward3];
  ld1d2 = d1->line[1-backward4];

  c2 = backward2 ? lc1c2->a : lc1c2->b;
  b2 = backward3 ? lb1b2->a : lb1b2->b;
  d2 = backward4 ? ld1d2->a : ld1d2->b;

  if (c2->x - c1->x != b1->x - a->x) return (0);
  if (c2->y - c1->y != b1->y - a->y) return (0);
  TESTALIGNED (a, b1, b2);
  if (d2->x - d1->x != b1->x - a->x) return (0);
  if (d2->y - d1->y != b1->y - a->y) return (0);

  //printf ("found a cross slide (%lf,%lf)\n", a->x, a->y);

  SWAPCOORDS (a, b1);
  a->line[i1] = lab1;
  a->line[i2] = lc1c2;
  a->line[i3] = lb1b2;
  a->line[i4] = ld1d2;

  MVARC_DANGER (lab1, line);

  CHAINFC4 (backward1, a, lab1, b1, line);
  CHAINFC2 (backward2, a, lc1c2);
  CHAINFC2 (backward3, a, lb1b2);
  CHAINFC2 (backward4, a, ld1d2);

  RMN (c1);
  RMN (d1);
  RML (lac1);
  RML (lad1);

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
  if ((a->type & V_CROSS) == 0)
  {
    backward1 = 1;
    a = line->b;
  }
  if ((a->type & V_CROSS) == 0) return (0);
  for (i1 = 0; i1 < 4; i1++)
  {
    if (a->line[i1] == line) break;
  }
  i3 = 3 - i1;
  lab1 = a->line[i3];
  b1 = lab1->b;
  if (b1 == a) {backward3 = 1; b1 = lab1->a;}
  if (b1->type & V_CROSS) return (0);
  lb1b2 = b1->line[1-backward3];
  b2 = backward3 ? lb1b2->a : lb1b2->b;
  if (b2->type & V_CROSS) return (0);

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
  if (c1->type & V_CROSS) return (0);
  if (d1->type & V_CROSS) return (0);

  lc1c2 = c1->line[1-backward2];
  ld1d2 = d1->line[1-backward4];

  c2 = backward2 ? lc1c2->a : lc1c2->b;
  d2 = backward4 ? ld1d2->a : ld1d2->b;
  if (c2->type & V_CROSS) return (0);
  if (d2->type & V_CROSS) return (0);
  TESTALIGNED (b2, b1, c2);
  TESTALIGNED (a, d1, d2);

  lc2c3 = c2->line[1-backward2];
  lb2b3 = b2->line[1-backward3];
  ld2d3 = d2->line[1-backward4];

  c3 = backward2 ? lc2c3->a : lc2c3->b;
  d3 = backward4 ? ld2d3->a : ld2d3->b;
  if (c3->type & V_CROSS) return (0);
  if (d3->type & V_CROSS) return (0);
  TESTALIGNED (c1, c2, c3);
  TESTALIGNED (b1, b2, d3);

  lc3c4 = c3->line[1-backward2];
  c4 = backward2 ? lc3c4->a : lc3c4->b;
  if (c4->type & V_CROSS) return (0);
  TESTALIGNED (a, b1, c4);

  //printf ("stairs found (%lf,%lf)\n", a->x, a->y);

  /* changing coordinates */

  SWAPCOORDS (a, b2);

  c3->x = b1->x;
  c3->y = b1->y;

  b1->x = d1->x;
  b1->y = d1->y;

  /* changing topological links */

  a->line[i1] = lab1;
  assert (a->line[i2] == lac1);
  a->line[i3] = lb2b3;
  a->line[i4] = ld2d3;

  CHAINFC6 (backward1, a, lab1, b1, lb1b2, b2, line);
  CHAINFC3 (backward2, a, lac1, c3)
  CHAINFC2 (backward3, a, lb2b3);
  CHAINFC2 (backward4, a, ld2d3);

  /* maintain arc references and refcount */
  MVARC_DANGER (lab1, line);
  MVARC_DANGER (lb1b2, line);

  /* rimozioni */
  RMN (c1);
  RMN (c2);
  RMN (d1);
  RMN (d2);
  RML (lc1c2);
  RML (lc2c3);
  RML (lad1);
  RML (ld1d2);

  return (1);
}

int
check_lower_plateau (struct polyline *contour, struct line *line)
{
  struct vertex *a, *b, *aback, *bprev, *bt, *btt;
  struct line *l, *lback, *lend, *l1, *l2, *l3;
  int abx, aby, px, py, qx, qy, bprevx, bprevy;
  int i, i1, i2, i3, i4, vec, vecg, count;
  int backward2;
  //int crossingfound = 0;

  a = line->a;
  if (a->type & V_CROSS) return (0);
  b = line->b;
  //if (b->type & V_CROSS) return (0);

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
    if (b->type & V_CROSS)
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
      if (bt->type & V_CROSS) return (0);
      if (backward2)
      {
        btt = bt->line[0]->a;
      } else {
        btt = bt->line[1]->b;
      }
      TESTALIGNED (b, bt, btt);
      //crossingfound = 1;
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
      if (vecg == vec && (bprev->type & V_CROSS)) return (0);
      if (vecg == vec)
      {
        if (b->type & V_CROSS) return (0);
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
    if (b->type & V_CROSS)
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
      assert ((bt->type & V_CROSS) == 0);
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
      //MVARC_DANGER (l2, l1);
      MVARC (l2, l1, l3->earc);
    } else {
      l = b->line[1];
    }
    b = l->b;
  }
  /* change topological info at start */
  line->a = lback->a;
  line->a->line[1] = line;
  RMN (a);
  RML (lback);
  /* change topological info at end */
  lend = b->line[1];
  l->b = lend->b;
  l->b->line[0] = l;
  RMN (b);
  RML (lend);

  //if (crossingfound)
  //{
    //printf ("crossing... (%lf,%lf)\n", line->b->x, line->b->y);
    //return (-1);
  //}
  return (1);
}

int
check_unwind (struct polyline *contour, struct line *lb)
{
  struct vertex *a, *b, *c, *d, *e, *f, *g, *h;
  struct line *la, *lc, *ld, *le, *lf, *lg, *lh, *li;
  int backward1 = 0;
  int backward2 = 0;
  int backward3 = 0;
  int backward4 = 0;
  int i1, i2, i3, i4;
  int dex, dey, gdx, gdy, vec;
  double temp;

  e = lb->a;
  d = lb->b;
  if ((e->type & V_CROSS) == 0)
  {
    backward1 = 1;
    e = lb->b;
    d = lb->a;
  }
  if ((e->type & V_CROSS) == 0) return (0);
  if (d->type & V_CROSS) return (0);
  la = d->line[1-backward1];
  g = backward1 ? la->a : la->b;
  if (g->type & V_CROSS) return (0);

  dex = e->x - d->x;
  dey = e->y - d->y;
  gdx = d->x - g->x;
  gdy = d->y - g->y;

  vec = dex*gdy - dey*gdx;
  if (vec == 0) return (0);

  for (i1 = 0; i1 < 4; i1++)
  {
    if (e->line[i1] == lb) break;
  }
  i3 = 3 - i1;
  i2 = i1i2[i1];
  if (vec < 0) i2 = i1i2[i3];
  i4 = 3 - i2;

  lc = e->line[i2];
  lg = e->line[i3];
  lh = e->line[i4];

  a = lc->b;
  if (a == e) {backward2 = 1; a = lc->a;}
  if (a->type & V_CROSS) return (0);
  ld = a->line[1-backward2];
  b = backward2 ? ld->a : ld->b;
  if (b->type & V_CROSS) return (0);
  TESTALIGNED (g, a, b);
  le = b->line[1-backward2];
  c = backward2 ? le->a : le->b;
  if (c->type & V_CROSS) return (0);
  TESTALIGNED (a, b, c);
  lf = c->line[1-backward2];
  h = backward2 ? lf->a : lf->b;
  if (h->type & V_CROSS) return (0);

  f = lg->b;
  if (f == e) {backward3 = 1; f = lg->a;}
  if (f->type & V_CROSS) return (0);
  li = f->line[1-backward3];

  TESTALIGNED (e, f, h);

  if (lh->b == e) backward4 = 1;

  //printf ("found an unwind at (%lf,%lf)\n", e->x, e->y);

  /* changing coordinates */

  SWAPCOORDS (e, f);

  /* changing topological links */

  assert (e->line[i1] == lb);
  e->line[i2] = lf;
  e->line[i3] = li;
  e->line[i4] = ld;

  CHAINFC6 (backward1, e, lb, b, lc, a, la);
  CHAINFC2 (backward2, e, lf);
  CHAINFC2 (backward3, e, li);
  CHAINFC4 (backward4, e, ld, f, lh);

  /* maintain arc references and refcount */

  MVARC (lc, la, ld->earc);
  MVARC (ld, lh, le->earc);

  /* removals */

  RMN (c);
  RMN (d);
  RML (le);
  RML (lg);

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
