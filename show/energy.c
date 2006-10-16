#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "showcontour.h"
#include "energy.h"

#define K1_COEFF 1.0     /* perimeter */
#define K2_COEFF 0.1     /* k^2 */
#define K3_COEFF 0.7     /* cross+cusp repulsion */
#define K4_COEFF 1.0     /* arc repulsion */
#define K5_COEFF 0.0     /* |k| */

#define ALLOW_REPULSION 1

#define NODE_SEP 0.4
#define ORTO_HEAVINESS 2.0
#define CUSP_PRECISION 0.8           /* 1 means perfect cusp (180 degrees) */
#define SMOOTH_K5 0.001

/*
 * a value too near 1 hase the effect of shooting out (sometimes) cusps towards
 * infinity, probably a side effect of the addition of nodes for elongating
 * segments.  The energy associated with cusps should be modified in order
 * to avoid such tendency
 */

#ifdef CHECK_GRADIENT
void check_gradient (struct polyline *contour);
#endif

//extern int test;
extern double time;
//extern double timerrep;
//extern double taurep;

void
energyinit (void)
{
  allowrepulsion = ALLOW_REPULSION;
  k1_coeff = K1_COEFF;
  k2_coeff = K2_COEFF;
  k3_coeff = K3_COEFF;
  k4_coeff = K4_COEFF;
  k5_coeff = K5_COEFF;
}

#ifdef CHECK_GRADIENT
void
check_gradient (struct polyline *contour)
{
  double baseenergy, normgradsq, energy, *gradx, *grady;
  double dx, saved, nder, der, normgrad, relerr;
  struct vertex *v;

  dx = 1e-7;
  compute_gradient (contour);  /* anticipate gradient computation */
  gradx = contour->gradx;
  grady = contour->grady;

  normgradsq = normsq (gradx, contour->numvertices) +
               normsq (grady, contour->numvertices);
  normgrad = sqrt (normgradsq);
  compute_repulsive_energy (contour);
  baseenergy = compute_energy (contour);
  for (v = contour->vertex; v; v = v->next)
  {
    saved = v->x;
    v->x += dx;
    compute_repulsive_energy (contour);
    energy = compute_energy (contour);
    v->x = saved;
    nder = (energy - baseenergy)/dx;
    der = gradx[v->tag];
    relerr = fabs ((der - nder)/normgrad);
    if (relerr > 0.1)
      printf ("der/nder = %lf\n", der/nder);
    saved = v->y;
    v->y += dx;
    compute_repulsive_energy (contour);
    energy = compute_energy (contour);
    v->y = saved;
    nder = (energy - baseenergy)/dx;
    der = grady[v->tag];
    relerr = fabs ((der - nder)/normgrad);
    if (relerr > 0.01)
      printf ("nodo %d, der/nder = %lf\n", v->tag, der/nder);
  }
}
#endif

//static double *rgradx = 0;
//static double *rgrady = 0;
//static double renergy = 0.0;
//static int rgraddim = 0;

/* qui ci sono i contributi delle tre energie */

double
compute_energy (struct polyline *contour)
{
  struct line *line;
  struct vertex *a, *p, *b, *v1, *v2;
  double dx, dy, dist;
  double alpha, la, lb;
  double energy1 = 0.0;
  double energy2 = 0.0;
  double energy3 = 0.0;
  double energy4 = 0.0;
  double energy5 = 0.0;
  int i, j;

  /* primo contributo, dovuto al perimetro */
 
  for (line = contour->line; line; line = line->next)
  {
    a = line->a;
    b = line->b;
    dx = a->x - b->x;
    dy = a->y - b->y;
    energy1 += sqrt (dx*dx + dy*dy);
  }

  energy1 *= k1_coeff;

  /* secondo contributo, dovuto alla k^2 */

  for (p = contour->vertex; p; p = p->next)
  {
    switch (p->type)
    {
      case V_REGULAR:
      case V_CUSP:
        a = p->line[0]->a;
        b = p->line[1]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        if (p->type == V_CUSP)
        {
          alpha -= CUSP_PRECISION*M_PI;
        }
        energy2 += 2*alpha*alpha/(la + lb);
        break;

      case V_CROSS:
        a = p->line[0]->a;
        if (a == p) a = p->line[0]->b;
        b = p->line[3]->a;
        if (b == p) b = p->line[3]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        energy2 += 2*alpha*alpha/(la + lb);

        a = p->line[1]->a;
        if (a == p) a = p->line[1]->b;
        b = p->line[2]->a;
        if (b == p) b = p->line[2]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        energy2 += 2*alpha*alpha/(la + lb);

        a = p->line[0]->a;
        if (a == p) a = p->line[0]->b;
        b = p->line[1]->a;
        if (b == p) b = p->line[1]->b;
        alpha = get_alpha (a, p, b, &la, &lb) - M_PI/2;
        energy2 += 2*ORTO_HEAVINESS*alpha*alpha/(la + lb);
        break;
    }
  }
  energy2 *= k2_coeff;

  /* terzo contributo: repulsione tra nodi/cuspidi */
  /*
   * la funzione energia utilizzata e'
   *  f(d) = (NODE_SEP - d)^2/d   se d < NODE_SEP
   */

  for (i = 0; i < contour->specnodesnum; i++)
  {
    v1 = contour->specnodes[i];
    for (j = i+1; j < contour->specnodesnum; j++)
    {
      v2 = contour->specnodes[j];
      dx = v1->x - v2->x;
      dy = v1->y - v2->y;
      dist = sqrt (dx*dx + dy*dy);
      if (dist < NODE_SEP)
      {
        energy3 += k3_coeff*repulsive_field (dist);
      }
    }
  }

  /* quarto contributo: repulsione tra tutti i nodi */

  if (allowrepulsion)
  {
    energy4 = k4_coeff*contour->renergy;
  }

  /* fifth contribution: integral of |k| */

  for (p = contour->vertex; p; p = p->next)
  {
    switch (p->type)
    {
      case V_REGULAR:
      case V_CUSP:
        a = p->line[0]->a;
        b = p->line[1]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        energy5 += k5_coeff*fabs(alpha);
        break;

      case V_CROSS:
        a = p->line[0]->a;
        if (a == p) a = p->line[0]->b;
        b = p->line[3]->a;
        if (b == p) b = p->line[3]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        energy5 += k5_coeff*fabs(alpha);

        a = p->line[1]->a;
        if (a == p) a = p->line[1]->b;
        b = p->line[2]->a;
        if (b == p) b = p->line[2]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        energy5 += k5_coeff*fabs(alpha);
        break;
    }
  }
 
  return (energy1 + energy2 + energy3 + energy4 + energy5);
}

/* here we compute the gradient of the energy */

void
compute_gradient (struct polyline *contour)
{
  struct line *line;
  struct vertex *a, *p, *b, *q, *v1, *v2;
  int i, j, tag1, tag2;
  double dx, dy, len, xel, yel, distsq, dist, force;
  double alpha, gcx, gcy, la, lb, sigma;
  double *gradx, *grady;
  double apx, apy, pqx, pqy, qbx, qby, vec1, vec2, pqsq;

  gradx = contour->gradx;
  grady = contour->grady;
  for (i = 0; i < contour->numvertices; i++)
    gradx[i] = grady[i] = 0.0;

  /* primo contributo, dovuto al perimetro */
 
  for (line = contour->line; line; line = line->next)
  {
    a = line->a;
    b = line->b;
    dx = a->x - b->x;
    dy = a->y - b->y;
    len = sqrt (dx*dx + dy*dy);
    xel = k1_coeff*dx/len;
    yel = k1_coeff*dy/len;
    gradx[a->tag] += xel;
    grady[a->tag] += yel;
    gradx[b->tag] -= xel;
    grady[b->tag] -= yel;
  }

  /* secondo contributo, dovuto alla k^2 */

  for (p = contour->vertex; p; p = p->next)
  {
    switch (p->type)
    {
      case V_REGULAR:
      case V_CUSP:
        a = p->line[0]->a;
        b = p->line[1]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        if (p->type == V_CUSP)
        {
          alpha -= CUSP_PRECISION*M_PI;
          //printf ("node with cusp...\n");
        }
        grad_curv (b, p, alpha, lb, la, &gcx, &gcy);
        gcx *= k2_coeff;
        gcy *= k2_coeff;
        gradx[b->tag] += gcx;
        grady[b->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;

        grad_curv (a, p, -alpha, la, lb, &gcx, &gcy);
        gcx *= k2_coeff;
        gcy *= k2_coeff;
        gradx[a->tag] += gcx;
        grady[a->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;
        break;

      case V_CROSS:
        a = p->line[0]->a;
        if (a == p) a = p->line[0]->b;
        b = p->line[3]->a;
        if (b == p) b = p->line[3]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        grad_curv (b, p, alpha, lb, la, &gcx, &gcy);
        gcx *= k2_coeff;
        gcy *= k2_coeff;
        gradx[b->tag] += gcx;
        grady[b->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;

        grad_curv (a, p, -alpha, la, lb, &gcx, &gcy);
        gcx *= k2_coeff;
        gcy *= k2_coeff;
        gradx[a->tag] += gcx;
        grady[a->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;

        a = p->line[1]->a;
        if (a == p) a = p->line[1]->b;
        b = p->line[2]->a;
        if (b == p) b = p->line[2]->b;
        alpha = get_alpha (a, p, b, &la, &lb);
        grad_curv (b, p, alpha, lb, la, &gcx, &gcy);
        gcx *= k2_coeff;
        gcy *= k2_coeff;
        gradx[b->tag] += gcx;
        grady[b->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;

        grad_curv (a, p, -alpha, la, lb, &gcx, &gcy);
        gcx *= k2_coeff;
        gcy *= k2_coeff;
        gradx[a->tag] += gcx;
        grady[a->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;

        a = p->line[0]->a;
        if (a == p) a = p->line[0]->b;
        b = p->line[1]->a;
        if (b == p) b = p->line[1]->b;
        alpha = get_alpha (a, p, b, &la, &lb) - M_PI/2;
        grad_curv (b, p, alpha, lb, la, &gcx, &gcy);
        gcx *= ORTO_HEAVINESS*k2_coeff;
        gcy *= ORTO_HEAVINESS*k2_coeff;
        gradx[b->tag] += gcx;
        grady[b->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;

        grad_curv (a, p, -alpha, la, lb, &gcx, &gcy);
        gcx *= ORTO_HEAVINESS*k2_coeff;
        gcy *= ORTO_HEAVINESS*k2_coeff;
        gradx[a->tag] += gcx;
        grady[a->tag] += gcy;
        gradx[p->tag] -= gcx;
        grady[p->tag] -= gcy;
        break;
    }
  }

  /* terzo contributo, distanza tra i nodi non regolari */

  for (i = 0; i < contour->specnodesnum; i++)
  {
    v1 = contour->specnodes[i];
    tag1 = v1->tag;
    for (j = i+1; j < contour->specnodesnum; j++)
    {
      v2 = contour->specnodes[j];
      tag2 = v2->tag;
      dx = v1->x - v2->x;
      dy = v1->y - v2->y;
      distsq = dx*dx + dy*dy;
      dist = sqrt (distsq);
      force = k3_coeff*repulsive_force (dist);
      gradx[tag1] += dx*force;
      gradx[tag2] -= dx*force;
      grady[tag1] += dy*force;
      grady[tag2] -= dy*force;
    }
  }

  /* quarto contributo, repulsione tra tutti i nodi */
  if (allowrepulsion)
  {
    // compute_repulsive_gradient (contour);
    for (a = contour->vertex; a; a = a->next)
    {
      gradx[a->tag] += k4_coeff*contour->rgradx[a->tag];
      grady[a->tag] += k4_coeff*contour->rgrady[a->tag];
    }
  }

  /* fifth contribution: integral of |k| */
 
  for (line = contour->line; line; line = line->next)
  {
    p = line->a;
    q = line->b;
    /* we need the angles at a and b */
    b = nextp(line, q)->b;
    a = prevp(line, p)->a;
    apx = p->x - a->x;
    apy = p->y - a->y;
    pqx = q->x - p->x;
    pqy = q->y - p->y;
    qbx = b->x - q->x;
    qby = b->y - q->y;
    vec1 = - pqy*apx + pqx*apy;
    vec2 = - qby*pqx + qbx*pqy;
    if (vec1 * vec2 >= 0) continue;   /* no inflection */
    if (SMOOTH_K5 == 0.0)
      sigma = 2.0; else
      sigma = -vec1*vec2/SMOOTH_K5;
    if (sigma > 2.0) sigma = 2.0;
    if (vec2 < 0) sigma *= -1;
    pqsq = pqx*pqx + pqy*pqy;
    gradx[p->tag] += k5_coeff*sigma*pqy/pqsq;
    grady[p->tag] -= k5_coeff*sigma*pqx/pqsq;
    gradx[q->tag] -= k5_coeff*sigma*pqy/pqsq;
    grady[q->tag] += k5_coeff*sigma*pqx/pqsq;
  }

  for (a = contour->vertex; a; a = a->next)
  {
    if (a->type == V_FIXED) gradx[a->tag] = grady[a->tag] = 0.0;
  }
  return;
}

/* this shouldn't be computed at each time step! */

double
compute_repulsive_energy (struct polyline *contour)
{
  struct line *s1, *s2;
  struct vertex *a, *b, *c, *d;
  double s1x, s1y, s2x, s2y, px, py;
  double dd, dsq, f1, f2, f3;
  double energy = 0;

  contour->rentime = contour->time;

  for (s1 = contour->line; s1; s1 = s1->next)
  {
    a = s1->a;
    b = s1->b;
    s1x = b->x - a->x;
    s1y = b->y - a->y;
    for (s2 = s1->next; s2; s2 = s2->next)
    {
      c = s2->a;
      d = s2->b;
      s2x = d->x - c->x;
      s2y = d->y - c->y;
      px = 0.5*(c->x + d->x - a->x - b->x);
      py = 0.5*(c->y + d->y - a->y - b->y);
      dsq = px*px + py*py;
      dd = sqrt (dsq);
      f1 = - s1y*px + s1x*py;
      f2 = - s2y*px + s2x*py;
      f3 = repulsive_field (dd) / dsq;
      energy += fabs(f1*f2*f3);
    }
  }
  contour->renergy = energy;
  return (energy);
}

void
compute_repulsive_gradient (struct polyline *contour)
{
  struct line *s1, *s2;
  struct vertex *a, *b, *c, *d;
  double s1x, s1y, s2x, s2y, px, py;
  double dd, dsq, f1, f2, f3;
  double fpond, e, ds1ex, ds1ey, ds2ex, ds2ey, hdpex, hdpey;
  int taga, tagb, tagc, tagd, i, signf1, signf2;
  double *rgradx, *rgrady;

  contour->rgradtime = contour->time;
  rgradx = contour->rgradx;
  rgrady = contour->rgrady;

  if (contour->rgraddim != contour->numvertices || rgradx == 0 || rgrady == 0)
  {
    if (rgradx) free (rgradx);
    if (rgrady) free (rgrady);
    contour->rgraddim = contour->numvertices;
    rgradx = contour->rgradx = (double *) malloc (contour->rgraddim * sizeof (double));
    rgrady = contour->rgrady = (double *) malloc (contour->rgraddim * sizeof (double));
  }

  for (i = 0; i < contour->rgraddim; i++)
  {
    rgradx[i] = rgrady[i] = 0.0;
  }
  
  for (s1 = contour->line; s1; s1 = s1->next)
  {
    a = s1->a;
    b = s1->b;
    s1x = b->x - a->x;
    s1y = b->y - a->y;
    taga = a->tag;
    tagb = b->tag;
    for (s2 = s1->next; s2; s2 = s2->next)
    {
      c = s2->a;
      d = s2->b;
      s2x = d->x - c->x;
      s2y = d->y - c->y;
      tagc = c->tag;
      tagd = d->tag;
      px = 0.5*(c->x + d->x - a->x - b->x);
      py = 0.5*(c->y + d->y - a->y - b->y);
      dsq = px*px + py*py;
      dd = sqrt (dsq);
      signf1 = signf2 = 1;
      f1 = - s1y*px + s1x*py;
      if (f1 < 0) signf1 = -1;
      f1 *= signf1;
      f2 = - s2y*px + s2x*py;
      if (f2 < 0) signf2 = -1;
      f2 *= signf2;
      f3 = repulsive_field (dd)/dsq;
      fpond = repulsive_force (dd);

      e = f1*f2*f3;
      ds1ex = signf1*f2*f3*py;
      ds1ey = -signf1*f2*f3*px;
      ds2ex = signf2*f1*f3*py;
      ds2ey = -signf2*f1*f3*px;
      hdpex = -signf1*f2*f3*s1y - signf2*f1*f3*s2y - (2*f3 - fpond)*f1*f2*px/dsq;
      hdpey =  signf1*f2*f3*s1x + signf2*f1*f3*s2x - (2*f3 - fpond)*f1*f2*py/dsq;

      hdpex *= 0.5;
      hdpey *= 0.5;

      rgradx[taga] -= ds1ex + hdpex;
      rgrady[taga] -= ds1ey + hdpey;
      rgradx[tagb] += ds1ex - hdpex;
      rgrady[tagb] += ds1ey - hdpey;
      rgradx[tagc] -= ds2ex - hdpex;
      rgrady[tagc] -= ds2ey - hdpey;
      rgradx[tagd] += ds2ex + hdpex;
      rgrady[tagd] += ds2ey + hdpey;
    }
  }
  return;
}

double
repulsive_field (double dist)
{
  if (dist < NODE_SEP)
  {
    return ((NODE_SEP - dist)*(NODE_SEP - dist)/dist);
  }
  return (0.0);
}

/*
 * this is actually the derivative of the field with respect to
 * dist divided by dist
 */

double
repulsive_force (double dist)
{
  double fct;

  if (dist < NODE_SEP)
  {
    fct = (NODE_SEP - dist)/dist;
    return (- fct*(2 + fct)/dist);
  }
  return (0.0);
}

double
get_alpha (struct vertex *a, struct vertex *p, struct vertex *b,
           double *lapt, double *lbpt)
{
  double alfa, apx, apy, pbx, pby;

  apx = p->x - a->x;
  apy = p->y - a->y;
  pbx = b->x - p->x;
  pby = b->y - p->y;

  alfa = atan2 (pby, pbx) - atan2 (apy, apx);
  while (alfa > M_PI) alfa -= 2*M_PI;
  while (alfa < -M_PI) alfa += 2*M_PI;

  *lapt = sqrt (apx*apx + apy*apy);
  *lbpt = sqrt (pbx*pbx + pby*pby);

  return (alfa);
}

void
grad_curv (struct vertex *b, struct vertex *p, double alpha, double lb, double la,
           double *gcxpt, double *gcypt)
{
  double bpx, bpy;

  bpx = b->x - p->x;
  bpy = b->y - p->y;

  *gcxpt = -2*bpy/lb - alpha*bpx/(la + lb);
  *gcypt =  2*bpx/lb - alpha*bpy/(la + lb);

  *gcxpt *= 2*alpha/(la + lb)/lb;
  *gcypt *= 2*alpha/(la + lb)/lb;

//printf ("gcx = %lf, gcy = %lf\n", *gcxpt, *gcypt);
}
