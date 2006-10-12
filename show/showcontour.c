/*
 * compile with "cc -o showcontour showcontour.c glutcontour.o parser.o -lglut"
 *
 * requires freeglut-devel package (on Fedora core)
 *
 * usage:  "./showcontour <example.morse"
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "showcontour.h"
#include "../parser.h"

//#define CHECK_GRADIENT 1

#define K1_COEFF 1.0     /* perimeter */
#define K2_COEFF 0.1     /* k^2 */
#define K3_COEFF 0.7     /* cross+cusp repulsion */
#define K4_COEFF 1.0     /* arc repulsion */
#define K5_COEFF 1.0     /* |k| */

#define ALLOW_REPULSION 1
#define ALLOW_NODE_REDEFINE 1
#define ALLOW_NODE_REDEFINE_AT_END 0

#define NODE_SEP 0.4
#define ORTO_HEAVINESS 2.0
#define CUSP_PRECISION 0.8           /* 1 means perfect cusp (180 degrees) */
/*
 * a value too near 1 hase the effect of shooting out (sometimes) cusps towards
 * infinity, probably a side effect of the addition of nodes for elongating
 * segments.  The energy associated with cusps should be modified in order
 * to avoid such tendency
 */

#define BUFSIZE 1000
#define REL_H 0.5
#define MAX_H 0.1

#define ME_TRAN 1
#define ME_TOP 2
#define ME_BOT 3
#define ME_CROSS 4
#define ME_NEWROW 5
#define ME_LASTROW 6

#define V_REGULAR 1
#define V_CUSP 2
#define V_CROSS 3
#define V_FIXED 4

#define TOP_LENGTH 0.25

void test_contour (struct polyline *contour);
void compute_gradient (struct polyline *contour);
double compute_energy (struct polyline *contour);
void compute_repulsive_gradient (struct polyline *contour);
double compute_repulsive_energy (struct polyline *contour);
double repulsive_field (double dist);
double repulsive_force (double dist);
double normsq (double *x, int dim);
void reorder_node_ptr (struct polyline *contour);
int settags (struct polyline *contour);
double getlen (struct line *line);
struct line *nextp (struct line *l, struct vertex *p);
struct line *prevp (struct line *l, struct vertex *p);
double get_alpha (struct vertex *a, struct vertex *p, struct vertex *b, 
       double *lapt, double *lbpt);
void grad_curv (struct vertex *b, struct vertex *p, double alpha, double lb, double la,
           double *gcxpt, double *gcypt);
void discretizepolyline (struct polyline *contour);
void specnodesinit (struct polyline *contour);
void redistributenodes (struct polyline *contour);
struct line *splitline (struct polyline *contour, struct line *line, double f);
void getmorseevent (struct morseevent *mev);
void getarcinfo (struct morseevent *morseevent);
void getoricusps (int *oript, int *cuspspt);
struct polyline *buildpolyline (void);
struct vertex *newvertex (struct polyline *contour, 
          double x, double y, int type);
struct line *newline (struct polyline *contour, 
                      struct vertex *a, 
                      struct vertex *b);
int inherit_orientation (struct line *line);
void insertcusps (struct polyline *contour);
int insert_cusps_on_arc (struct line *l);
void removenode (struct polyline *contour, struct vertex *v);
void removeline (struct polyline *contour, struct line *l);
#ifdef CHECK_GRADIENT
void check_gradient (struct polyline *contour);
#endif

struct polyline *contour;
static double time;
static double tau;
static double taurn;
static double taurep;
static double curenergy = 0.0;
static double k1_coeff = K1_COEFF;
static double k2_coeff = K2_COEFF;
static double k3_coeff = K3_COEFF;
static double k4_coeff = K4_COEFF;
static double k5_coeff = K5_COEFF;

static double timerrn;
static double timerrep;
static int allowrepulsion = ALLOW_REPULSION;

static int test = 0;

int
main (int argc, char *argv[])
{
  struct line *l;
  struct vertex *v;
  int tok, count, vertexnum, iarg;
  char *grident;

  grident = grinit(&argc, argv);
  for (iarg = 1; iarg < argc; iarg++)
  {
    if (*argv[iarg] == '-')    /* this is an option */
    {
      if (strcmp (argv[iarg], "--grident") == 0)
      {
        printf ("%s\n", grident);
        return (0);
      }
      if (strcmp (argv[iarg], "--test") == 0) {
        test = 1;
      } else {
        printf ("Invalid option: %s\n", argv[iarg]);
        return (1);
      }
      continue;
    }
    printf ("Invalid argument: %s\n", argv[iarg]);
    return (1);
  }
  if (test == 0) allowrepulsion = 0;
  tok = gettoken (stdin);
  if (tok != TOK_MORSE) return (0);
  tok = gettoken (stdin);
  if (tok != TOK_LBRACE) return (0);

  contour = buildpolyline ();

  tok = gettoken (stdin);
  if (tok != TOK_RBRACE) exit (1);

  reorder_node_ptr (contour);
  discretizepolyline (contour);
  insertcusps (contour);
  specnodesinit (contour);

//test_contour (contour);
  vertexnum = settags (contour);
  contour->gradx = (double *) malloc (vertexnum * sizeof (double));
  contour->grady = (double *) malloc (vertexnum * sizeof (double));

  count = 0;
  for (v = contour->vertex; v; v = v->next) count++;
//  printf ("ci sono %d vertici\n", count);
  count = 0;
  for (l = contour->line; l; l = l->next)
  {
    count++;
//    printf ("lung %lf\n", getlen (l));
  }
//  printf ("ci sono %d archi\n", count);

  tau = contour->h * contour->h;
  taurn = taurep = tau;
  tau = tau*tau;
//  printf ("tau = %lf, taurn = %lf, taurep = %lf\n", tau, taurn, taurep);
  time = 0.0;
  timerrn = time + taurn;
  timerrep = -1.0;          /* must compute immediately */

  //evolve (contour, incrtime);

  grmain();
  return (0);
}

void
test_contour (struct polyline *contour)
{
  struct vertex *a, *b, *c;
  struct line *l1, *l2;

  contour->vertex = 0;
  contour->line = 0;

  a = newvertex (contour, -0.5, -0.5, V_FIXED);
  b = newvertex (contour, 0.5, -0.5, V_REGULAR);
  c = newvertex (contour, 0.5, 0.5, V_FIXED);
  l1 = newline (contour, a, b);
  l2 = newline (contour, b, c);
  a->line[0] = 0;
  a->line[1] = l1;
  b->line[0] = l1;
  b->line[1] = l2;
  c->line[0] = l2;
  c->line[1] = 0;

  contour->h = 0.001;
}

void
insertcusps (struct polyline *contour)
{
  struct vertex *p;
  struct line *l;
  int i, iss1;

  /* first deal with the non-s1 arcs */

  for (p = contour->vertex; p; p = p->next)
  {
    if (p->type != V_CROSS) continue;
    for (i = 0; i < 4; i++)
    {
      l = p->line[i];
      if (l->a != p) continue;   /* we want a 'leaving' arc */
      if (l->cusps <= 0) continue;
      iss1 = insert_cusps_on_arc (l);
      assert (iss1 == 0);
    }
  }

  /* now the s1's */

  for (p = contour->vertex; p; p = p->next)
  {
    if (p->type == V_CROSS) continue;
    l = p->line[1];
    if (l-> cusps <= 0) continue;
    iss1 = insert_cusps_on_arc (l);
    assert (iss1);
  }

  /* now check that no cusps remain to be places */

  for (l = contour->line; l; l = l->next)
  {
    assert (l->cusps == 0);
  }
}

int
insert_cusps_on_arc (struct line *l)
{
  struct line *wl;
  struct vertex *p;
  int cusps, iss1, count, ccount, nnodes;

  iss1 = 0;
  count = 0;                /* count number of nodes */
  cusps = l->cusps;
  wl = l;
  do {
    wl->cusps = 0;
    count++;
    if ((p = wl->b)->type == V_CROSS) break;
  } while (wl = p->line[1], wl != l);

  if (p->type != V_CROSS) iss1 = 1;

  nnodes = count/(cusps + 1) + 1;

  count = ccount = 0;
  wl = l;
  do {
    count++;
    if ((p = wl->b)->type == V_CROSS) break;
    if (count >= nnodes)
    {
      count = 0;
      ccount++;
      assert (p->type == V_REGULAR);
      p->type = V_CUSP;
    }
  } while (wl = p->line[1], wl != l);

  assert (ccount == cusps);
//printf ("inserted %d cusps\n", cusps);
  return (iss1);
}

int
settags (struct polyline *contour)
{
  struct vertex *v;
  int count = 0;

  for (v = contour->vertex; v; v = v->next)
    v->tag = count++;

  contour->numvertices = count;
  return (count);
}

static struct timeval timer;

void
settimer (double incrtime)
{
  int seconds, useconds;

  gettimeofday(&timer, 0);

  seconds = (int) incrtime;
  useconds = (int) 1000000*(incrtime - seconds);

  timer.tv_usec += useconds;
  timer.tv_sec += seconds;
  if (timer.tv_usec > 1000000)
  {
    timer.tv_usec -= 1000000;
    timer.tv_sec++;
  }
}

int
checktimer (void)
{
  struct timeval now;

  gettimeofday (&now, 0);
  return (timercmp (&now, &timer, < ));
}

static int allownodered = ALLOW_NODE_REDEFINE;
static int allownoderedatend = ALLOW_NODE_REDEFINE_AT_END;

void
tryredistributenodes (struct polyline *contour)
{
  static int count = 0;

  if (time < timerrn) return;
  count++;

  //printf ("would redistribute nodes: %d...\n", count);
  if (allownodered) redistributenodes (contour);
  timerrn += taurn;
}

void
tryrepulsiveenergy (struct polyline *contour)
{
  static int count = 0;

  if (time < timerrep) return;
  count++;

  //printf ("would compute repulsive energy: %d...\n", count);
  compute_repulsive_energy (contour);
  compute_repulsive_gradient (contour);
#ifdef CHECK_GRADIENT
  check_gradient (contour);
#endif
  timerrep += taurep;
}

double
evolve (struct polyline *contour, double incrtime)
{
  double newenergy, decrease, predicted_decrease;
  struct vertex *v;
  int gradient_is_ok;
  int timesteps = 0;
  int timebsteps = 0;
  int allowbackstep = 0;
  int allowstepcontrol = 0;

  settimer (incrtime);
  gradient_is_ok = 0;
  while (checktimer())
  {
    tryredistributenodes (contour);
    tryrepulsiveenergy (contour);
    time += tau;
    timesteps++;
    if (curenergy == 0.0) curenergy = compute_energy (contour);
    if (! gradient_is_ok) compute_gradient (contour);

    for (v = contour->vertex; v; v = v->next)
    {
      v->x -= tau*contour->gradx[v->tag];
      v->y -= tau*contour->grady[v->tag];
    }
    gradient_is_ok = 0;
    newenergy = compute_energy (contour);
    decrease = curenergy - newenergy;
 /* perhaps something is wrong, since this happens too often */
    if (allowbackstep && decrease < 0.0)
    {
//      fprintf (stderr, "WARNING: energy is increasing!\n");
      /* torno indietro di uno step e diminuisco tau */
      for (v = contour->vertex; v; v = v->next)
      {
        v->x += tau*contour->gradx[v->tag];
        v->y += tau*contour->grady[v->tag];
      }
      time -= tau;
      timesteps--;
      timebsteps++;
      tau /= 4;
      gradient_is_ok = 1;
      continue;
    }
    curenergy = newenergy;
    predicted_decrease = tau * (
      normsq (contour->gradx, contour->numvertices) +
      normsq (contour->grady, contour->numvertices) );
    //printf ("predicted decrease/actual decrease = %lf/%lf = %lf\n",
    //   predicted_decrease, decrease, predicted_decrease/decrease);

    if (allowstepcontrol)
    {
 /* automatic step-size control does not work yet! */
      if (fabs(predicted_decrease - decrease) < 0.05*predicted_decrease)
      {
//        printf ("good prediction, increasing tau\n");
        tau *= 1.1;
      }

      if (fabs(predicted_decrease - decrease) > 0.2*predicted_decrease)
      {
//        printf ("bad prediction, decreasing tau\n");
        tau /= 2;
      }
    }
  }
//  printf ("timesteps: %d, backsteps: %d, time = %lf, energy = %lf\n", 
//          timesteps, timebsteps, time, curenergy);
  if (allownoderedatend) redistributenodes (contour);
  return (time);
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

double
normsq (double *x, int dim)
{
  double norm = 0.0;
  int i;

  for (i = 0; i < dim; i++)
  {
    norm += x[i]*x[i];
  }
  return (norm);
}

static double *rgradx = 0;
static double *rgrady = 0;
static double renergy = 0.0;
static int rgraddim = 0;

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
    energy4 = k4_coeff*renergy;
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
  int i, j, tag1, tag2, sigma;
  double dx, dy, len, xel, yel, distsq, dist, force;
  double alpha, gcx, gcy, la, lb;
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
      gradx[a->tag] += k4_coeff*rgradx[a->tag];
      grady[a->tag] += k4_coeff*rgrady[a->tag];
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
    sigma = 2;
    if (vec2 < 0) sigma = -2;
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
  renergy = energy;
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

  if (rgraddim != contour->numvertices || rgradx == 0 || rgrady == 0)
  {
    if (rgradx) free (rgradx);
    if (rgrady) free (rgrady);
    rgraddim = contour->numvertices;
    rgradx = (double *) malloc (rgraddim * sizeof (double));
    rgrady = (double *) malloc (rgraddim * sizeof (double));
  }

  for (i = 0; i < rgraddim; i++)
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

struct line *
nextp (struct line *l, struct vertex *p)
{
  int i;

  if (p->type != V_CROSS) return (p->line[1]);
  for (i = 0; i < 4; i++)
  {
    if (p->line[i] == l) return (p->line[3-i]);
  }
  assert (0);
}

struct line *
prevp (struct line *l, struct vertex *p)
{
  int i;

  if (p->type != V_CROSS) return (p->line[0]);
  for (i = 0; i < 4; i++)
  {
    if (p->line[i] == l) return (p->line[3-i]);
  }
  assert (0);
}

void
reorder_node_ptr (struct polyline *contour)
{
  struct vertex *v;
  struct line *ltemp;

  for (v = contour->vertex; v; v = v->next)
  {
    if (v->type == V_CROSS) continue;   /* non tocco l'ordine */
    if (v->line[0]->b != v)
    {
      ltemp = v->line[0];
      v->line[0] = v->line[1];
      v->line[1] = ltemp;
    }
  }
}

double
getlen (struct line *line)
{
  double dx, dy;

  dx = line->a->x - line->b->x;
  dy = line->a->y - line->b->y;

  return (sqrt (dx*dx + dy*dy));
}

void
discretizepolyline (struct polyline *contour)
{
  double newh;
  struct vertex *a, *b;
  struct line *line, *nline;
  double diffx, diffy, len, dt;
  int i, numsub;

  newh = REL_H * contour->h;
  if (newh > MAX_H) newh = MAX_H;
  contour->h = newh;
  for (line = contour->line; line; line = line->next)
  {
    a = line->a;
    b = line->b;
    diffx = a->x - b->x;
    diffy = a->y - b->y;
    len = sqrt (diffx*diffx + diffy*diffy);

    numsub = len/newh + 1;
    dt = len/numsub;
    for (i = 1; i < numsub; i++)
    {
      nline = splitline (contour, line, dt/len);
      len -= dt;
      line = nline;
    }
  }
}

void
specnodesinit (struct polyline *contour)
{
  struct vertex *v;
  int i;

  contour->specnodesnum = 0;

  for (v = contour->vertex; v; v = v->next)
  {
    if (v->type != V_REGULAR) contour->specnodesnum++;
  }
  contour->specnodes = (struct vertex **) malloc 
    (contour->specnodesnum*sizeof (struct vertex *));
  for (i = 0, v = contour->vertex; v; v = v->next)
  {
    if (v->type != V_REGULAR)
    {
      contour->specnodes[i++] = v;
    } 
  }
}

void
redistributenodes (struct polyline *contour)
{
  struct line *line, *nline, *aline;
  struct vertex *a, *b, *v, *next;
  double diffx, diffy, len, alen, h, dt;
  int numsub, i, vertexnum, forward, backidx;
  int addednodes = 0, removednodes = 0;

  /* contour->h is the target step size */

  h = contour->h;
  for (line = contour->line; line; line = line->next)
  {
    a = line->a;
    b = line->b;
    diffx = a->x - b->x;
    diffy = a->y - b->y;
    len = sqrt (diffx*diffx + diffy*diffy);

    if (len < 0.5*h)
    {   /* can I derefine? */
      v = line->b;
      forward = 1;
      if (v->type != V_REGULAR) {v = line->a; forward = 0;}
      if (v->type != V_REGULAR) continue;
      aline = v->line[forward];
      assert (aline != line);
      next = aline->a;
      if (forward) next = aline->b;
      assert (next != v);
      backidx = 1 - forward;
      if (next->type == V_CROSS)
      {
        for (backidx = 0; backidx < 4; backidx++)
          if (next->line[backidx] == aline) break;
      }
      diffx = v->x - next->x;
      diffy = v->y - next->y;
      alen = sqrt (diffx*diffx + diffy*diffy);
      if (alen + len < 0.9*h)   /* can derefine */
      { /* aline will be removed! */
        removednodes++;
        assert (next->line[backidx] == aline);
        if (forward) line->b = next;
          else line->a = next;
        next->line[backidx] = line;
        removenode (contour, v);
        removeline (contour, aline);
      }
    } else {
      /* possibly refine, only if len >= 2*h */
      numsub = len/h;
      dt = len/numsub;
      for (i = 1; i < numsub; i++)
      {
        addednodes++;
        nline = splitline (contour, line, dt/len);
        len -= dt;
        line = nline;
      }
    }
  }
  vertexnum = settags (contour);
  free (contour->gradx);
  free (contour->grady);
  contour->gradx = (double *) malloc (vertexnum * sizeof (double));
  contour->grady = (double *) malloc (vertexnum * sizeof (double));
  curenergy = 0.0;
  // if (addednodes) printf ("added %d nodes\n", addednodes);
  // if (removednodes) printf ("removed %d nodes\n", removednodes);
}

struct line *
splitline (struct polyline *contour, struct line *line, double f)
{
  struct vertex *a, *b, *p;
  struct line *nline;
  double x, y;
  int nl, i;

  a = line->a;
  b = line->b;

  x = (1 - f)*a->x + f*b->x;
  y = (1 - f)*a->y + f*b->y;

  p = newvertex (contour, x, y, V_REGULAR);
  nline = (struct line *) malloc (sizeof (struct line));
  nline->orientation = line->orientation;
  nline->cusps = line->cusps;
  nline->next = line->next;
  nline->a = p;
  nline->b = b;
  nl = 2;
  if (b->type == V_CROSS) nl = 4;
  for (i = 0; i < nl; i++)
  {
    if (b->line[i] == line) b->line[i] = nline;
  }

  line->next = nline;
  line->b = p;

  p->line[0] = line;
  p->line[1] = nline;

  return (nline);
}

void
getmorseevent (struct morseevent *morseevent)
{
  int tok;

  morseevent->ori = morseevent->ori2 = 0;
  morseevent->cusps = morseevent->cusps2 = 0;

  /* ho gia letto la graffa aperta */
  tok = gettokens (stdin);

  switch (tok)
  {
    case TOK_SEMICOLON:
      tok = gettokens (stdin);
      if (tok == TOK_RBRACE)
      {
        morseevent->type = ME_LASTROW;
      } else {
        morseevent->type = ME_NEWROW;
      }
      ungettoken (tok);
      break;

    case KEY_HAT:
    case KEY_A:
      morseevent->type = ME_TOP;
      getarcinfo (morseevent);
      break;

    case KEY_U:
    case KEY_V:
    case KEY_UNDERSCORE:
      morseevent->type = ME_BOT;
      getarcinfo (morseevent);
      break;

    case KEY_SLASH:
    case KEY_BSLASH:
    case KEY_BACKQUOTE:
    case TOK_LPAREN:
    case TOK_RPAREN:
    case KEY_PIPE:
    case KEY_I:
      morseevent->type = ME_TRAN;
      getarcinfo (morseevent);
      break;

    case KEY_X:
      morseevent->type = ME_CROSS;
      getarcinfo (morseevent);
      break;
  }
  return;
}

void
getarcinfo (struct morseevent *morseevent)
{
  getoricusps (&morseevent->ori, &morseevent->cusps);
  if (morseevent->type == ME_CROSS)
    getoricusps (&morseevent->ori2, &morseevent->cusps2);
}

void
getoricusps (int *oript, int *cuspspt)
{
  int tok, prevd;
  int require_rbr = 1;
  int depthind = 0;

  *oript = *cuspspt = 0;
  tok = gettokens (stdin);
  if (tok == ISNUMBER || tok == KEY_LEFT ||
      tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN)
  {
    ungettoken (tok);
    tok = TOK_LBRACKET;
    require_rbr = 0;
  }
  if (tok != TOK_LBRACKET)
  {
    ungettoken (tok);
    return;
  }
  tok = gettokens (stdin);
  if (tok == TOK_RBRACKET) return;
  if (tok == KEY_LEFT || tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN)
  {
    if (tok == KEY_LEFT || tok == KEY_DOWN) *oript = 1;
    if (tok == KEY_RIGHT || tok == KEY_UP) *oript = -1;
    tok = gettokens (stdin);
  }
  if (tok == TOK_COMMA || tok == ISNUMBER)
  {
    if (tok == ISNUMBER) ungettoken (tok);
    prevd = 0;
    while ((tok = gettokens (stdin)) == ISNUMBER ||
            tok == TOK_PLUS || tok == TOK_MINUS)
    {
      switch (tok)
      {
        case ISNUMBER:
        prevd = gettokennumber ();
        depthind++;
        break;
        case TOK_PLUS:
        ++prevd;
        depthind++;
        break;
        case TOK_MINUS:
        --prevd;
        depthind++;
        break;
      }
    }
  }
  if (depthind > 0) *cuspspt = depthind - 1;
  if (require_rbr == 0)
  {
    ungettoken (tok);
    tok = TOK_RBRACKET;
  }
  if (tok != TOK_RBRACKET)
  {
    fprintf (stderr, "Error: right paren expected: %d\n", tok);
    return;
  }
}

struct line *
newline (struct polyline *contour, struct vertex *a, struct vertex *b)
{
  struct line *line;

  line = (struct line *) malloc (sizeof (struct line));
  line->a = a;
  line->b = b;
  line->orientation = 0;
  line->next = contour->line;
  contour->line = line;
  return (line);
}

struct vertex *
newvertex (struct polyline *contour, double x, double y, int type)
{
  struct vertex *v;
  int numarcs, i;

  numarcs = 2;
  if (type == V_CROSS) numarcs = 4;
  v = (struct vertex *) malloc (sizeof (struct vertex) 
                         + numarcs*sizeof (struct line *));
  v->x = x;
  v->y = y;
  v->next = contour->vertex;
  v->type = type;
  for (i = 0; i < numarcs; i++) v->line[i] = 0;
  contour->vertex = v;
  return (v);
}

void
removeline (struct polyline *contour, struct line *line)
{
  struct line *l;

  if (contour->line == line)
  {
    contour->line = line->next;
    free (line);
    return;
  } else {
    for (l = contour->line; l; l = l->next)
    {
      if (l->next == line)
      {
        l->next = line->next;
        free (line);
        return;
      }
    }
  }
  assert (0);
}

void
removenode (struct polyline *contour, struct vertex *node)
{
  struct vertex *n;

  if (contour->vertex == node)
  {
    contour->vertex = node->next;
    free (node);
    return;
  } else {
    for (n = contour->vertex; n; n = n->next)
    {
      if (n->next == node)
      {
        n->next = node->next;
        free (node);
        return;
      }
    }
  }
  assert (0);
}

struct polyline *
buildpolyline (void)
{
  double dx, dy, x, y;
  struct morseevent morseevent;
  struct vertex *danglingnodes[BUFSIZE];
  struct vertex *prevdanglingnodes[BUFSIZE];
  struct vertex *v1, *v2, *v3, *v4, *v5, *v;
  struct line *line;
  struct polyline *contour;
  int goon, i, k, prevdangnodes, dangind, prevdangind;
  int numrows = 0, maxrowlen = 0, oriented;

  contour = (struct polyline *) malloc (sizeof (struct polyline));
  contour->vertex = 0;
  contour->line = 0;
  contour->gradx = contour->grady = 0;
  contour->specnodes = 0;
  contour->specnodesnum = 0;

  dx = dy = 1.0;    /* aggiusto alla fine */
  i = 0;            /* conta gli eventi su ciascuna riga */
  numrows = 0;            /* conta il numero di righe */
  prevdangnodes = prevdangind = dangind = 0;
  y = x = 0.0;
  goon = 1;
  while (goon)
  {
    i++;
    x += dx;
    getmorseevent (&morseevent);
    switch (morseevent.type)
    {
      case ME_LASTROW:
        goon = 0;
        i--;
        if (i > maxrowlen) maxrowlen = i;
        if (prevdangnodes != prevdangind)
        {
          fprintf (stderr, "dangling nodes: prev = %d, current = %d\n",
                  prevdangnodes, dangind);
        }

        if (dangind != 0) fprintf (stderr, "Dangling nodes remain...\n");
        break;

      case ME_NEWROW:
        y += dy;
        x = 0.0;
        numrows++; i--;
        if (i > maxrowlen) maxrowlen = i;
        i = 0;
        if (prevdangnodes != prevdangind)
        {
          fprintf (stderr, "dangling nodes: prev = %d, current = %d\n",
                  prevdangnodes, dangind);
        }
        for (k = 0; k < dangind; k++)
        {
          prevdanglingnodes[k] = danglingnodes[k];
        }
        prevdangnodes = dangind;
        prevdangind = dangind = 0;
        break;

      case ME_TOP:
        v1 = newvertex (contour, x - TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
        v2 = newvertex (contour, x + TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
        line = newline (contour, v2, v1);
        line->orientation = morseevent.ori;
        line->cusps = morseevent.cusps;
        v1->line[0] = line;
        v2->line[0] = line;
        danglingnodes[dangind++] = v1;
        danglingnodes[dangind++] = v2;
        break;

      case ME_BOT:
        v1 = newvertex (contour, x - TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
        v2 = newvertex (contour, x + TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
        line = newline (contour, v2, v1);
        line->orientation = morseevent.ori;
        line->cusps = morseevent.cusps;
        v1->line[0] = line;
        v2->line[0] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v1);
        v1->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v2);
        v2->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        break;

      case ME_CROSS:
        v1 = newvertex (contour, x, y, V_CROSS);
        v2 = newvertex (contour, x - TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
        v3 = newvertex (contour, x + TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
        v4 = newvertex (contour, x - TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
        v5 = newvertex (contour, x + TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
        line = newline (contour, v2, v1);
        v1->line[0] = line;
        v2->line[0] = line;
        line = newline (contour, v3, v1);
        v1->line[1] = line;
        v3->line[0] = line;
        line = newline (contour, v1, v4);
        line->orientation = morseevent.ori;
        line->cusps = morseevent.cusps;
        v1->line[2] = line;
        v4->line[0] = line;
        line = newline (contour, v1, v5);
        line->orientation = morseevent.ori2;
        line->cusps = morseevent.cusps2;
        v1->line[3] = line;
        v5->line[0] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v2);
        v2->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v3);
        v3->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        danglingnodes[dangind++] = v4;
        danglingnodes[dangind++] = v5;
        break;

      case ME_TRAN:
        v1 = newvertex (contour, x, y - TOP_LENGTH*dy, V_REGULAR);
        v2 = newvertex (contour, x, y + TOP_LENGTH*dy, V_REGULAR);
        line = newline (contour, v1, v2);
        line->orientation = morseevent.ori;
        line->cusps = morseevent.cusps;
        v1->line[0] = line;
        v2->line[0] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v1);
        v1->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        danglingnodes[dangind++] = v2;
        break;

      default:
        fprintf (stderr, "Invalid morse event\n");
    }
  }

  /* rinormalizzazione ascisse e ordinate per stare in [-1,1]x[-1,1] */
  dx = 2.0/(maxrowlen + 1);
  dy = 2.0/numrows;
  contour->h = dx;
  if (dy < dx) contour->h = dy;

  for (v = contour->vertex; v; v = v->next)
  {
    v->x = dx*v->x - 1.0;
    v->y = 1.0 - dy*v->y;
    assert (v->line[0]);
    assert (v->line[1]);
    if (v->type == V_CROSS) assert (v->line[2] && v->line[3]);
  }

  /* orientazione degli archi */

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (line = contour->line; line; line = line->next)
    {
      if (line->orientation != 0 && inherit_orientation (line)) goon = 1;
    }
  }

  oriented = 1;
  for (line = contour->line; line; line = line->next)
  {
    if (line->orientation == 0)
    {
      fprintf (stderr, "nonoriented arc\n");
      oriented = 0;
    }
  }

  if (oriented == 0) fprintf (stderr, "Cannot orient arc\n");

  return (contour);
}

/* inherit orientation to directly adjacent arcs */

int
inherit_orientation (struct line *line)
{
  struct line *wline;
  struct vertex *v, *vtemp;
  int goon;

  if (line->orientation < 0)
  {
    vtemp = line->a;
    line->a = line->b;
    line->b = vtemp;
    line->orientation *= -1;
  }

  goon = 0;
  v = line->a;  /* going back */
  if (v->type != V_CROSS)
  {
    wline = v->line[0];
    if (wline == line) wline = v->line[1];
    assert (wline && wline != line);
    if (wline->a == v)
    {
      vtemp = wline->a;
      wline->a = wline->b;
      wline->b = vtemp;
      wline->orientation *= -1;
    }
    if (wline->orientation < 0) fprintf (stderr, "Conflicting orientations\n");
    if (wline->orientation == 0)
    {
      goon = 1;
      wline->orientation = 1;
      wline->cusps = line->cusps;
    }
  }

  v = line->b;  /* going forward */

  if (v->type != V_CROSS)
  {
    wline = v->line[0];
    if (wline == line) wline = v->line[1];
    if (wline->b == v)
    {
      vtemp = wline->a;
      wline->a = wline->b;
      wline->b = vtemp;
      wline->orientation *= -1;
    }
    if (wline->orientation < 0) fprintf (stderr, "Conflicting orientations\n");
    if (wline->orientation == 0)
    {
      goon = 1;
      wline->orientation = 1;
      wline->cusps = line->cusps;
    }
  }

  return (goon);
}
