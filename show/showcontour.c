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
#include "energy.h"
#include "../parser.h"

#define ALLOW_NODE_REDEFINE 1
#define ALLOW_NODE_REDEFINE_AT_END 0

#define BUFSIZE 1000
#define REL_H 0.5
#define MAX_H 0.1

#define ME_TRAN 1
#define ME_TOP 2
#define ME_BOT 3
#define ME_CROSS 4
#define ME_NEWROW 5
#define ME_LASTROW 6

#define TOP_LENGTH 0.25

void parseargs (int argc, char *argv[]);
void test_contour (struct polyline *contour);
double normsq (double *x, int dim);
void reorder_node_ptr (struct polyline *contour);
int settags (struct polyline *contour);
double getlen (struct line *line);
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

struct polyline *contour;
static double time;
static double tau;
static double taurn;
double taurep;
static double curenergy = 0.0;

static double timerrn;
double timerrep;

static int test = 0;
static char *grident;

int
main (int argc, char *argv[])
{
  struct line *l;
  struct vertex *v;
  int tok, count, vertexnum;

  grident = grinit(&argc, argv);
  energyinit ();
  parseargs (argc, argv);
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
parseargs (int argc, char *argv[])
{
  int iarg;

  for (iarg = 1; iarg < argc; iarg++)
  {
    if (*argv[iarg] == '-')    /* this is an option */
    {
      if (strcmp (argv[iarg], "--grident") == 0)
      {
        printf ("%s\n", grident);
        exit (0);
      }
      if (strcmp (argv[iarg], "--test") == 0) {
        test = 1;
      } else if (strcmp (argv[iarg], "--k5") == 0) {
        k5_coeff = atof (argv[++iarg]);
      } else {
        printf ("Invalid option: %s\n", argv[iarg]);
        exit (1);
      }
      continue;
    }
    printf ("Invalid argument: %s\n", argv[iarg]);
    exit (1);
  }
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
