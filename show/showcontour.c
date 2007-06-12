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
#include "morseevent.h"
#include "doptimize.h"
#include "energy.h"
#include "xfigexport.h"
#include "../parser.h"

#define ALLOW_NODE_REDEFINE 1
#define ALLOW_BACKSTEP 0
#define ALLOW_STEPCONTROL 1

//#define KICK_OUT_TIME 50.0    /* k^2 only after this time, relative to dx^2 */
#define KICK_OUT_TIME 0.0

#define BUFSIZE 1000
#define REL_H 0.5
//#define REL_H 0.25
#define MAX_H 0.1

#define TOP_LENGTH 0.25
#define RENORMALIZE 1

void kick_in (struct polyline *contour);
void kick_out (struct polyline *contour);
void activate_timer (int event, double time);
FILE *parseargs (int argc, char *argv[]);
void test_contour (struct polyline *contour);
void reorder_node_ptr (struct polyline *contour);
int settags (struct polyline *contour);
double getlen (struct line *line);
void renormalize (struct polyline *contour);
void discretizepolyline (struct polyline *contour);
void specnodesinit (struct polyline *contour);
struct line *splitline (struct polyline *contour, struct line *line, double f);
struct polyline *buildpolyline (FILE *filein);
struct vertex *newvertex (struct polyline *contour, 
          double x, double y, int type);
struct line *newline (struct polyline *contour, 
                      struct vertex *a, 
                      struct vertex *b);
int inherit_orientation (struct line *line);
int inherit_earc (struct line *line);
void setarcinfo (struct polyline *contour);
int setarcinfo1 (struct line *l);
int enforce_crossover_info (struct polyline *contour);
struct earc *mergearcinfo (struct earc *arc1, struct earc *arc2);
void check_contour (struct polyline *contour);
void insertcusps (struct polyline *contour);
int insert_cusps_on_arc (struct line *l);

struct polyline *contour;
//static double time;
static double tau;
static double taurn;
static double taurep;
static double curenergy = 0.0;

//static double timerrn;
//double timerrep;
//static double timerkickout;

char *xfigproblem = 0;
char *title = 0;
static double skiptime = 0.0;
static double skiprtime = 0.0;
static int test = 0;
static int dodoptimize = 1;
#ifdef HAVE_GTK
static int gr_id = GO_GTK;
#else
#ifdef HAVE_GLUT
static int gr_id = GO_GLUT;
#else
static int gr_id = GO_NULL;
#endif
#endif

//static char *grident;

int
main (int argc, char *argv[])
{
  struct line *l;
  struct vertex *v;
  int tok, count, vertexnum;
  double kick_out_time;
  FILE *filein;
  char *grident;

  energyinit ();
  filein = parseargs (argc, argv);
  if (filein == 0) filein = stdin;
  if (test == 0) allowrepulsion = 0;
  grident = grinit(&argc, argv);

  tok = gettoken (filein);
  if (tok != TOK_MORSE) return (0);
  tok = gettoken (filein);
  if (tok != TOK_LBRACE) return (0);

  contour = buildpolyline (filein);

  tok = gettoken (filein);
  if (tok != TOK_RBRACE) exit (1);

  reorder_node_ptr (contour);
  if (dodoptimize) doptimize (contour);

#ifdef RENORMALIZE
  renormalize (contour);
#endif
  discretizepolyline (contour);
  setarcinfo (contour);        /* initialize arc structure */
  enforce_crossover_info (contour);

  check_contour (contour);

  insertcusps (contour);       /* also sets the d value of each segment */
  specnodesinit (contour);
  init_rarc (contour);

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

  contour->time = 0.0;
  tau = contour->h * contour->h;
  taurn = taurep = tau;
  tau = tau*tau;
  if (KICK_OUT_TIME > 0.0)
  {
    kick_out_time = KICK_OUT_TIME*contour->baseh*contour->baseh;
    printf ("kick_in, out at time %lf\n", kick_out_time);
    activate_timer (EVENT_KICKOUT, kick_out_time);
    kick_in (contour);
  }
//  printf ("tau = %lf, taurn = %lf, taurep = %lf\n", tau, taurn, taurep);
  activate_timer (EVENT_REDISTRIBUTENODES, contour->time + taurn);
  activate_timer (EVENT_REPULSIVEENERGY, -1.0);    /* must compute immediately */

  if (skiptime > 0.0 || skiprtime > 0.0)
  {
    if (skiptime <= 0.0) skiptime = 100.0;
    if (skiprtime <= 0.0) skiprtime = 3600.0;
    evolve (contour, skiprtime, skiptime);
  }
  grmain();
  return (0);
}

struct timerevent *timerlist = 0;
struct timerevent dummytimer;

void
activate_timer (int event, double time)
{
  struct timerevent *newtimer;

  if (timerlist == 0)
  {
    timerlist = &dummytimer;
    timerlist->next = 0;
  }
  newtimer = (struct timerevent *) malloc (sizeof (struct timerevent));
  newtimer->next = timerlist->next;
  timerlist->next = newtimer;
  newtimer->event = event;
  newtimer->time = time;
}

static double k2_coeff_saved;
static double tau_saved;
static double h_saved;
static int kicked_in = 0;
static int immediate_exit = 0;

void
kick_in (struct polyline *contour)
{
  if (kicked_in) 
  {
    //fprintf (stderr, "Warning: already kicked in\n");
    return;
  }
  kicked_in = 1;
  k2_coeff_saved = k2_coeff;
  k2_coeff = 0;
  tau_saved = tau;
  tau = sqrt (tau);
  h_saved = contour->h;
  contour->h /= 2;
}

void
kick_out (struct polyline *contour)
{
  if (! kicked_in)
  {
    //fprintf (stderr, "Warning: already kicked out\n");
    return;
  }
printf ("kick_out\n");
  kicked_in = 0;
  k2_coeff = k2_coeff_saved;
  tau = tau_saved;
  contour->h = h_saved;
}

FILE *
parseargs (int argc, char *argv[])
{
  int iarg, gr_lid;
  FILE *filein = 0;

  grparser (&argc, argv);
  for (iarg = 1; iarg < argc; iarg++)
  {
    if (*argv[iarg] == '-')    /* this is an option */
    {
      if (strcmp (argv[iarg], "--grident") == 0)
      {
        switch (gr_id)
        {
          case GO_NULL:
          printf ("%s\n", "null");
          break;

          case GO_GTK:
          printf ("%s\n", "gtk");
          break;

          case GO_GLUT:
          printf ("%s\n", "glut");
          break;

          case GO_XFIG:
          printf ("%s\n", "xfig");
          break;

          default:
          printf ("unknown\n");
          break;
        }
        exit (0);
      }
      if (strcmp (argv[iarg], "--test") == 0) {
        test = 1;
      } else if (strcmp (argv[iarg], "--nodoptimize") == 0) {
        dodoptimize = 0;
      } else if (strcmp (argv[iarg], "--xfigfile") == 0) {
        ++iarg;
        xfigproblem = argv[iarg];
      } else if (strcmp (argv[iarg], "--skip") == 0) {
        skiptime = atof(argv[++iarg]);
      } else if (strcmp (argv[iarg], "--skiprtime") == 0) {
        skiprtime = atof(argv[++iarg]);
      } else if (strcmp (argv[iarg], "--ge") == 0) {
        ++iarg;
        gr_lid = 0;
        if (strcmp (argv[iarg], "null") == 0) gr_lid = GO_NULL;
#ifdef HAVE_GTK
        if (strcmp (argv[iarg], "gtk") == 0) gr_lid = GO_GTK;
#endif
#ifdef HAVE_GLUT
        if (strcmp (argv[iarg], "glut") == 0) gr_lid = GO_GLUT;
#endif
        if (strcmp (argv[iarg], "xfig") == 0) gr_lid = GO_XFIG;
        if (gr_lid) gr_id = gr_lid;
          else printf ("Invalid ge value: %s\n", argv[iarg]);
      } else if (strcmp (argv[iarg], "--title") == 0) {
        title = argv[++iarg];
        if (xfigproblem == 0) xfigproblem = argv[iarg];
      } else if (strcmp (argv[iarg], "--k1") == 0) {
        k1_coeff = atof (argv[++iarg]);
      } else if (strcmp (argv[iarg], "--k2") == 0) {
        k2_coeff = atof (argv[++iarg]);
      } else if (strcmp (argv[iarg], "--k3") == 0) {
        k3_coeff = atof (argv[++iarg]);
      } else if (strcmp (argv[iarg], "--k4") == 0) {
        k4_coeff = atof (argv[++iarg]);
      } else if (strcmp (argv[iarg], "--k5") == 0) {
        k5_coeff = atof (argv[++iarg]);
      } else if (strcmp (argv[iarg], "--help") == 0) {
        printf ("usage: %s [--help] [--nodoptimize] [--pause]\n", argv[0]);
        printf ("      [--steps <n>][--title <title>]\n");
        printf ("      [--ge <graphic engine>][--skip <simtime>]\n");
        printf ("      [--skiprtime <seconds>]\n");
        printf ("      [--xfigspecial][--onlyvisible]\n");
        printf ("\nvalid graphic engines (option --ge): ");
#ifdef HAVE_GTK
        printf ("gtk ");
#endif
#ifdef HAVE_GLUT
        printf ("glut ");
#endif
        printf ("xfig\n");
        exit (0);
      } else {
        printf ("Invalid option: %s\n", argv[iarg]);
        exit (1);
      }
      continue;
    }
    if (filein == 0)
    {
      filein = fopen (argv[iarg], "r");
      if (filein == 0)
      {
        perror ("Error opening morse file");
        exit (10);
      }
    } else {
      printf ("Invalid argument: %s\n", argv[iarg]);
      exit (1);
    }
  }
  return (filein);
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
check_contour (struct polyline *contour)
{
  struct line *line, *first, *last, *l;
  struct earc *arc;
  int count, found;

  for (line = contour->line; line; line = line->next)
  {
    arc = line->earc;
    if (arc->loop)
    {
      /* this is an S1 */
      assert (arc->first == 0 && arc->last == 0);
      first = last = arc->loop;
      l = first;
      count = 0;
      found = 0;
      do {
        count++;
        assert (l->earc == arc);
        assert ((l->b->type & V_CROSS) == 0);
        if (l == line) found = 1;
      } while (l = l->b->line[1], l != first);
      assert (count == arc->numsegments);
      assert (found);
    } else {
      first = arc->first;
      last = arc->last;
      assert (first && last);
      assert (arc->loop == 0);
      assert (first->a->type & V_CROSS);
      assert (last->b->type & V_CROSS);
      found = 0;
      if (last == line) found = 1;
      count = 1;
      for (l = first; l != last; l = l->b->line[1])
      {
        count++;
        assert (l->earc == arc);
        assert ((l->b->type & V_CROSS) == 0);
        if (l == line) found = 1;
      }
      assert (count == arc->numsegments);
      assert (found);
    }
  }
}

void
setarcinfo (struct polyline *contour)
{
  struct vertex *p;
  struct line *l;
  int i, iss1;

  /* first deal with the non-s1 arcs */

  for (p = contour->vertex; p; p = p->next)
  {
    if ((p->type & V_CROSS) == 0) continue;
    for (i = 0; i < 4; i++)
    {
      l = p->line[i];
      if (l->a != p) continue;   /* we want a 'leaving' arc */
if (l->earc->loop)
{
  printf ("Questo e' un loop(bis)\n");
  assert (0);
}
      iss1 = setarcinfo1 (l);
      assert (iss1 == 0);
    }
  }

  /* now the s1's */

  for (p = contour->vertex; p; p = p->next)
  {
    if (p->type & V_CROSS) continue;
    l = p->line[1];
    if (l->earc->first) continue;
    iss1 = setarcinfo1 (l);
    assert (iss1);
  }

  /* now check that all arcs are OK */

  for (l = contour->line; l; l = l->next)
  {
    assert (l->earc->loop || (l->earc->first && l->earc->last));
  }
}

int
setarcinfo1 (struct line *l)
{
  struct earc *firstarc;
  struct line *wl;
  struct vertex *p;
  int count;

  firstarc = l->earc;
  if (l->a->type & V_CROSS)
    firstarc->first = l;
  else
    firstarc->loop = l;

  p = 0;

  count = 0;
  wl = l;
  do
  {
    count++;
    if (wl->earc != firstarc)
    {
      wl->earc = mergearcinfo (firstarc, wl->earc);
      if (wl->earc) wl->earc->refcount++;
    }
    if ((p = wl->b)->type & V_CROSS) break;
  } while (wl = p->line[1], wl != l);

  if (p && p->type & V_CROSS) firstarc->last = wl;

  firstarc->numsegments = count;
  if (firstarc->loop) return (1);
  return (0);
}

struct earc *
mergearcinfo (struct earc *arc1, struct earc *arc2)
{
  int conflict = 0;

  //fprintf (stderr, "duplicate arc info\n");
  //if (arc1->orientation != arc2->orientation) conflict++;
  if (arc1->cusps != arc2->cusps) conflict++;

  /* should check also all other information */
  if (conflict > 0)
  {
    fprintf (stderr, "conflicting information for same arc!\n");
    exit (111);
  }
  arc1->numsegments += arc2->numsegments;
  arc2->refcount--;
  if (arc2->refcount <= 0) free (arc2);
  return (arc1);
}

int
enforce_crossover_info (struct polyline *contour)
{
  struct vertex *v;
  int incrnwse, incrnesw, goon;
  int *d0, *d1, *d2, *d3;

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (v = contour->vertex; v; v = v->next)
    {
      if (v->type & V_CROSS)
      {
        if ((v->type & V_CROSS_NWSE) ||
            (v->type & V_CROSS_NESW))
        {
          incrnwse = incrnesw = 0;
          if (v->type & V_CROSS_NWSE) incrnesw = 2;
          if (v->type & V_CROSS_NESW) incrnwse = 2;
          if (v->line[0]->b == v) incrnesw *= -1;
          if (v->line[1]->a == v) incrnwse *= -1;
          d0 = v->line[0]->earc->d;
          d3 = v->line[3]->earc->d;
          if (d0 && d0[0] > MINF && (d3 == 0 || d3[0] <= MINF))
          {
            goon = 1;
            if (d3) fprintf (stderr, "caso non gestito ancora\n");
            d3 = v->line[3]->earc->d = (int *) malloc (sizeof (int));
            if (v->line[0]->a == v)  /* serve il primo valore */
            {
              d3[0] = d0[0] + incrnwse;
            } else {
              d3[0] = d0[v->line[0]->earc->cusps] + incrnwse;
            }
            //fprintf (stderr, "eredito nw->se (%d)\n", incrnwse);
          }
          d1 = v->line[1]->earc->d;
          d2 = v->line[2]->earc->d;
          if (d1 && d1[0] > MINF && (d2 == 0 || d2[0] <= MINF))
          {
            goon = 1;
            if (d2) fprintf (stderr, "caso non gestito ancora\n");
            d2 = v->line[2]->earc->d = (int *) malloc (sizeof (int));
            if (v->line[1]->a == v)  /* serve il primo valore */
            {
              d2[0] = d1[0] + incrnesw;
            } else {
              d2[0] = d1[v->line[1]->earc->cusps] + incrnesw;
            }
            //fprintf (stderr, "eredito ne->sw (%d)\n", incrnesw);
          }
        }
      }
    }
  }
  return (1);
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
    if ((p->type & V_CROSS) == 0) continue;
    for (i = 0; i < 4; i++)
    {
      l = p->line[i];
      if (l->a != p) continue;   /* we want a 'leaving' arc */
      assert (l->earc->cuspsinserted == 0);
      l->earc->cuspsinserted = 1;
      if (l->earc->cusps <= 0) continue;
      iss1 = insert_cusps_on_arc (l);
      assert (iss1 == 0);
    }
  }

  /* now the s1's */

  for (p = contour->vertex; p; p = p->next)
  {
    if (p->type & V_CROSS) continue;
    l = p->line[1];
    if (l->earc->cuspsinserted) continue;
    l->earc->cuspsinserted = 1;
    if (l->earc->cusps <= 0) continue;
    iss1 = insert_cusps_on_arc (l);
    assert (iss1);
  }

  /* now check that no cusps remain to be placed
   * and also place depth value on arcs with no cusps
   */

  for (l = contour->line; l; l = l->next)
  {
    if (l->earc->cusps == 0)
    {
      if (l->earc->d) l->d = l->earc->d[0];
        else l->d = 0;    /* can happen if we don't have huffman labelling */
    }
    assert (l->earc->cuspsinserted);
  }
}

int
insert_cusps_on_arc (struct line *l)
{
  struct line *wl;
  struct vertex *p;
  struct earc *arc;
  int cusps, iss1, count, ccount, nnodes, nums;

  iss1 = 0;
  if ((l->a->type & V_CROSS) == 0) iss1 = 1;
  arc = l->earc;
  cusps = arc->cusps;
  // nnodes = arc->numsegments/(cusps + 1) + 1;
  nums = arc->numsegments;
  nnodes = (nums + 1)/(cusps + 1);

  count = ccount = 0;
  wl = l;
  do {
    count++;
    if (arc->d) wl->d = arc->d[ccount];
    if ((p = wl->b)->type & V_CROSS) break;
    if (count >= nnodes)
    {
      ccount++;
      nums -= count;
      count = 0;
      // nnodes = arc->numsegments*(ccount + 1)/(cusps + 1) + 1;
      nnodes = (nums + 1)/(cusps + 1 - ccount);
      assert (p->type == V_REGULAR);
      p->type = V_CUSP;
      //if (ccount >= cusps) break;
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
static double simtime;

void
settimer (double incrtime, double incsimtime)
{
  int seconds, useconds;

  simtime = -1.0;
  if (incsimtime > 0.0) simtime = contour->time + incsimtime;
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

  if (immediate_exit)
  {
    immediate_exit = 0;
    return (0);
  }
  gettimeofday (&now, 0);
  if (simtime > 0.0 && contour->time > simtime) return (0);
  return (timercmp (&now, &timer, < ));
}

static int allownodered = ALLOW_NODE_REDEFINE;

void
check_timers (struct polyline *contour)
{
  struct timerevent *timer, *prevtimer, *nexttimer;
  static int rncount = 0;
  static int repcount = 0;
  int tc = 0;

  for (timer = timerlist->next; timer; timer = nexttimer)
  {
    tc++;
    nexttimer = timer->next;
    if (contour->time >= timer->time)
    {
      switch (timer->event)
      {
        case EVENT_REDISTRIBUTENODES:
        activate_timer (EVENT_REDISTRIBUTENODES, contour->time + taurn);
        rncount++;
        //printf ("would redistribute nodes: %d...(%lf, %lf)\n", rncount, contour->time, taurn);
        if (allownodered) redistributenodes (contour);
        break;

        case EVENT_REPULSIVEENERGY:
        activate_timer (EVENT_REPULSIVEENERGY, contour->time + taurep);
        repcount++;
        //printf ("would compute repulsive energy: %d...(%lf, %lf)\n", repcount, contour->time, taurep);
        compute_repulsive_energy (contour);
        compute_repulsive_gradient (contour);
#ifdef CHECK_GRADIENT
        check_gradient (contour);
#endif
        break;

        case EVENT_KICKOUT:
        kick_out(contour);
        break;
      }
      for (prevtimer = timerlist; prevtimer; prevtimer = prevtimer->next)
      {
        if (prevtimer->next == timer) prevtimer->next = timer->next;
      }
      free (timer);
    }
  }
}

double
evolve (struct polyline *contour, double incrtime, double incsimtime)
{
  double newenergy, decrease, predicted_decrease;
  struct vertex *v;
  int gradient_is_ok, gradient_is_reliable;
  int timesteps = 0;
  int timebsteps = 0;
  int allowbackstep = ALLOW_BACKSTEP;
  int allowstepcontrol = ALLOW_STEPCONTROL;

  settimer (incrtime, incsimtime);
  gradient_is_ok = 0;
  while (checktimer())
  {
    check_timers (contour);
    if (curenergy == 0.0) curenergy = compute_energy (contour);
    if (! gradient_is_ok) compute_gradient (contour);
    timesteps++;
    for (v = contour->vertex; v; v = v->next)
    {
      v->x -= tau*contour->gradx[v->tag];
      v->y -= tau*contour->grady[v->tag];
    }
    gradient_is_reliable = 1;
    if (k4_coeff && contour->time != contour->rgradtime) 
        gradient_is_reliable = 0;
    contour->time += tau;
    gradient_is_ok = 0;
    if (k4_coeff && gradient_is_reliable) compute_repulsive_energy (contour);
    newenergy = compute_energy (contour);
    decrease = curenergy - newenergy;
    if (allowbackstep && decrease < 0.0)
    {
//      fprintf (stderr, "WARNING: energy is increasing!\n");
      /* torno indietro di uno step e diminuisco tau */
      for (v = contour->vertex; v; v = v->next)
      {
        v->x += tau*contour->gradx[v->tag];
        v->y += tau*contour->grady[v->tag];
      }
      contour->time -= tau;
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

    if (allowstepcontrol && gradient_is_reliable)
    {
 /* automatic step-size control does not work yet! */
      //if (fabs(predicted_decrease - decrease) < 0.05*predicted_decrease)
      if (decrease > 0.95*predicted_decrease)
      {
        tau *= 1.1;
        //printf ("good prediction, increasing tau %lg\n", tau);
      }

      //if (fabs(predicted_decrease - decrease) > 0.2*predicted_decrease)
      if (decrease < 0.4*predicted_decrease)
      {
        if (tau > 1e-12) 
        {
          tau /= 2;
          //printf ("bad prediction, decreasing tau, %lg\n", tau);
        } else {
          toggle_motion (1);
          immediate_exit = 1;
          printf ("stopping motion due to a very small time step\n");
        }
      }
    }
  }
//  printf ("timesteps: %d, backsteps: %d, time = %lf, energy = %lf\n", 
//          timesteps, timebsteps, time, curenergy);
  return (contour->time);
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

  if ((p->type & V_CROSS) == 0) return (p->line[1]);
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

  if ((p->type & V_CROSS) == 0) return (p->line[0]);
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
    if (v->type & V_CROSS) continue;   /* non tocco l'ordine */
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
  contour->baseh = contour->h;
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
  struct vertex *a, *b, *v, *next, *prev;
  double diffx, diffy, len, alen, h, dt;
  int numsub, i, vertexnum, forward, backidx;
  int addednodes = 0, removednodes = 0;
  int ns;

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
      prev = line->a;
      forward = 1;
      if (v->type != V_REGULAR)
      {
        v = line->a;
        prev = line->b;
        forward = 0;
      }
      if (v->type != V_REGULAR) continue;
      aline = v->line[forward];
      assert (aline != line);
      next = aline->a;
      if (forward) next = aline->b;
      assert (next != v);
      backidx = 1 - forward;
      if (next->type & V_CROSS)
      {
        for (backidx = 0; backidx < 4; backidx++)
          if (next->line[backidx] == aline) break;
      }
      ns = aline->earc->numsegments;
      if (ns < 4)
      {
        kick_out(contour);
        //fprintf (stderr, "arc with less than 4 segments\n");
      }
      diffx = v->x - next->x;
      diffy = v->y - next->y;
      alen = sqrt (diffx*diffx + diffy*diffy);
      if (alen + len < 0.9*h && ns >= 4)   /* can derefine */
      { /* aline will be removed! */
        removednodes++;
        line->earc->numsegments--;
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
        line->earc->numsegments++;
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
  //if (addednodes) printf ("added %d nodes\n", addednodes);
  //if (removednodes) printf ("removed %d nodes\n", removednodes);
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
  nline->d = line->d;
  nline->earc = line->earc;
  nline->earc->refcount++;
  nline->rarc = 0;
  if (line->rarc)
  {
    nline->rarc = line->rarc;
    nline->rarc->numsegments++;
    if (line->rarc->last == line) line->rarc->last = nline;
  }
  nline->next = line->next;
  nline->a = p;
  nline->b = b;
  nl = 2;
  if (b->type & V_CROSS) nl = 4;
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
init_rarc (struct polyline *contour)
{
  struct line *line, *l, *last;
  struct vertex *b;
  struct earc *earc;
  struct rarc *rarc;
  int icus, count, canbreak;

  for (line = contour->line; line; line = line->next)
  {
    if (line->rarc) continue;
    earc = line->earc;
    if (earc->first != line && earc->loop != line) continue;
    /* let's walk along the arc */
    last = earc->last;
    if (last == 0)
    {
      last = earc->loop;
      last = last->a->line[0];
    }
    icus = 0;
    rarc = (struct rarc *) malloc (sizeof (struct rarc));
    if (earc->d) rarc->d = earc->d[icus];
      else rarc->d = 0;
    rarc->first = line;
    rarc->loop = rarc->last = 0;
    rarc->numsegments = 0;
    l = line;
    count = 0;
    while (1)
    {
      b = l->b;
      count++;
      l->rarc = rarc;
      rarc->numsegments++;
      canbreak = 1;
      if (b->type == V_CUSP && l != last)
      {
        if (earc->loop && icus == 0)
        {  /* this is the first cusp of a loop
            * we want to start here... 
            */
          last = l;
          canbreak = 0;
          count = 0;
        } else {
          rarc->last = l;
          rarc = (struct rarc *) malloc (sizeof (struct rarc));
        }
        icus++;
        rarc->d = earc->d[icus];
        rarc->first = b->line[1];
        rarc->loop = rarc->last = 0;
        rarc->numsegments = 0;
      }
      if (l == last && canbreak) break;
      l = b->line[1];
    }
    rarc->last = l;
    if (earc->loop && earc->cusps == 0)
    {
      assert (l->b->line[1] == rarc->first);
      rarc->loop = rarc->first;
      rarc->first = rarc->last = 0;
    }
    assert (count == earc->numsegments);
    assert (icus == earc->cusps);
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
  line->earc = 0;
  line->rarc = 0;
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
  if (type & V_CROSS) numarcs = 4;
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
  struct rarc *rarc;

  if (contour->line == line)
  {
    contour->line = line->next;
  } else {
    for (l = contour->line; l; l = l->next)
    {
      if (l->next == line)
      {
        l->next = line->next;
        break;
      }
    }
  }
  line->earc->refcount--;
  if (line->earc->refcount <= 0) free (line->earc);
  if ((rarc = line->rarc))
  {
    rarc->numsegments--;
    if (rarc->loop == line)
      rarc->loop = line->b->line[1];  /* there cannot be crosses for a loop */
    if (rarc->first == line)
      rarc->first = line->b->line[1];  /* there cannot be a cross here */
    if (rarc->last == line)
      rarc->last = line->a->line[0];  /* there cannot be a cross here */

    assert (rarc->first != line || rarc->numsegments == 0);
    assert (rarc->last != line || rarc->numsegments == 0);
    assert (rarc->loop != line || rarc->numsegments == 0);
    if (rarc->numsegments == 0) free (rarc);
  }
  free (line);
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
buildpolyline (FILE *filein)
{
  double dx, dy, dx1, dy1, dx12, dy12, x, y, xstart, ystart;
  struct morseevent morseevent;
  struct vertex *danglingnodes[BUFSIZE];
  struct vertex *prevdanglingnodes[BUFSIZE];
  struct vertex *v1, *v2, *v4, *v5, *v, *dangv;
  struct line *line, *dangline;
  struct polyline *contour;
  int goon, i, ii, k, prevdangnodes, dangind, prevdangind;
  int oriented, dangori;
  int counttrans = 0, after, type;

  contour = (struct polyline *) malloc (sizeof (struct polyline));
  contour->vertex = 0;
  contour->line = 0;
  contour->gradx = contour->grady = 0;
  contour->specnodes = 0;
  contour->specnodesnum = 0;
  contour->renergy = 0;
  contour->rgradx = contour->rgrady = 0;
  contour->rgraddim = 0;
  contour->rentime = contour->rgradtime = -1.0;

  contour->h = dx = dy = 1.0;    /* aggiusto alla fine */
  i = 0;            /* conta gli eventi su ciascuna riga */
  prevdangnodes = prevdangind = dangind = 0;
  y = x = ystart = xstart = 0.0;
  dx1 = dx12 = dy1 = dy12 = 0.0;
  goon = 1;
  after = 0;
  while (goon)
  {
    i++;
    x += dx;
    y -= dy;
    getmorseevent (filein, &morseevent, 1);
    if (morseevent.type != ME_TRAN)
    {
      switch (morseevent.type)
      {
        case ME_TOP:
        case ME_BOT:
        case ME_CROSS:
        case ME_CROSS_NWSE:
        case ME_CROSS_NESW:
          assert (after == 0);
          after = 1;
          break;
      }
      switch (morseevent.type)
      {
        case ME_TOP:
          dx1 = dx12 = 0.0;
          dy1 = dy;
          dy12 = 2*dy;
          break;
        case ME_BOT:
          dx1 = dx;
          dx12 = 2*dx;
          dy1 = dy12 = 0.0;
          break;
        case ME_CROSS:
        case ME_CROSS_NWSE:
        case ME_CROSS_NESW:
          dx1 = dx12 = dx;
          dy1 = 0.0;
          dy12 = dy;
          break;

        case ME_NEWROW:
        case ME_LASTROW:
          if (after)
          {
            after = 0;
            dx1 = dx - dx1;
            dy1 = dy - dy1;
            dx12 = 2*dx - dx12;
            dy12 = 2*dy - dy12;
          } else {
            counttrans = 0;
            morseevent.type = ME_RESET;
          }
          break;
      }
      // printf ("counted %d consecutive trans\n", counttrans);
      prevdangind -= counttrans;
      x -= counttrans*dx;
      y += counttrans*dy;
      for (ii = 0; ii < counttrans; ii++)
      {
        v1 = newvertex (contour, x + dx1, y + dy1, V_REGULAR);
        v2 = newvertex (contour, x + dx12, y + dy12, V_REGULAR);
        line = newline (contour, v1, v2);
        v1->line[0] = line;
        v2->line[0] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v1);
        v1->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        danglingnodes[dangind++] = v2;
        x += dx;
        y -= dy;
      }
      counttrans = 0;
    }
    switch (morseevent.type)
    {
      case ME_LASTROW:
        goon = 0;
        i--;
        if (prevdangnodes != prevdangind)
        {
          fprintf (stderr, "dangling nodes: prev = %d, current = %d\n",
                  prevdangnodes, dangind);
        }

        if (dangind != 0) fprintf (stderr, "Dangling nodes remain...\n");
        break;

      case ME_NEWROW:
        xstart += 2*dx - dx12;
        ystart += 2*dy - dy12;
        x = xstart;
        y = ystart;
        i--;
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
        x -= dx;
        y += dy;
        v = newvertex (contour, x + dx, y, V_REGULAR);
        v1 = newvertex (contour, x + dx, y + dy, V_REGULAR);
        v2 = newvertex (contour, x + 2*dx, y, V_REGULAR);
        line = newline (contour, v, v1);
        line->orientation = morseevent.ori;
        line->earc = morseevent.arc;
        if (line->earc) line->earc->refcount++;
        v1->line[0] = line;
        v->line[0] = line;
        line = newline (contour, v, v2);
        v2->line[0] = line;
        v->line[1] = line;
        danglingnodes[dangind++] = v1;
        danglingnodes[dangind++] = v2;
        break;
      case ME_BOT:
        x += dx;
        y -= dy;
        v = newvertex (contour, x, y + dy, V_REGULAR);
        line = newline (contour, prevdanglingnodes[prevdangind], v);
        v->line[0] = prevdanglingnodes[prevdangind++]->line[1] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v);
        line->orientation = morseevent.ori;
        line->earc = morseevent.arc;
        if (line->earc) line->earc->refcount++;
        v->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        break;

      case ME_CROSS_NWSE:
      case ME_CROSS_NESW:
        //fprintf (stderr, "warning: crossover modifier not implemented\n");
      case ME_CROSS:
        type = V_CROSS;
        if (morseevent.type == ME_CROSS_NWSE) type += V_CROSS_NWSE;
        if (morseevent.type == ME_CROSS_NESW) type += V_CROSS_NESW;
        x += dx;
        y -= dy;
        v1 = newvertex (contour, x, y + dy, type);
        v4 = newvertex (contour, x, y + 2*dy, V_REGULAR);
        v5 = newvertex (contour, x + dx, y + dy, V_REGULAR);
        line = newline (contour, v1, v4);
        line->orientation = morseevent.ori;
        line->earc = morseevent.arc;
        if (line->earc) line->earc->refcount++;
        v1->line[2] = line;
        v4->line[0] = line;
        line = newline (contour, v1, v5);
        line->orientation = morseevent.ori2;
        line->earc = morseevent.arc2;
        if (line->earc) line->earc->refcount++;
        v1->line[3] = line;
        v5->line[0] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v1);
        v1->line[0] = prevdanglingnodes[prevdangind++]->line[1] = line;
        line = newline (contour, prevdanglingnodes[prevdangind], v1);
        v1->line[1] = prevdanglingnodes[prevdangind++]->line[1] = line;
        danglingnodes[dangind++] = v4;
        danglingnodes[dangind++] = v5;
        break;

      case ME_TRAN:
        counttrans++;
        /* we want here only to deal with
         * the (possible) arc information, attaching it
         * to dangling and already defined arcs
         * all other work has been moved in the above
         * section cumulatively
         */
        dangv = prevdanglingnodes[prevdangind++];
        dangline = dangv->line[0];
        dangori = dangline->orientation;
        if (dangline->b != dangv) dangori *= -1;
        if (morseevent.arc) morseevent.arc->refcount++;
        if (morseevent.ori*dangori != 0
            && morseevent.ori != dangori)
          fprintf (stderr, "Conflicting orientation!\n");
        if (dangline->earc == 0) dangline->earc = morseevent.arc;
          else if (morseevent.arc) 
            dangline->earc = mergearcinfo (dangline->earc, morseevent.arc);
        if (dangline->orientation == 0) 
          dangline->orientation = morseevent.ori;
        if (dangline->b != dangv) dangline->orientation *= -1;

        break;

      case ME_RESET:
        x = xstart;
        y = ystart;
        prevdangind = dangind = 0;
        break;

      default:
        fprintf (stderr, "Invalid morse event\n");
    }
  }

  for (v = contour->vertex; v; v = v->next)
  {
    assert (v->line[0]);
    assert (v->line[1]);
    if (v->type & V_CROSS) assert (v->line[2] && v->line[3]);
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

  /* inherit arc info */

  goon = 1;
  while (goon)
  {
    while (goon)
    {
      goon = 0;
      for (line = contour->line; line; line = line->next)
      {
        if (line->earc != 0 && inherit_earc (line)) goon = 1;
      }
    }
    for (line = contour->line; line; line = line->next)
    {
      if (line->earc == 0)
      {
        line->earc = (struct earc *) malloc (sizeof (struct earc));
        line->earc->cusps = line->earc->cuspsinserted = 0;
        line->earc->first = line->earc->last = line->earc->loop = 0;
        line->earc->refcount = 1;
        line->earc->d = 0;
        goon = 1;
        break;
        // fprintf (stderr, "arc without earc info\n");
      }
    }
  }

  if (oriented == 0) fprintf (stderr, "Cannot orient arc\n");

  return (contour);
}

void
renormalize (struct polyline *contour)
{
  struct vertex *v;
  double xmax, ymax, xmin, ymin, dx, dy;

  /* compute containing rectangle */
  xmax = xmin = contour->vertex->x;
  ymax = ymin = contour->vertex->y;
  for (v = contour->vertex; v; v = v->next)
  {
    if (xmax < v->x) xmax = v->x;
    if (xmin > v->x) xmin = v->x;
    if (ymax < v->y) ymax = v->y;
    if (ymin > v->y) ymin = v->y;
  }
  //printf ("(%lf, %lf), (%lf, %lf)\n", xmin, ymin, xmax, ymax);
  /* rinormalizzazione ascisse e ordinate per stare in [-1,1]x[-1,1] */
  dx = 2.0/(xmax - xmin);
  dy = 2.0/(ymax - ymin);
  contour->h = dx;
  if (dy < dx) contour->h = dy;

  for (v = contour->vertex; v; v = v->next)
  {
    v->x = dx*v->x - 1.0;
    v->y = 1.0 - dy*v->y;
  }
}

/* inherit orientation to directly adjacent arcs */
/* also inherit across nodes! */

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
  if ((v->type & V_CROSS) == 0)
  {
    wline = v->line[0];
    if (wline == line) wline = v->line[1];
  } else wline = prevp (line, v);
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
  }

  v = line->b;  /* going forward */

  if ((v->type & V_CROSS) == 0)
  {
    wline = v->line[0];
    if (wline == line) wline = v->line[1];
  } else wline = nextp (line, v);
  assert (wline && wline != line);
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
  }

  return (goon);
}

/* inherit earc info to directly adjacent arcs */
/* we can assume that arcs are all oriented (but not the nodes) */

int
inherit_earc (struct line *line)
{
  struct line *wline;
  struct vertex *v;
  int goon;

  goon = 0;
  v = line->a;  /* going back */
  if ((v->type & V_CROSS) == 0)
  {
    wline = v->line[0];
    if (wline == line) wline = v->line[1];
    assert (wline && wline != line);
    if (wline->earc == 0)
    {
      goon = 1;
      wline->earc = line->earc;
      wline->earc->refcount++;
    }
  }

  v = line->b;  /* going forward */

  if ((v->type & V_CROSS) == 0)
  {
    wline = v->line[0];
    if (wline == line) wline = v->line[1];
    assert (wline && wline != line);
    if (wline->earc == 0)
    {
      goon = 1;
      wline->earc = line->earc;
      wline->earc->refcount++;
    }
  }

  return (goon);
}

char *
grinit (int *argcpt, char *argv[])
{
  static char idxfig[] = "xfig";
  static char idnull[] = "null";
  static char idunknown[] = "unknown";

  switch (gr_id)
  {
    case GO_XFIG:
    return (idxfig);

    case GO_NULL:
    return (idnull);

#ifdef HAVE_GTK
    case GO_GTK:
    return (gtk_grinit (argcpt, argv));
#endif
#ifdef HAVE_GLUT
    case GO_GLUT:
    return (glut_grinit (argcpt, argv));
#endif
    default:
    return (idunknown);
  }
}

extern int steps;

int
grmain (void)
{
  int i;
  static double incrtime = 0.25;

  switch (gr_id)
  {
    case GO_XFIG:
    xfig_export0 (contour, xfigproblem, &grflags);
    if (steps >= 10000) steps = 0;
    if (steps)
    {
      for (i = 0; i < steps; i++)
      {
        evolve (contour, incrtime, 0.0);
        xfig_export0 (contour, xfigproblem, &grflags);
      }
    }
    return (0);

    case GO_NULL:
    return (0);

#ifdef HAVE_GTK
    case GO_GTK:
    return (gtk_grmain ());
#endif
#ifdef HAVE_GLUT
    case GO_GLUT:
    return (glut_grmain ());
#endif
    default:
    return (0);
  }
}

void
toggle_motion (int motion)
{
  switch (gr_id)
  {
    case GO_XFIG:
    return;

    case GO_NULL:
    return;

#ifdef HAVE_GTK
    case GO_GTK:
    return (gtk_toggle_motion (motion));
#endif
#ifdef HAVE_GLUT
    case GO_GLUT:
    return (glut_toggle_motion (motion));
#endif
    default:
    return;
  }
}
