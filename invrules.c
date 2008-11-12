#include <assert.h>
#include <strings.h>
#include "contour.h"

/*
 * definizione regole di trasformazione inverse per superfici isotope
 */

extern int debug;

/* local prototypes */
int c_mergearcs (struct sketch *s, struct region *r,
	struct arc *a1, struct arc *a2, int a1l, int a2l);
int c_createwrinkle (struct sketch *s, struct region *r, int stratum);
int common_work_mergearcs (struct sketch *s, 
			   struct border *bp1, struct border *bp2, 
                           int a1pos, int a2pos, int rule);

/*
 * creazione di un wrinkle (inv C1)
 */

static int countwrrules = 0;
static int applywr = 0;
static int applywrc = 0;

int
list_strata (struct sketch *s)
{
  int res;
  struct region *r = 0;
  extern int quiet;
  extern struct tagged_data user_data;

  if (user_data.mrnum > 0) {
    r = findregion (s, user_data.region[0]);
    if (r == 0) {fprintf (stderr, "Cannot find region %d\n",
      user_data.region[0]); return (0);}
  }

  if (r && quiet == 0) printf ("Restricted to region %d\n", r->tag);
  countwrrules = 0;
  applywr = applywrc = 0;
  res = c_createwrinkle (s, r, -1);
  if (quiet) printf ("\n");
  return (res);
}

int
c_createwrinkle (struct sketch *s, struct region *r, int stratum)
{
  extern int quiet;
  int res;

  if (r == 0) {
    for (r = s->regions; r; r = r->next) {
      if (quiet == 0 && applywr == 0) printf ("Region %d:\n", r->tag);
      res = c_createwrinkle (s, r, -1);
      if (res) return (res);
    }
    return (0);
  }

  if (stratum < 0) {
    for (stratum = 0; stratum < r->f; stratum++) {
      res = c_createwrinkle (s, r, stratum);
      if (res) return (res);
    }
    return (0);
  }

  assert (stratum >= 0);
  res = apply_createwrinkle (s, r, stratum, 1);

  if (res) {
    countwrrules++;
    if (applywr == 0) {
      if (quiet == 0)
        printf ("-r %d --stratum %d (", r->tag, stratum);
      printf ("INVC1:%d", countwrrules);
      if (quiet == 0) printf (")\n");
        else printf (" ");
    }
    if (applywr && countwrrules == applywrc) {
      res = apply_createwrinkle (s, r, stratum, 0);
      return (res);
    } else {
      return (0);
    }
  }

  return (0);
}

int
apply_createwrinkle (struct sketch *s, struct region *r, 
        int stratum, int test)
{
  int res;
  extern int verbose;

  assert (s && r && stratum >= 0);
  if (stratum >= r->f) return (0);
  if (test) return (1);
printf ("NOT IMPLEMENTED\n");
return (0);
}

/*
 * Questa funzione serve ad applicare una mossa inversa
 * per il merge di due archi (inversa di N1-N4 oppure C2)
 */

#define INV_N1 1
#define INV_N2 2
#define INV_N3 3
#define INV_N3bis 4
#define INV_N4 5
#define INV_C2 6

#define NUM_INVNC 7

static char *invmergerules[] = {
  "",
  "INVN1",
  "INVN2",
  "INVN3",
  "INVN3bis",
  "INVN4",
  "INVC2",
  0};

static int countmarules[NUM_INVNC];
static int applyma = 0;
static int applymac = 0;

int
lookup_mergearcs (char *rule)
{
  int i;

  for (i = 1; invmergerules[i]; i++)
  {
    if (strcasecmp (rule, invmergerules[i]) == 0) return (i);
  }
  return (0);
}

int
rule_mergearcs (struct sketch *s, int rule, int rcount)
{
  int i, res;
  struct region *r = 0;
  extern struct tagged_data user_data;

  if (debug) printf ("Chiamato rule_mergearcs con rule: %s rcount %d\n",
    invmergerules[rule], rcount);

  if (user_data.mrnum > 0) {
    r = findregion (s, user_data.region[0]);
    if (r == 0) {fprintf (stderr, "Cannot find region %d\n",
      user_data.region[0]); return (0);}
  }

  for (i = 1; i < NUM_INVNC; i++) countmarules[i] = 0;
  applyma = rule;
  applymac = rcount;
  res = c_mergearcs (s, r, 0, 0, -1, -1);
  return (res);
}

int
list_mergearcs (struct sketch *s)
{
  int i, res;
  struct region *r = 0;
  struct arc *a1 = 0;
  struct arc *a2 = 0;
  int a1l = -1;
  int a2l = -1;
  extern int quiet;
  extern struct tagged_data user_data;

  if (user_data.mrnum > 0) {
    r = findregion (s, user_data.region[0]);
    if (r == 0) {fprintf (stderr, "Cannot find region %d\n",
      user_data.region[0]); return (0);}
  }

  if (r && quiet == 0) printf ("Restricted to region %d\n", r->tag);
  for (i = 1; i < NUM_INVNC; i++) countmarules[i] = 0;
  applyma = applymac = 0;
  res = c_mergearcs (s, r, a1, a2, a1l, a2l);
  if (quiet) printf ("\n");
  return (res);
}

int
c_mergearcs (struct sketch *s, struct region *r,
	struct arc *a1, struct arc *a2, int a1l, int a2l)
{
  struct borderlist *bl;
  struct border *bp;
  int imax, i, res, rule;
  extern int quiet;

  if (r == 0) {
    for (r = s->regions; r; r = r->next) {
      if (quiet == 0 && applyma == 0) printf ("Region %d:\n", r->tag);
      res = c_mergearcs (s, r, 0, 0, -1, -1);
      if (res) return (res);
    }
    return (0);
  }

  if (a1 == 0 || (a1l >= 0 && a2 == 0)) {
    for (bl = r->border; bl; bl = bl->next)
    {
      if (bl->sponda == 0) continue;
      bp = bl->sponda;
      do {
        res = c_mergearcs (s, r, 
		(a1)?a1:bp->info, (a1)?bp->info:0, a1l, -1);
	if (res) return (res);
	bp = bp->next;
      } while (bp != bl->sponda);
    }
    return (0);
  }

  assert (s && r && a1);
  if (a1l < 0) {
    imax = a1->dvalues;
    for (i = 0; i < imax; i++) {
      res = c_mergearcs (s, r, a1, 0, i, -1);
      if (res) return (res);
    }
    return (0);
  }

  assert (a2 && a1l >= 0);
  if (a2l < 0) {
    imax = a2->dvalues;
    for (i = 0; i < imax; i++) {
      res = c_mergearcs (s, r, a1, a2, a1l, i);
      if (res) return (res);
    }
    return (0);
  }

  assert (a2l >= 0);
  rule = apply_mergearcs (s, r, a1, a2, a1l, a2l, 1);
  if (rule) {
    countmarules[rule]++;
    if (applyma == 0) {
      if (quiet == 0)
        printf ("-r %d -a %d:%d -a %d:%d (", 
          r->tag, a1->tag, a1l, a2->tag, a2l);
      printf ("%s:%d", invmergerules[rule], countmarules[rule]);
      if (quiet == 0) printf (")\n");
        else printf (" ");
    }
    if (rule == applyma && countmarules[rule] == applymac) {
      res = apply_mergearcs (s, r, a1, a2, a1l, a2l, 0);
      return (res);
    } else {
      return (0);
    }
  }

  return (0);
}

int
apply_mergearcs (struct sketch *s, struct region *r,
	struct arc *a1, struct arc *a2, int a1l, int a2l, int test)
{
  struct borderlist *bl;
  struct border *bp;
  struct border *bp1 = 0;
  struct border *bp2 = 0;
  int d1 = -1000;
  int d2 = -1000;
  int deltad, res;
  extern int verbose;

  assert (s && r && a1 && a2 && a1l >= 0 && a2l >= 0);

  /* consistency check about the number of dvalues */
  if (a1l >= 0 && a1l < a1->dvalues) d1 = a1->depths[a1l];
  if (a2l >= 0 && a2l < a2->dvalues) d2 = a2->depths[a2l];

  if (d1 < 0) fprintf (stderr, "invalid subarc (%d) for first arc\n", a1l);
  if (d2 < 0) fprintf (stderr, "invalid subarc (%d) for second arc\n", a2l);
  if (d1 > d2) { /* ensure d1 <= d2 */
    if (test == 0) fprintf (stderr, 
	"d value of first arc cannot exceed d value of the second\n");
    return (0);
  }

  /* check if given arcs bound the given region */
  for (bl = r->border; bl; bl = bl->next)
  {
    if (bl->sponda == 0) continue;
    bp = bl->sponda;
    do {
      if (bp->info == a1) bp1 = bp;
      if (bp->info == a2) bp2 = bp;
      bp = bp->next;
    } while (bp != bl->sponda);
  }
  if (bp1 == 0) 
    fprintf (stderr, "Arc %d does not bound given region\n", a1->tag);
  if (bp2 == 0) 
    fprintf (stderr, "Arc %d does not bound given region\n", a2->tag);

  if (! (bp1 && bp2 && d1 >= 0 && d2 >= d1)) return (0);

  if (debug) {
    if (bp1->border == bp2->border) {
      fprintf (stderr, 
	"Arcs belong to the same c.c. of the boundary of region.\n");
      fprintf (stderr, "The region will be splitted in two\n");
    } else {
      fprintf (stderr,
	"Arcs belong to different c.c. of the boundary of region.\n");
      fprintf (stderr, "The number of c.c. will decrease by one\n");
    }
  }

  assert (d1 <= d2);
  deltad = d2 - d1;
  if (bp1->orientation < 0 && bp2->orientation < 0)
  {
    if (verbose) {
      fprintf (stderr, "Negative orientations: can apply inv N1\n");
      fprintf (stderr, "with first arc above second arc\n");
    }
    if (test) return (INV_N1);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N1);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation > 0 && deltad >= 2)
  {
    if (verbose) {
      fprintf (stderr, "Positive orientations, delta d >= 2,\n");
      fprintf (stderr, "can apply inv N4\n");
    }
    if (test) return (INV_N4);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N4);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation > 0 && deltad == 1)
  {
    if (verbose) {
      fprintf (stderr, "Positive orientations, delta d = 1,\n");
      fprintf (stderr, "can apply inv C2\n");
    }
    if (test) return (INV_C2);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_C2);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation > 0 && deltad == 0)
  {
    if (test == 0) {
      fprintf (stderr, "Positive orientations with delta d = 0\n");
      fprintf (stderr, "   Cannot merge arcs\n");
    }
    return (0);
  }

  if (bp1->orientation > 0 && bp2->orientation < 0 && deltad == 0)
  {
    if (verbose) {
      fprintf (stderr, "Pos-neg orientations with delta d = 0\n");
      fprintf (stderr, "can apply inv N3\n");
    }
    /* if (test) return (INV_N3bis);
     */
    if (test) return (0); /* do not report a duplicate of INV_N3 */
    res = common_work_mergearcs (s, bp2, bp1, a2l, a1l, INV_N3);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation < 0 && deltad == 1)
  {
    if (test == 0) {
      fprintf (stderr, "Pos-neg orientations with delta d = 1\n");
      fprintf (stderr, "   Cannot merge arcs\n");
    }
    return (0);
  }

  if (bp1->orientation > 0 && bp2->orientation < 0 && deltad >= 2)
  {
    if (verbose) {
      fprintf (stderr, "Pos-neg orientations with delta d >= 2\n");
      fprintf (stderr, "can apply inv N2\n");
    }
    if (test) return (INV_N2);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N2);
    return (res);
  }

  if (bp1->orientation < 0 && bp2->orientation > 0)
  {
    if (verbose) {
      fprintf (stderr, "Neg-pos orientations\n");
      fprintf (stderr, "can apply inv N3\n");
    }
    if (test) return (INV_N3);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N3);
    return (res);
  }

  return (0);
}

int
common_work_mergearcs (struct sketch *s, 
		struct border *bp1, struct border *bp2, 
		int a1l, int a2l, int rule)
{
  int res, changes;

  pinch_arcs (&bp1, a1l, &bp2, a2l, s, (rule == INV_C2)?(-1):0);

  if (debug) printf ("dopo pinch_arcs\n");
  if (debug) printsketch (s);

  assert (bp1 != bp2);

  if (rule == INV_C2) {
    bp1->next->info->depths[0] = bp2->next->info->depths[1];
    bp2->next->info->depths[0] = bp1->next->info->depths[1];
    taglia_nodo (bp1->next, s, 0, 0);
    res = 1;
  } else {
    res = aggiungi_losanga (bp1, bp2, s);
    changes = adjust_isexternalinfo (s);
    if (debug && changes) 
      printf ("%d changes in adjust_isexternalinfo\n", changes);
    if (debug) printf ("dopo aggiungi_losanga\n");
    if (debug) printsketch (s);
  }
  if (res)
  {
    checkconsistency (s);
    postprocesssketch (s);
    canonify (s);
  }
  return (res);
}

/*
 * contrario di rimuovi_losanga:
 *
 * \bp2  /    \   /
 *  \   /      \ /
 *   \ /        X
 *    X   ==>  / \
 *   / \       \ /
 *  /   \       X
 * /  bp1\     / \
 *
 * si conviene che il valore piu' basso di d
 * venga assegnato alla parte indicata da bp1
 */

int
aggiungi_losanga (struct border *bp1, struct border *bp2,
		struct sketch *s)
{
  struct region *newr;
  struct borderlist *newbl;
  struct border *bp1nt, *bp2nt, *newbp1, *newbp2;
  struct border *newbpt1, *newbpt2;
  struct arc *newa1, *newa2;
  int o1, o2, d1, d2;

  o1 = bp1->orientation;
  assert (o1 == 1 || o1 == -1);
  o2 = bp2->orientation;
  assert (bp2->next->orientation == o1);
  assert (bp1->next->orientation == o2);

  /*
   * la creazione della losanga richiede numerose
   * nuove strutture dati: due archi, una regione
   * e 4 nuove sponde.
   */

  newa1 = newarc (s);
  newa2 = newarc (s);
  newa1->depths = (int *) malloc (sizeof (int));
  newa2->depths = (int *) malloc (sizeof (int));
  newa1->depthsdim = newa2->depthsdim = 1;  /* niente cuspidi */
  newa1->cusps = newa2->cusps = 0;
  newa1->dvalues = newa2->dvalues = 1;
  newa1->endpoints = newa2->endpoints = 2;
  if (o1 > 0) d1 = bp1->info->depths[bp1->info->dvalues - 1];
    else d1 = bp1->info->depths[0];
  if (o2 > 0) d2 = bp2->info->depths[bp2->info->dvalues - 1];
    else d2 = bp2->info->depths[0];
  /* calcolo i due valori di d. Arco 1 passa davanti */
  newa1->depths[0] = d1;   /* questo e' giusto perche' passa davanti */
  newa2->depths[0] = d2 - 2*o1;  /* varia di 2 con segno che dipende */

  newr = newregion (s);
  newbl = newborderlist (newr);
  newbp1 = newborder (newbl);
  newbl->sponda = newbp1;
  newbp2 = newborder (newbl);
  newbp2->next = newbp1;  /* o anche = newbp1->next */
  newbp1->next = newbp2;

  bp1nt = gettransborder (bp1->next);  /* sponda opposta al next */
  bp2nt = gettransborder (bp2->next);
  newbpt1 = newborder (bp1nt->border); /* e' la sponda opposta a bpn1 */
  newbpt1->next = bp1nt->next;  /* da inserire dopo bp1nt */
  bp1nt->next = newbpt1;
  newbp1->info = newbpt1->info = newa1;
  newbp1->orientation = -o1;
  newbpt1->orientation = o1;
  newbpt2 = newborder (bp2nt->border); /* e' la sponda opposta a bpn1 */
  newbpt2->next = bp2nt->next;
  bp2nt->next = newbpt2;
  newbp2->info = newbpt2->info = newa2;
  newbp2->orientation = -o2;
  newbpt2->orientation = o2;

  /* infine sistemo i puntatori alle regioni right/left */
  if (o1 > 0) {
    newa1->regionright = newbp1;
    newa1->regionleft = newbpt1;
  } else {
    newa1->regionleft = newbp1;
    newa1->regionright = newbpt1;
  }
  if (o2 > 0) {
    newa2->regionright = newbp2;
    newa2->regionleft = newbpt2;
  } else {
    newa2->regionleft = newbp2;
    newa2->regionright = newbpt2;
  }
  return (1);
}

