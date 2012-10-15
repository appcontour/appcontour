#include <assert.h>
#include "contour.h"

extern int debug;

int
arcmult (struct arc *arc)
{
  return ((arc->regionleft->border->region->f +
           arc->regionright->border->region->f)/2);
}

int
sketchcmp (struct sketch *s1, struct sketch *s2)
{
  struct region *r1, *r2;
  int diff, res;

  /* set di criteri estetici */
  /* primo: numero archi invisibili contati con la d */
  diff = count_hidden_arcs(s1) - count_hidden_arcs(s2);
  if (diff > 0) return (1);   /* s1 ha piu' archi nascosti */
  if (diff < 0) return (-1);  /* s2 ha piu' archi nascosti */

  /* set di criteri tecnici */
  /* primo criterio: numero di regioni */
  for (r1 = s1->regions, r2 = s2->regions;
       r1 && r2; r1 = r1->next, r2 = r2->next);
  if (r1 != 0) return (1);   /* s1 ha piu' regioni */
  if (r2 != 0) return (-1);  /* s2 ha piu' regioni */

  /* secondo criterio: confronto lessicografico tra ciascuna regione */
  /* per ultima confronto la regione esterna! */

  for (r1 = s1->regions->next, r2 = s2->regions->next;
       r1 && r2; r1 = r1->next, r2 = r2->next)
  {
    if ((res = regioncmp (r1, r2)) != 0) return (res);
  }

  if ((res = regioncmp (s1->regions, s2->regions)) != 0) return (res);

  return (0);
}

int
regioncmp (struct region *r1, struct region *r2)
{
  struct borderlist *h1, *h2;
  int res;

  /* primo criterio: numero di buchi */
  for (h1 = r1->border->next, h2 = r2->border->next;
       h1 && h2; h1 = h1->next, h2 = h2->next);
  if (h1) return (1);
  if (h2) return (-1);

  /* secondo/terzo criterio: confronto dei bordi esterni,
   * confronto lessicografico dei bordi dei buchi
   */

  for (h1 = r1->border, h2 = r2->border;
       h1 && h2; h1 = h1->next, h2 = h2->next)
  {
    if ((res = bordercmp (h1->sponda, h2->sponda)) != 0) return (res);
  }

  /* quarto criterio: valore di f (e' una informazione che non dipende
   * dalla rappresentazione)
   * (Negli esempi testati questo test non cambia nulla, ma in generale
   * non so.  Nel dubbio lo lascio).
   */

  if (r1->f > r2->f) return (1);
  if (r1->f < r2->f) return (-1);

  return (0);
}

/*
 * confronto di una lista di regioni
 */

int
subregionlistcmp (struct region *r1, struct region *r2, int *iperm1, int *iperm2)
{
  int res;

  if (r1 == 0 && r2 == 0) return (0);
  if (r2 == 0) return (-1);
  if (r1 == 0) return (1);

  if ( (res = subregioncmp (r1, r2, iperm1, iperm2)) != 0 ) return (res);
  return (subregionlistcmp (r1->next, r2->next, iperm1, iperm2));
}

/*
 * confronto regioni subordinato ad un ordinamento degli archi
 * basato sul loro tag
 */

int
subregioncmp (struct region *r1, struct region *r2, int *iperm1, int *iperm2)
{
  struct borderlist *h1, *h2;
  int res;

  /* primo criterio: numero di buchi */
  for (h1 = r1->border->next, h2 = r2->border->next;
       h1 && h2; h1 = h1->next, h2 = h2->next);
  if (h1) return (1);
  if (h2) return (-1);

  /* secondo/terzo criterio: confronto dei bordi esterni,
   * confronto lessicografico dei bordi dei buchi
   */

  for (h1 = r1->border, h2 = r2->border;
       h1 && h2; h1 = h1->next, h2 = h2->next)
  {
    if ((res = subbordercmp (h1->sponda, h2->sponda, iperm1, iperm2)) != 0) return (res);
  }

  /* quarto criterio: valore di f (e' una informazione che non dipende
   * dalla rappresentazione)
   * (Negli esempi testati questo test non cambia nulla, ma in generale
   * non so.  Nel dubbio lo lascio).
   */

  if (r1->f > r2->f) return (1);
  if (r1->f < r2->f) return (-1);

  return (0);
}

int
bordercmp (struct border *b1, struct border *b2)
{
  struct border *bp1, *bp2;
  int res;

  /* criterio zero: bordo vuoto vince (bordo esterno della regione esterna) */

  if (b1 == 0 && b2 == 0) return (0);

  if (b1 == 0 || b2 == 0)
  {
    if (b2 != 0) return (-1);
    return (1);
  }

  /* primo criterio, conto il numero di archi */

  assert (b1 && b2);
  bp1 = b1;
  bp2 = b2;

  do {
    bp1 = bp1->next;
    bp2 = bp2->next;
  } while (bp1 != b1 && bp2 != b2);

  if (bp1 != b1) return (1);
  if (bp2 != b2) return (-1);

  /* secondo criterio, confronto lessicografico tra i tratti di bordo */
  if ((res = singlebordercmp (b1, b2)) != 0) return (res);
  do {
    bp1 = bp1->next;
    bp2 = bp2->next;
    if ((res = singlebordercmp (bp1, bp2)) != 0) return (res);
  } while (bp1 != b1 && bp2 != b2);

  return (0);
}

/*
 * come prima, ma subordinatamente ai "tag" degli archi
 */

int
subbordercmp (struct border *b1, struct border *b2, int *iperm1, int *iperm2)
{
  struct border *bp1, *bp2;
  int res;

  /* criterio zero: bordo vuoto vince (bordo esterno della regione esterna) */

  if (b1 == 0 && b2 == 0) return (0);

  if (b1 == 0 || b2 == 0)
  {
    if (b2 != 0) return (-1);
    return (1);
  }

  /* primo criterio, conto il numero di archi */

  assert (b1 && b2);
  bp1 = b1;
  bp2 = b2;

  do {
    bp1 = bp1->next;
    bp2 = bp2->next;
  } while (bp1 != b1 && bp2 != b2);

  if (bp1 != b1) return (1);
  if (bp2 != b2) return (-1);

  /* secondo criterio, confronto lessicografico tra i tratti di bordo */
  if ((res = subsinglebordercmp (b1, b2, iperm1, iperm2)) != 0) return (res);
  do {
    bp1 = bp1->next;
    bp2 = bp2->next;
    if ((res = subsinglebordercmp (bp1, bp2, iperm1, iperm2)) != 0) return (res);
  } while (bp1 != b1 && bp2 != b2);

  return (0);
}

int
singlebordercmp (struct border *b1, struct border *b2)
{
  int res;

  /* primo criterio, confronto i due archi */
  if ((res = arccmp (b1->info, b2->info, 1)) != 0) return (res);

  /* se gli archi sono uguali guardo l'orientazione */
  if (b1->orientation == b2->orientation) return (0);
  if (b1->orientation > b2->orientation) return (1);
  return (-1);
}

/*
 * confronto subordinato all'ordinamento archi dato da "tag"
 */

int
subsinglebordercmp (struct border *b1, struct border *b2, int *iperm1, int *iperm2)
{
  int tag1, tag2;

  if (b1 == b2) return (0);
  tag1 = b1->info->tag;
  tag2 = b2->info->tag;
  if (iperm1) tag1 = iperm1[tag1];
  if (iperm2) tag2 = iperm2[tag2];

  /* primo criterio, confronto i due archi */
  if (tag1 != tag2) return (tag1 < tag2)?(-1):1;

  /* se gli archi sono uguali guardo l'orientazione */
  if (b1->orientation == b2->orientation) return (0);
  if (b1->orientation > b2->orientation) return (1);
  return (-1);
}

int
arccmp (struct arc *a1, struct arc *a2, int depth)
{
  int i, f1, f2;
  int res;

  /* short-cirtuit: if the pointer are equal they are equal! */
  if (a1 == a2) return (0);

  /* primo criterio, numero di estremi */
  if (a1->endpoints > a2->endpoints) return (1);
  if (a1->endpoints < a2->endpoints) return (-1);

  /* secondo criterio, numero di cuspidi */
  if (a1->cusps > a2->cusps) return (1);
  if (a1->cusps < a2->cusps) return (-1);

  /* terzo criterio, confronto lessicografico delle profondita' */

  assert (a1->dvalues == a2->dvalues);
  for (i = 0; i < a1->dvalues; i++)
  {
    if (a1->depths[i] > a2->depths[i]) return (1);
    if (a1->depths[i] < a2->depths[i]) return (-1);
  }

  /*
   * versions < 1.3.0 of appcontour corresponds to
   * depth = 0
   */

  if (depth-- <= 0) return (0);
  /* quarto criterio, confronto dei valori di f sulla sinistra */

  f1 = a1->regionleft->border->region->f;
  f2 = a2->regionleft->border->region->f;

  if (f1 < f2) return (-1);
  if (f1 > f2) return (1);

  if (depth-- <= 0) return (0);
  /* altri criteri? */
  res = borderarccmp (a1->regionleft, a2->regionleft, depth);
  if (res != 0) return (res);
  res = borderarccmp (a1->regionright, a2->regionright, depth);
  return (res);
}

int
borderarccmp (struct border *b1, struct border *b2, int depth)
{
  struct border *b1pt, *b2pt;
  int l1 = 0;
  int l2 = 0;
  int res;
  struct borderlist *bl1, *bl2;
  int isext1, isext2, isextr1, isextr2;

  /*
   * intrinsic comparison of two borders taking into account also
   * the starting pointer
   */

  /*
   * primo criterio: lunghezza del bordo
   */
  for (b1pt = b1->next; b1pt != b1; b1pt = b1pt->next) l1++;
  for (b2pt = b2->next; b2pt != b2; b2pt = b2pt->next) l2++;

  if (l1 != l2) return ( (l1 < l2)?(-1):1 );

  /*
   * secondo criterio: controllo quale sta nel bordo esterno
   */

  bl1 = b1->border;
  bl2 = b2->border;
  isext1 = (bl1->region->border == bl1);
  isext2 = (bl2->region->border == bl2);

  if (isext1 != isext2)
  {
    if (isext1) return (-1);
    return (1);
  }

  isextr1 = (bl1->region->border == 0);
  isextr2 = (bl2->region->border == 0);

  if (isextr1 != isextr2)
  {
    if (isextr1) return (-1);
    return (1);
  }

  if (depth-- <= 0) return (0);

  /*
   * controllo in profondita' su tutti gli archi
   */

  for (b1pt = b1->next, b2pt = b2->next; b1pt != b1; b1pt = b1pt->next, b2pt = b2pt->next)
  {
    if (b1pt->orientation != b2pt->orientation)
    {
      return ( (b1pt->orientation < 0)?(-1):1 );
    }
    res = arccmp (b1pt->info, b2pt->info, depth);
    if (res) return (res);
  }

  return (0);
}

int
count_hidden_arcs (struct sketch *s)
{
  struct arc *a;
  int k, count = 0;

  for (a = s->arcs; a; a = a->next)
  {
    if (a->dvalues <= 0) return (1000);
    for (k = 0; k < a->dvalues; k++)
    {
      count += a->depths[k];
    }
  }
  return (count);
}

/*
 * complete canonification process.  It is based
 * on a comparison among all subcanonifications
 * subordinated to every possible reorderings
 * of indistinguishable arcs
 */

/*
 * local prototypes
 */

int succ (int *ivec, int n);
int succvinc (int *ivec, int n, int *ipart, int partnum);

void
supercanonify (struct sketch *sketch)
{
  struct arc *arc;
  int i, *iperm, *ipermopt, *ipart;
  int res, count;
  int arcsnum = 0;
  int partnum = 0;
  int arccmpdepth = 5;
  struct sketch sketch2;

  if (sketch->arcs == 0) return;

  if (debug) printf ("canonify arcs...\n");
  if (sketch->isempty) return;
  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    canonifyarc (arc);
  }

  if (debug) printf ("sort arcs...\n");
  sortarcs (sketch, arccmpdepth);

  for (arc = sketch->arcs; arc; arc = arc->next) arcsnum++;

  iperm = (int *) malloc ( (arcsnum + 1) * sizeof(int) );
  ipermopt = (int *) malloc ( (arcsnum + 1) * sizeof(int) );
  ipart = (int *) malloc ( (arcsnum + 1) * sizeof(int) );

  for (i = 1; i <= arcsnum; i++) iperm[i] = ipermopt[i] = i;

  ipart[partnum] = 1;
  if (debug) printf ("computing the equivalence classes for the %d arcs\n", arcsnum);
  for (arc = sketch->arcs; arc->next; arc = arc->next)
  {
    if (arccmp (arc, arc->next, arccmpdepth) == 0)
    {
      ipart[partnum]++;
    } else {
      partnum++;
      ipart[partnum] = 1;
    }
  }
  partnum++;

  if (debug)
  {
    printf ("partitions:");
    for (i = 0; i < partnum; i++)
    {
      printf (" %d", ipart[i]);
    }
    printf ("\n");
  }
  /*
   * duplicate the regions structure of sketch
   */

  sketch2.arcs = sketch->arcs;
  sketch2.regions = dupregions (sketch->regions);
  if (debug > 2)
  {
    printf ("original sketch:\n");
    printsketch (sketch);
    printf ("duplicated sketch:\n");
    printsketch (&sketch2);
  }

  /*
   * cycle on all permutations consistent with the partitioning
   * warning: vector iperm is used from index 1 to arcsnum!
   */

  count = 1;
  while (succvinc (iperm, arcsnum, ipart, partnum))
  {
    if (debug) printf ("permutation number %d\n", count);
    count++;
    subcanonify (&sketch2, iperm);
    res = subregionlistcmp (sketch2.regions, sketch->regions, iperm, ipermopt);
    if (res < 0)
    {
      /* new optimal region description */
      for (i = 1; i <= arcsnum; i++) ipermopt[i] = iperm[i];
      subcanonify (sketch, ipermopt);
    }
  }
  if (debug) printf ("checked %d consistent permutations\n", count);

  remregions (sketch2.regions);

  if (debug) printf ("removed duplicated regions\n");

  /*
   * retagging regions and reordering
   */

  for (arc = sketch->arcs; arc; arc = arc->next) arc->tag = ipermopt[arc->tag];
  sketch->arcs = reorderbytags (sketch->arcs);
}

struct arc *
reorderbytags (struct arc *a)
{
  struct arc *anext, *apt;

  if (a == 0 || a->next == 0) return (a);
  anext = reorderbytags (a->next);

  if (a->tag < anext->tag)
  {
    a->next = anext;
    return (a);
  }
  for (apt = anext; apt; apt = apt->next)
  {
    if (apt->next == 0)
    {
      apt->next = a;
      a->next = 0;
      return (anext);
    }
    if (a->tag < apt->next->tag)
    {
      a->next = apt->next;
      apt->next = a;
      return (anext);
    }
  }
  assert (0);
  return (0);  
}

int
succvinc (int *ivec, int n, int *ipart, int partnum)
{
  if (partnum <= 1) return succ (ivec, n);
  if (succvinc (ivec + ipart[0], n - ipart[0], ipart + 1, partnum - 1)) return (1);
  if (succ(ivec, ipart[0]) == 0) return (0);
  return (1);
}

int
succ (int *ivec, int n)
{
  int k, newval;
  int lower, upper, middle;

  if (n <= 1) return (0);
  if (succ(ivec+1, n-1)) return (1);
  /* sono finite le permutazioni che non toccano
   * il primo elemento, in particolare viene ripristinato
   * l'ordine iniziale, ora aumento il primo elemento
   */
  if (ivec[1] > ivec[n])
  {
    /* ho finito... ripristino la permutazione iniziale */
    newval = ivec[1];
    for (k = 1; k < n; k++) ivec[k] = ivec[k+1];
    ivec[n] = newval;
    return (0);
  }
  lower = 1;
  upper = n;
  while (upper - lower > 1)
  {
    middle = (lower + upper)/2;
    if (ivec[middle] > ivec[1])
    {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  newval = ivec[upper];
  ivec[upper] = ivec[1];
  ivec[1] = newval;
  return (1);
}

/*
 * end of cycle on permutations functions
 */

struct region *
dupregions (struct region *r)
{
  struct region *newr;

  if (r == 0) return (0);

  newr = (struct region *) malloc (sizeof (struct region));
  newr->f = r->f;
  newr->tag = r->tag;
  newr->strati = 0;
  newr->border = dupborderlist (r->border, newr);
  newr->next = dupregions (r->next);
  return (newr);
}

struct borderlist *
dupborderlist (struct borderlist *bl, struct region *r)
{
  struct borderlist *newbl;

  if (bl == 0) return (0);

  newbl = (struct borderlist *) malloc (sizeof (struct borderlist));
  newbl->region = r;
  newbl->sponda = dupborder (bl->sponda, newbl);
  newbl->next = dupborderlist (bl->next, r);
  return (newbl);
}

struct border *
dupborder (struct border *b, struct borderlist *bl)
{
  struct border *bnext, *newbnext, *bpt;

  if (b == 0) return (0);

  bnext = b->next;
  b->next = 0;       /* opening the circular list */
  newbnext = dupborder_opened (bnext, bl);
  b->next = bnext;   /* reclosing original circular list */
  for (bpt = newbnext; bpt; bpt = bpt->next)
  {
    if (bpt->next == 0)
    {
      bpt->next = newbnext;
      return (bpt);
    }
  }
  assert (0);
  return (0);
}

struct border *
dupborder_opened (struct border *b, struct borderlist *bl)
{
  struct border *newb;

  if (b == 0) return (0);

  newb = (struct border *) malloc (sizeof (struct border ));
  newb->border = bl;
  newb->info = b->info;
  newb->orientation = b->orientation;
  newb->next = dupborder_opened (b->next, bl);
  return (newb);
}

void
remregions (struct region *rl)
{
  struct region *r;
  struct borderlist *bl;
  struct border *bstart, *b;

  /*
   * first we must remove all pointers to the arc structures
   */

  for (r = rl; r; r = r->next)
  {
    if (debug > 1) printf ("region %d, zeroing info pointers\n", r->tag);
    for (bl = r->border; bl; bl = bl->next)
    {
      if (bl->sponda == 0) continue;  /* external region */
      bstart = bl->sponda;
      bstart->info = 0;
      
      for (b = bstart->next; b != bstart; b = b->next)
      {
        b->info = 0;
      }
    }
  }

  if (debug) printf ("zeroed all info pointers\n");

  /*
   * now we can use freeregion
   */

  freeregion (rl);
}

/*
 * canonification process subordinate to a given
 * ordering of the arcs.
 * The ordering of the arcs is induced by their "tag"
 * value, if "iperm" is nonnull, then it works as a
 * permutation vector on the "tag" value:
 *   newtag = iperm[tag]
 *
 */

void
subcanonify (struct sketch *sketch, int *iperm)
{
  struct region *r;
  struct borderlist *h;
  int tag;

  if (debug) printf ("subcanonify region borders...\n");
  for (r = sketch->regions; r; r = r->next)
  {
    for (h = r->border; h; h = h->next)
    {
      if (h->sponda) h->sponda = subcanonifyborder (h->sponda, iperm);
    }
  }

  if (debug) printf ("sort holes for each region...\n");
  for (r = sketch->regions; r; r = r->next)
  {
    r->border->next = subsortholelist (r->border->next, iperm);
  }

  if (debug) printf ("sort bounded regions...\n");
  if (sketch->regions->next)
  {
    sketch->regions->next = subsortregionlist (sketch->regions->next, iperm);
    for (tag = 0, r = sketch->regions; r; r = r->next)
    {
      r->tag = tag++;
    }
  }
}

void
canonify (struct sketch *sketch)
{
  struct arc *arc;
  struct region *r;
  struct borderlist *h;
  int tag;

  if (debug) printf ("canonify arcs...\n");
  if (sketch->isempty) return;
  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    canonifyarc (arc);
  }

  if (debug) printf ("sort arcs...\n");
  sortarcs (sketch, 1);

  if (debug) printf ("canonify region borders...\n");
  for (r = sketch->regions; r; r = r->next)
  {
    for (h = r->border; h; h = h->next)
    {
      if (h->sponda) h->sponda = canonifyborder (h->sponda);
    }
  }

  if (debug) printf ("sort holes for each region...\n");
  for (r = sketch->regions; r; r = r->next)
  {
    r->border->next = sortholelist (r->border->next);
  }

  if (debug) printf ("sort bounded regions...\n");
  if (sketch->regions->next)
  {
    sketch->regions->next = sortregionlist (sketch->regions->next);
    for (tag = 0, r = sketch->regions; r; r = r->next)
    {
      r->tag = tag++;
    }
  }

  sketch->arcs = sortequivarcs (sketch->arcs);
  /* rinumero gli archi */
  tag = 1;
  for (arc = sketch->arcs; arc; arc = arc->next) arc->tag = tag++;
}

void
canonifyarc (struct arc *arc)
{
  static int buf[200];
  int i, j, n, iopt;

  if (debug) printf ("looking for loops with cusps...\n");
  if (arc->endpoints > 0) return;
  if (arc->cusps < 1) return;
  //if (arc->dvalues < 1) return;
  if (debug) printf ("arc %d with %d endpoints and %d d values\n",
    arc->tag, arc->endpoints, arc->dvalues);
  if (arc->dvalues > 97) {fprintf (stderr, "too many cusps..."); return;}
  n = arc->dvalues;
  if (debug) printf ("rotating modulo %d\n", n);
  for (i = 0; i < n; i++) buf[i] = buf[i+n] = arc->depths[i];
  iopt = 0;
  for (i = 1; i < n; i++)
  { /* confronta la rotazione i con la rotazione iopt
     * in modo lessicografico 
     */
    for (j = 0; j < n; j++)
    {
      if (buf[i+j] > buf[iopt+j])   /* nuova rotazione peggiore */
        break;
      if (buf[i+j] < buf[iopt+j])   /* nuova rotazione migliore */
      {
        iopt = i;
        break;
      }
    }
  }
  for (i = 0; i <= arc->dvalues; i++) arc->depths[i] = buf[iopt + i];
  /* ricorda che il primo e ultimo valore DEVONO essere duplicati
   * per archi senza estremi
   */
}

/*
 * un bordo di una regione (o di un buco di una regione)
 * e' definito a meno del punto di partenza, quindi si tratta
 * di determinare il punto di partenza ottimale che renda
 * il bordo risultante lessicograficamente vantaggioso 
 * (in termini della funzione bordercmp che confronta due
 * bordi tra di loro)
 * notasi che il confronto tra bordi dipende anche dall'esito del
 * confronto tra archi, quindi e' fondamentale aver prima canonificato
 * gli archi senza "endpoints"
 */

struct border *
canonifyborder (struct border *b)
{
  struct border *bp, *bopt;

  if (b->next == b) return (b);  /* una sola componente! */

  /* l'idea e' di chiamare ripetutamente la funzione di
   * confronto tra bordi passando via via puntatori in punti
   * diversi, alla ricerca del puntatore che fornisce il 
   * risultato ottimale
   */

  bopt = b;
  bp = b->next;
  do {
    if (bordercmp (bp, bopt) < 0)
      /* ho trovato un aggancio migliore */
      bopt = bp;
    bp = bp->next; 
  } while (bp != b);
  return (bopt);
}

/*
 * la canonificazione del bordo subordinata ad un ordinamento degli
 * archi indotto dal loro "tag" (o una permutazione)
 */

struct border *
subcanonifyborder (struct border *b, int *iperm)
{
  struct border *bp, *bopt;

  if (b->next == b) return (b);  /* una sola componente! */

  /* l'idea e' di chiamare ripetutamente la funzione di
   * confronto tra bordi passando via via puntatori in punti
   * diversi, alla ricerca del puntatore che fornisce il 
   * risultato ottimale
   */

  bopt = b;
  bp = b->next;
  do {
    if (subbordercmp (bp, bopt, iperm, iperm) < 0)
      /* ho trovato un aggancio migliore */
      bopt = bp;
    bp = bp->next; 
  } while (bp != b);
  return (bopt);
}

/*
 * per comodita' di visualizzazione conviene riordinare
 * gli archi (e rinumerarli) in ordine di complessita'
 * crescente.  Questo non ha un reale effetto sulle funzioni
 * di confronto, ma e' consigliabile per un veloce confronto
 * visivo da parte dell'utente.
 * nota che l'ordinamento e' indefinito per archi non
 * distinguibili in base al numero di estremi e di cuspidi e 
 * ai valori delle profondita'
 */

void
sortarcs (struct sketch *s, int depth)
{
  struct arc *arc;
  int tag = 1;

  if (s->isempty) return;
 /*
  * usiamo un metodo naive con complessita' n^2!
  */
  s->arcs = sortarclist (s->arcs, depth);

  for (arc = s->arcs; arc; arc = arc->next)
    arc->tag = tag++;
}

struct arc *
sortarclist (struct arc *arc, int depth)
{
  struct arc *a, *aopt, *aprev, *aoptprev;

  assert (arc);
  if (arc->next == 0) return (arc);

  aopt = arc;
  aoptprev = 0;
  aprev = arc;
  for (a = arc->next; a; a = a->next)
  {
    if (arccmp (a, aopt, depth) < 0) {aopt = a; aoptprev = aprev;}
    aprev = a;
  }

  if (debug) printf ("arco ottimale trovato: %d\n", aopt->tag);
  if (aopt == arc)
  {
    arc->next = sortarclist (arc->next, depth);
    return (arc);
  }
  /* altrimenti aoptprev punta ad aopt */
  assert (aoptprev);
  aoptprev->next = aopt->next;     /* rimuovo aopt dalla lista */
  aopt->next = sortarclist (arc, depth);  /* inserisco aopt in testa  */
  return (aopt);
}

/*
 * ora vogliamo definire un criterio per riordinare
 * tra di loro gli archi indistinguibili.  Scegliamo
 * di guardare quale arco compare per primo nella
 * descrizione delle regioni.  Quindi la canonificazione
 * delle regioni deve essere effettuata prima!
 */

struct arc *
sortequivarcs (struct arc *arc)
{
 /*
  * usiamo un metodo naive con complessita' n^2!
  */

  if (arc->next == 0) return (arc);        /* ultimo arco, ovviamente e' ordinato */
  arc->next = sortequivarcs (arc->next);   /* ordino tutti i successivi */
  return (mergeequivarcs (arc, arc->next));
}

struct arc *
mergeequivarcs (struct arc *arc, struct arc *rest)
{
  struct border *b1, *b2, *b;
  struct borderlist *bl1, *bl2, *bl;
  struct region *r1, *r2;

  arc->next = rest;       /* nell'eventualita' che arc < rest */
  if (rest == 0) return (arc);
  if (arccmp (arc, rest, 1) != 0) return (arc);
                   /* l'arco si distingue dai rimanenti, non lo muovo */
  /* ora devo confrontare arc e rest e vedere quale compare prima con
   * orientazione negativa nell'elenco delle regioni.  Non e' difficile
   * poiche' le regioni interessate sono quelle alla destra dell'arco
   */
  b1 = arc->regionright;
  b2 = rest->regionright;
  assert (b1 != b2);
  bl1 = b1->border;
  bl2 = b2->border;
  r1 = bl1->region;
  r2 = bl2->region;
  if (r1->tag < r2->tag) return (arc);
  if (r1->tag == r2->tag)
  {
    /* devo controllare su quale componente connessa stanno */
    if (bl1 != bl2)
    {
      for (bl = r1->border; bl; bl = bl->next)
      {
        if (bl == bl1) return (arc);
        if (bl == bl2) break;
      }
    } else {
      /* stanno nella stessa componente connessa */
      for (b = bl1->sponda;; b = b->next)
      {
        if (b == b1) return (arc);
        if (b == b2) break;
      }
    }
  }
  rest->next = mergeequivarcs (arc, rest->next);
  return (rest);
}

struct borderlist *
sortholelist (struct borderlist *hl)
{
  struct borderlist *h, *hopt, *hprev, *hoptprev;

  if (hl == 0 || hl->next == 0) return (hl);

  hopt = hl;
  hoptprev = 0;
  hprev = hl;
  for (h = hl->next; h; h = h->next)
  {
    if (bordercmp (h->sponda, hopt->sponda) < 0) 
    {hopt = h; hoptprev = hprev;}
    hprev = h;
  }

  if (hopt == hl)
  {
    hl->next = sortholelist (hl->next);
    return (hl);
  }
  /* altrimenti hoptprev punta ad hopt */
  assert (hoptprev);
  hoptprev->next = hopt->next;     /* rimuovo aopt dalla lista */
  hopt->next = sortholelist (hl);  /* inserisco aopt in testa  */
  return (hopt);
}

/*
 * sorting subordinato ad un ordinamento degli archi basato sul tag
 */

struct borderlist *
subsortholelist (struct borderlist *hl, int *iperm)
{
  struct borderlist *h, *hopt, *hprev, *hoptprev;

  if (hl == 0 || hl->next == 0) return (hl);

  hopt = hl;
  hoptprev = 0;
  hprev = hl;
  for (h = hl->next; h; h = h->next)
  {
    if (subbordercmp (h->sponda, hopt->sponda, iperm, iperm) < 0) 
    {hopt = h; hoptprev = hprev;}
    hprev = h;
  }

  if (hopt == hl)
  {
    hl->next = subsortholelist (hl->next, iperm);
    return (hl);
  }
  /* altrimenti hoptprev punta ad hopt */
  assert (hoptprev);
  hoptprev->next = hopt->next;     /* rimuovo aopt dalla lista */
  hopt->next = subsortholelist (hl, iperm);  /* inserisco aopt in testa  */
  return (hopt);
}

struct region *
sortregionlist (struct region *region)
{
  struct region *r, *ropt, *rprev, *roptprev;

  if (region == 0 || region->next == 0) return (region);

  ropt = region;
  roptprev = 0;
  rprev = region;
  for (r = region->next; r; r = r->next)
  {
    if (regioncmp (r, ropt) < 0) 
    {ropt = r; roptprev = rprev;}
    rprev = r;
  }

  if (ropt == region)
  {
    region->next = sortregionlist (region->next);
    return (region);
  }
  /* altrimenti hoptprev punta ad hopt */
  assert (roptprev);
  roptprev->next = ropt->next;           /* rimuovo ropt dalla lista */
  ropt->next = sortregionlist (region);  /* inserisco ropt in testa  */
  return (ropt);
}

/*
 * riordinamento subordinato ad un ordinamento degli archi
 * basato sul loro "tag"
 */

struct region *
subsortregionlist (struct region *region, int *iperm)
{
  struct region *r, *ropt, *rprev, *roptprev;

  if (region == 0 || region->next == 0) return (region);

  ropt = region;
  roptprev = 0;
  rprev = region;
  for (r = region->next; r; r = r->next)
  {
    if (subregioncmp (r, ropt, iperm, iperm) < 0) 
    {ropt = r; roptprev = rprev;}
    rprev = r;
  }

  if (ropt == region)
  {
    region->next = subsortregionlist (region->next, iperm);
    return (region);
  }
  /* altrimenti hoptprev punta ad hopt */
  assert (roptprev);
  roptprev->next = ropt->next;           /* rimuovo ropt dalla lista */
  ropt->next = subsortregionlist (region, iperm);  /* inserisco ropt in testa  */
  return (ropt);
}

void
postprocesssketch (struct sketch *sketch)
{
  struct region *region, *extregion;
  struct borderlist *hole;
  struct arc *arc;
  int tag;
  extern int dorecomputef;
  extern int doretagregions;
  extern int finfinity;

  /*
   * sfortunatamente free_connected_components distrugge l'informazione
   * sulle componenti connesse precedentemente eliminate durante un
   * extractcc, quindi non e' il caso di liberare gli strati.
   * d'altronde le sistemazioni di postprocess non dovrebbero influire
   * negativamente sulla struttura topologica del manifold 3D
   */
  //free_connected_components (sketch);
  if (debug) printf ("1: porta la regione esterna in prima posizione\n");

  if (sketch->regions->border->sponda)
  {
    /* la regione esterna non si trova in prima posizione, devo
     * cercarla e portarla in prima posizione
     */
    assert (sketch->regions->next);
    for (region = sketch->regions; region->next; region = region->next)
    {
      if (region->next->border->sponda == 0)
      {
        extregion = region->next;
        region->next = extregion->next;
        extregion->next = sketch->regions;
        sketch->regions = extregion;
        break;
      }
    }
  }

  if (debug && doretagregions) printf ("2: rinumero gli archi e le regioni\n");
  if (debug && !doretagregions) printf ("2: rinumero gli archi\n");

  for (tag = 0, region = sketch->regions; region; tag++, region = region->next)
          if (doretagregions) region->tag = tag;
  sketch->regioncount = tag;
  for (tag = 1, arc = sketch->arcs; arc; arc = arc->next)
          arc->tag = tag++;
  sketch->arccount = tag - 1;

  if (debug) printf ("3: definisco regionleft e regionright per gli archi\n");

  for (region = sketch->regions; region; region = region->next)
  {
    for (hole = region->border; hole; hole = hole->next)
    {
      if (hole->sponda) defineregionleftright (hole->sponda);
    }
  }

  if (debug) printf ("4: conto gli estremi di ciascun arco\n");

  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    arc->endpoints = 0;
    if (arc->regionleft->next != arc->regionleft) arc->endpoints++;
    if (arc->regionright->next != arc->regionright) arc->endpoints++;
    if (sketch->huffman_labelling && arc->cusps != arc->depthsdim - 1)
      fprintf (stderr, "depthsdim not equal cusps + 1, this is not a problem\n");
    if (arc->endpoints == 0 && arc->depths)
    { /* in questo caso il numero di cuspidi e' uguale al numero
       * di valori di d, a meno che non siano zero.  
       * L'utente deve pero' assegnare un valore
       * ulteriore di d, uguale al primo...
       */
      arc->dvalues = arc->cusps;
      if (arc->cusps == 0) arc->dvalues = 1;
      if (arc->depths[arc->cusps] != arc->depths[0])
      {
        fprintf (stderr, "Warning: the last value of d for a loop (arc %d) MUST duplicate\n", arc->tag);
        fprintf (stderr, "the first one; adjusting...\n");
      }
      arc->depths[arc->cusps] = arc->depths[0];
    } else {
      // arc->cusps = arc->depthsdim - 1;
      arc->dvalues = arc->cusps + 1;
    }
  }

  if (debug && dorecomputef) printf ("5: definizione dei valori di f\n");

  if (dorecomputef)
  {
    extregion = sketch->regions;
    if (sketch->extregion) extregion = sketch->extregion;
    assert (extregion->border->sponda == 0);
    computefvalue (sketch, extregion, finfinity);
  }

  if (debug) printf ("6: controllo correttezza bordo esterno\n");

  assert (adjust_isexternalinfo (sketch) == 0);
}

void
computefvalue (struct sketch *s, struct region *extregion, int finfinity)
{
  int goon, fleft, fright;
  struct region *r;
  struct arc *a;

  for (r = s->regions; r; r = r->next)
  {
    r->f = F_UNDEF;
  }
  assert (extregion->border->sponda == 0); /* this MUST be the external region */
  extregion->f = finfinity;    /* valore 0 nella regione esterna */
  goon = 1;
  while (goon)
  {
    goon = 0;
    for (a = s->arcs; a; a = a->next)
    {
      fleft = a->regionleft->border->region->f;
      fright = a->regionright->border->region->f;
      if (fleft == F_UNDEF && fright == F_UNDEF) continue;
      if (fleft != F_UNDEF && fright != F_UNDEF) continue;
      goon = 1;
      if (fright != F_UNDEF) a->regionleft->border->region->f = fright + 2;
        else a->regionright->border->region->f = fleft - 2;
    }
  }
}

/*
 * in alcune circostanze (dopo particolari operazioni topologiche)
 * viene persa l'informazione su quale componente connessa del bordo
 * di una regione sia quella esterna.
 * questa funzione serve a ripristinare tale informazione e restituisce
 * il numero di modifiche che ha dovuto fare rispetto all'informazione
 * (non affidabile) preesistente.
 */

#define ISEXTERNAL_UNDEF (-1)

int
adjust_isexternalinfo (struct sketch *sketch)
{
  struct region *r;
  struct borderlist *bl, *blprev;
  int goon, changes, found;

  found = 0;
  for (r = sketch->regions; r; r = r->next)
  {
    for (bl = r->border; bl; bl = bl->next) 
    {
      if (bl->sponda == 0) {bl->isexternal = 1; found++; sketch->extregion = r;}
        else bl->isexternal = ISEXTERNAL_UNDEF;
    }
  }
  assert (found == 1);

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (r = sketch->regions; r; r = r->next)
    {
      if (iei_process_region (r) > 0) goon = 1;
    }
  }

  changes = 0;
  /* controllo di aver trovato tutte le informazioni */
  for (r = sketch->regions; r; r = r->next)
  {
    if (r->border->isexternal == 1) continue;
    blprev = r->border;
    found = 0;
    for (bl = r->border->next; bl; bl = bl->next)
    {
      if (bl->isexternal == 1)
      {
        changes++;           /* devo scambiare i bordi */
        blprev->next = bl->next;
        bl->next = r->border;
        r->border = bl;
        found = 1;
        break;
      }
      blprev = bl;
    }
    assert (found == 1);
  }
  assert (sketch->extregion->border->sponda == 0);
  return (changes);
}

int
iei_process_region (struct region *region)
{
  int found = 0;
  int changes = 0;
  struct borderlist *bl;
  struct border *bp, *btrans;

  for (bl = region->border; bl; bl = bl->next)
    if (bl->isexternal == 1) found++;

  assert (found <= 1);
  if (! found) return (0);

  for (bl = region->border; bl; bl = bl->next)
  {
    if (bl->isexternal == ISEXTERNAL_UNDEF)
    {
      changes++;
      bl->isexternal = 0;
    }
    if ((bp = bl->sponda) == 0) continue;
    do {
      btrans = gettransborder(bp);
      if (btrans->border->isexternal == ISEXTERNAL_UNDEF)
      {
        changes++;
        btrans->border->isexternal = 1;
      }
      bp = bp->next;
    } while (bp != bl->sponda);
  }
  return (changes);
}

void
defineregionleftright (struct border *border)
{
  struct border *b;
  struct arc *arc;

  b = border;
  do
  {
    arc = b->info;
    if (abs(b->orientation) == 2)
    {
      fprintf (stderr, "Implicitly orienting arc %d\n", arc->tag);
      b->orientation /= 2;
    }
    if (arc->depths == 0)
    {
      fprintf (stderr, "Dummy definition of depth for arc %d\n", arc->tag);
      arc->depths = (int *) malloc (0);
      arc->depthsdim = 0;
      arc->cusps = 0;
    }
    if (b->orientation > 0) arc->regionleft = b;
      else arc->regionright = b;
    b = b->next;
  } while (b != border);
}

void
printsketch (struct sketch *sketch)
{
  struct arc *a;
  struct region *r;
  struct borderlist *rl;
  int i, notfirst;
  char chleft, chright;

  if (sketch == 0) {printf ("sketch {}\n"); return;}
  printf ("sketch {\n");
  for (a = sketch->arcs; a; a = a->next)
  {
    chleft = chright = ' ';
    printf ("Arc %d: ", a->tag);
    switch (a->endpoints)
    {
      case 0:
        chleft = '(';
        chright = ')';
        break;
      case 1:
        chleft = '[';
        chright = ')';
        break;
      case 2:
        chleft = '[';
        chright = ']';
        break;
    }

    if (a->depths)
    {
      assert (a->depthsdim >= 0);
      printf ("%c", chleft);
      notfirst = 0;
      if (a->depthsdim > 0)
      {
        for (i = 0; i <= a->cusps; i++)
        {
          if (notfirst) printf (" ");
          notfirst = 1;
          printf ("%d", a->depths[i]);
        }
      } else {
        for (i = 0; i < a->cusps; i++)
        {
          if (notfirst) printf (" ");
          notfirst = 1;
          printf ("c");
        }
      }
      printf ("%c;\n", chright);
    } else printf (" [no information]\n");
  }
  for (r = sketch->regions; r; r = r->next)
  {
    printf ("Region %d ", r->tag);
    if (r->f != F_UNDEF) printf ("(f =%2d):", r->f);
           else          printf ("         :");
    for (rl = r->border; rl; rl = rl->next)
    {
      if (rl->sponda) printborder (rl->sponda, r);
        else printf (" ()");
    }
    printf (";\n");
  }
  printf ("}\n");
}

void
printborder (struct border *b, struct region *r)
{
  int ori, notfirst, tag;

  struct border *bp;
  bp = b;
  notfirst = 0;
  printf (" (");
  do
  {
    assert (bp->border->region == r);
    ori = bp->orientation;
    if (abs(ori) == 2)
    {
      fprintf (stderr, "warning: implicit orientation\n");
      ori /= 2;
    }
    assert (ori == 1 || ori == -1);
    if (notfirst) printf (" ");
    notfirst = 1;
    tag = 9999;
    if (bp->info) tag = bp->info->tag;
    printf ("%ca%d", (ori == 1)?'+':'-', tag);
    bp = bp->next;
  } while (bp != b);
  printf (")");
}

/*
 * compute the change in the d value across the
 * node ahead (ori > 0) or behind (ori < 0)
 */

int
get_d_increase_across_node (struct arc *arc, int ori)
{
  struct border *b;
  struct arc *arcnext;
  int dhere, dthere;

  dhere = arc->depths[0];
  if (ori > 0) dhere = arc->depths[arc->dvalues-1];
  b = arc->regionright;
  if (ori > 0) b = arc->regionleft;
  b = gettransborder (b->next);
  arcnext = b->next->info;
  dthere = arcnext->depths[0];
  if (ori < 0) dthere = arcnext->depths[arcnext->dvalues-1];
  if (debug > 1) printf ("arc %d, arcnext %d, increase %d\n", 
                arc->tag, arcnext->tag, dthere - dhere);
  return (dthere - dhere);
}

/*
 * restituisce il bordo che corrisponde a b
 * ma visto dalla regione adiacente
 */

struct border *
gettransborder (struct border *b)
{
  struct arc *arc;

  arc = b->info;
  if (arc->regionleft == b) return (arc->regionright);
  return (arc->regionleft);
}

/*
 * restituisce il precedente della lista circolare
 */

struct border *
prevborder (struct border *b)
{
  struct border *bp;

  bp = b;
  do {
    if (bp->next == b) return (bp);
    bp = bp->next;
  } while (bp != b);
  assert (0);
  return (0);
}

/*
 * b1 e b2 stanno sullo stesso s1?
 */

int
findinborder (struct border *b1, struct border *b2)
{
  struct border *b;

  if (b1 == b2) return (1);      /* ovvio */
  if (b1->next == b1 || b2->next == b2) return (0);
                                 /* un bordo si riduce ad un solo elemento */

  for (b = b1->next; b != b1; b = b->next)
  {
    if (b == b2) return (1);
  }
  return (0);
}

/*
 * incollo due archi consecutivi
 */

struct arc *
mergearcs (struct arc *arc1, struct arc *arc2, struct sketch *sketch)
{
  return (mergearcsc (arc1, arc2, 0, sketch));
}

struct arc *
mergearcsc (struct arc *arc1, struct arc *arc2, int dincr,
            struct sketch *sketch)
{
  int newdim, i, sign;
  int *newdepths;
  struct border *bl, *br;

  assert (arc1->regionleft->info == arc1);
  assert (arc2->regionleft->info == arc2);
  assert (arc1->regionright->info == arc1);
  assert (arc2->regionright->info == arc2);
  arc1->regionleft->info = 0;
  arc2->regionleft->info = 0;
  arc1->regionright->info = 0;
  arc2->regionright->info = 0;
  sign = 1;
  if (dincr < 0) {sign = -1; dincr = -dincr;}
  if (arc1 == arc2)   /* struttura a goccia */
  {
    if (debug) printf ("struttura a goccia\n");
    assert (arc1->endpoints == 1);
    arc1->endpoints = 0;
    if (sketch->huffman_labelling)
      assert (arc1->depths[0] == arc1->depths[arc1->dvalues-1] + sign*dincr);
    if (dincr != 0)          /* devo aggiungere un po' di cuspidi */
    {
      if (sketch->huffman_labelling) {
        arc1->depthsdim = arc1->cusps + 1 + dincr;
        newdepths = (int *) malloc ((arc1->depthsdim)*sizeof(int));
        for (i = 0; i <= arc1->cusps; i++) newdepths[i] = arc1->depths[i];
        for (i = 0; i < dincr; i++) newdepths[arc1->cusps + 1 + i] = 
                         newdepths[arc1->cusps + i] + sign;
        free (arc1->depths);
        arc1->depths = newdepths;
      }
      arc1->cusps += dincr;
      arc1->dvalues += dincr;
    }
    arc1->dvalues--;
    assert (arc1->dvalues == arc1->cusps);
    return (arc1);
  }

  /* due archi diversi */
  /* mi assicuro che le sponde non siano ancora eliminate! */
  assert (arc1->endpoints > 0 && arc2->endpoints > 0);
  if (sketch->huffman_labelling) {
    assert (arc1->depths[arc1->dvalues-1] + sign*dincr == arc2->depths[0]);
    newdim = arc1->cusps + arc2->cusps + 1 + dincr;
    newdepths = (int *) malloc (newdim * sizeof (int));
    for (i = 0; i < arc1->dvalues; i++) newdepths[i] = arc1->depths[i];
    for (i = 0; i < dincr; i++) 
      newdepths[i + arc1->dvalues] = newdepths[i + arc1->dvalues - 1] + sign;
    for (i = 0; i < arc2->dvalues; i++) 
      newdepths[i + arc1->dvalues + dincr - 1] = arc2->depths[i];
    arc1->depthsdim = newdim;
    free (arc1->depths);
    arc1->depths = newdepths;
  }
  arc1->cusps += arc2->cusps + dincr;
  arc1->dvalues += arc2->dvalues + dincr - 1;
  arc2->regionleft = arc2->regionright = 0;
  removearc (arc2, sketch);

  /* ora devo definire endpoints... */
  bl = arc1->regionleft;
  br = arc1->regionright;
  arc1->endpoints = 2;
  if (bl == bl->next->next || br == br->next->next) arc1->endpoints = 1;
  assert (bl != bl->next || br != br->next);
  return (arc1);
}

/*
 * operazione taglia/incolla riferita a bordi di regione
 * opera a livello piu' alto di topo_change, ma ha bisogno
 * di sapere se si tratta di erosione o di dilatazione
 */


int
topo_change_g (struct border *b1, struct border *b2, int type, 
               struct sketch *sketch)
{
  extern int heisemberg;
  struct region *r1, *r2, *newr;
  struct borderlist *bl1, *bl2, *newbl;
  int res = 0;

  bl1 = b1->border;
  bl2 = b2->border;
  r1 = bl1->region;
  r2 = bl2->region;

  if (type == TC_UNKNOWN)   /* cerco di capire in che caso sono */
  {
    if (r1 != r2) type = TC_DILATION;
  }
  assert (type == TC_EROSION || type == TC_DILATION);

  if ((type == TC_EROSION && bl1 != bl2) ||
      (type == TC_DILATION && bl1 == bl2))      /* 1 regione -> 1 regione */
  {
    return (topo_change (b1, b2));
  }

  if (type == TC_EROSION)                       /* 1 regione -> 2 regioni */
  {
    assert (bl1 == bl2);
    /* DOVE METTO I BUCHI? non c'e' modo di decidere in base alle sole informazioni
     * topologiche in modo sensato
     * qualunque scelta porta ad una configurazione sensata, quindi qui opero una
     * scelta arbitraria, consistente nel lascare tutti i buchi nella regione
     * originale, mentre b2 sta sulla nuova regione, che non ha buchi
     */
    if (bl1->region->border->next && heisemberg < 0)
    {   /* ci sono buchi, avviso l'utente! */
      fprintf (stderr, "Warning: Heisenberg rules!\n");
      fprintf (stderr, "  splitting a region with holes requires a degree of ");
      fprintf (stderr, "arbitraryness,\n  I will put all holes in one of the regions.\n");
      fprintf (stderr, "  Use option \"--ti <int>\" or C2::<int> (e.g.) to control ");
      fprintf (stderr, "the behaviour\n");
    }
    res = topo_change (b1, b2);
    /* e' stata creato un nuovo bordo per b2, che devo recuperare */
    newbl = b2->border;
    assert (b1->border == bl1 && newbl != bl1);
    newbl = extractborderlist (newbl);
    newr = newregion (sketch);
    newr->border = newbl;
    newbl->region = newr;
    if (heisemberg >= 0) /* move "holes" according to heisemberg */
    {
      transfer_isles (bl1, newbl, heisemberg);
    }
  } else {           /* type = TC_DILATION:   2 regioni -> 1 regione */
    assert (bl1 != bl2);
    assert (r1 == regionunion (r1, r2, sketch));
    res = topo_change (b1, b2);
  }
  return (res);
}

/*
 * operazione taglia/incolla sui bordi delle regioni
 * si vuole cambiare la topologia in modo che il 
 * successivo di b1 diventi il successivo di b2
 * e viceversa.
 * prerequisito: b1 e b2 devono appartenere alla stessa
 * regione.
 * Viene anche sistemata la lista delle componenti
 * connesse in modo opportuno (aggiunta o eliminazione
 * di una componente connessa.
 */

int
topo_change (struct border *b1, struct border *b2)
{
  struct borderlist *bl1, *bl2, *newbl;
  struct region *r;

  assert (b1 != b2);
  bl1 = b1->border;
  bl2 = b2->border;
  r = bl1->region;
  assert (r == bl2->region);

  if (bl1 == bl2)   /* b1 e b2 stanno nella stessa   */
                                  /* componente connessa del bordo */
  {
    newbl = newborderlist (r);
    topo_change_l (b1, b2);
    redefineborder (b2, newbl);
    bl1->sponda = b1;
    newbl->sponda = b2;
  } else {                        /* due componenti diverse */
    redefineborder (b2, bl1);
    topo_change_l (b1, b2);
    assert (bl2 == extractborderlist (bl2));
    bl2->sponda = 0;
    freeborderlist (bl2);
  }
  return (1);
}

void
topo_change_l (struct border *b1, struct border *b2)
{
  struct border *btemp;

  btemp = b1->next;
  b1->next = b2->next;
  b2->next = btemp;
}

struct borderlist *
extractborderlist (struct borderlist *bl)
{
  struct borderlist *blscan;
  struct region *r;

  r = bl->region;
  if (r->border == bl)
  {
    r->border = bl->next;
    bl->next = 0;
    return (bl);
  }
  for (blscan = r->border; blscan; blscan = blscan->next)
  {
    if (blscan->next == bl)
    {
      blscan->next = bl->next;
      bl->next = 0;
      return (bl);
    }
  }
  return (0);
}

/*
 * unione di due regioni (che alla fine produrranno una
 * regione con piu' componenti connesse)
 * dopo questa operazione la situazione non e' canonica
 * relativamente alla connessione e al bordo esterno.
 * il puntatore restituito sara' quello di r1.
 */

struct region *
regionunion (struct region *r1, struct region *r2, struct sketch *sketch)
{
  struct borderlist *bl;

  if (debug > 2) printf ("entering region union\n");
  assert (r1 != r2);
  if (debug > 2)
  {
    printf ("regione r1: \n");
    for (bl = r1->border; bl; bl = bl->next)
    {
      if (bl->sponda) printborder (bl->sponda, r1);
        else printf (" ()");
    }
    printf (";\n");
    printf ("regione r2: \n");
    for (bl = r2->border; bl; bl = bl->next)
    {
      if (bl->sponda) printborder (bl->sponda, r2);
        else printf (" ()");
    }
    printf (";\n");
  }
  redefineregion (r2, r1);
  /* riunisco le componenti connesse */
  for (bl = r1->border; bl; bl = bl->next)
  {
    if (bl->next == 0)
    {
      bl->next = r2->border;
      r2->border = 0;
      break;
    }
  }
  removeregion (r2, sketch);
  if (debug > 2)
  {
    printf ("regione r1finale: \n");
    for (bl = r1->border; bl; bl = bl->next)
    {
      if (bl->sponda) printborder (bl->sponda, r1);
        else printf (" ()");
    }
    printf (";\n");
  }
  // if (debug) printf ("exiting region union\n");
  return (r1);
}

/*
 * ridefinisco tutti i puntatori della regione r1 in modo che
 * puntino alla regione r2
 */

void 
redefineregion (struct region *r1, struct region *r2)
{
  struct borderlist *bl;

  if (r1 == r2) return;   /* nulla da fare in questo caso */

  for (bl = r1->border; bl; bl = bl->next) bl->region = r2;
}

/*
 * ridefinisco tutti i puntatori del bordo b1 in modo che
 * puntino al bordo bl2
 */

void 
redefineborder (struct border *b1, struct borderlist *bl2)
{
  struct border *bp;
  struct borderlist *bl1;

  bl1 = b1->border;
  if (bl1 == bl2) return;   /* nulla da fare in questo caso */

  bp = b1; 
  do {
    assert (bp->border == bl1);
    bp->border = bl2;
    bp = bp->next;
  } while (bp != b1);
}

/*
 * bl1 e bl2 sono c.c. del bordo di due regioni r1 e r2 distinte;
 * "transfer_isles" trasferisce intere "isole" da una regione all'altra
 * in un senso o nell'altro in base al valore di "hei"; ogni bit di
 * "hei" uguale a 1 indica la volonta' di trasferire una isola.
 *
 * Attenzione: questa descrizione e' intuitiva se bl1 e bl2 sono il bordo
 * esterno delle rispettive regioni, ma l'operazione viene fatta anche
 * se (ad esempio) bl1 e' il bordo di un buco di r1, nel qual caso
 * il bordo esterno viene trattato come se fosse un'isola, mentre bl1
 * viene trattato come fosse il bordo esterno. In questo caso
 * l'esito finale non rispettera' la convenzione sulla precedenza del
 * bordo esterno nella lista delle c.c. del bordo delle regioni.
 *
 * Questa routine *non* effettua controlli di consistenza, che sono
 * quindi responsabilita' del chiamante.
 */

void
transfer_isles (struct borderlist *bl1, struct borderlist *bl2, int hei)
{
  struct region *r1, *r2;
  struct borderlist *holes, *hpt, *hptn;
  int countholes = 0;

  r1 = bl1->region;
  r2 = bl2->region;

  if (r1 == r2) return;   /* stessa regione, non faccio nulla */

  bl1 = extractborderlist (bl1);
  bl2 = extractborderlist (bl2);

  holes = r1->border;
  if (holes == 0)
  {
    holes = r2->border;
  } else {
    hpt = holes;
    while (hpt->next) hpt = hpt->next;
    hpt->next = r2->border;
  }
  r1->border = bl1;
  r2->border = bl2;
  hpt = holes;
  while (hpt)
  {
    countholes++;
    hptn = hpt->next;
    if (hei & 1)
    {
      if (hpt->region == r1) hpt->region = r2;
        else hpt->region = r1;
    }
    hei >>= 1;
    if (hpt->region == r1)
    {
      hpt->next = bl1->next;
      bl1->next = hpt;
    } else {
      hpt->next = bl2->next;
      bl2->next = hpt;
    }
    hpt = hptn;
  }
  if (hei) fprintf (stderr, 
     "more changes requested than present islands (%d)\n", countholes);
}

/*
 * controlla la consistenza interna di tutti i puntatori
 * di uno sketch
 */

int
checkconsistency (struct sketch *sketch)
{
  struct region *r;
  struct borderlist *bl;
  struct border *bp;
  struct arc *arc;

  for (r = sketch->regions; r; r = r->next)
  {
    assert (r);
    for (bl = r->border; bl; bl = bl->next)
    {
      assert (bl);
      if (bl->region != r)
        printf ("bl->region: %d, r: %d\n", bl->region->tag, r->tag);
      assert (bl->region == r);
      if (bl->sponda)
      {
        bp = bl->sponda;
        do {
          assert (bp->border == bl);
          arc = bp->info;
          if (arc == 0)
          {
            printf ("Warning: missing pointer from sponda to arc\n");
          }
          if (arc && arc->tag != INFINITY_ARC && arc->regionleft && arc->regionright)
          {
            assert (arc->regionleft == bp || arc->regionright == bp);
            assert (arc->regionleft != arc->regionright);
            assert (arc->regionleft->orientation != arc->regionright->orientation);
          }
          bp = bp->next;
        } while (bp != bl->sponda);
      }
    }
  }
  return (1);
}

struct sketch *
newsketch ()
{
  struct sketch *s;

  s = (struct sketch *) malloc (sizeof (struct sketch));
  s->arcs = 0;
  s->regions = 0;
  s->extregion = 0;   /* if zero: external region is the first in the list */
  s->arccount = 0;
  s->regioncount = 0;
  s->cc_tagged = 0;
  s->huffman_labelling = 0;
  s->isempty = 0;  /* we assume that empty is an exceptional situation */
  return (s);
}

struct region *
newregion (struct sketch *s)
{
  struct region *r, *er;

  r = (struct region *) malloc (sizeof (struct region));
  r->next = 0;
  r->f = 0;
  r->strati = 0;
  r->border = 0;
  if ((er = s->regions))
  {
    r->next = er->next;
    er->next = r;
  } else {
    s->regions = r;
  }
  r->tag = s->regioncount++;  /* perche' sono numerate da zero */
  return (r);
}

struct borderlist *
newborderlist (struct region *region)
{
  struct borderlist *bl;

  bl = (struct borderlist *) malloc (sizeof (struct borderlist));
  bl->region = region;
  bl->sponda = 0;

  if (region->border)   /* c'e' gia il bordo esterno */
  {
    bl->next = region->border->next;
    region->border->next = bl;
  } else {
    bl->next = region->border;  /* (=0) e' il bordo esterno */
    region->border = bl;
  }
  return (bl);
}

/* a differenza di newborderlist, il nuovo bordo e' aggiunto
 * in coda alla lista
 */

struct borderlist *
newborderlist_tail (struct region *region)
{
  struct borderlist *bl, *bl2;

  bl = (struct borderlist *) malloc (sizeof (struct borderlist));
  bl->region = region;
  bl->sponda = 0;
  bl->next = 0;
  if (region->border)
  {
    for (bl2 = region->border; bl2->next; bl2 = bl2->next);
    assert (bl2->next == 0);
    bl2->next = bl;
  } else {
    region->border = bl;
  }
  return (bl);
}

struct border *
newborder (struct borderlist *bl)
{
  struct border *b;

  b = (struct border *) malloc (sizeof (struct border));
  b->next = b;
  b->border = bl;
  return (b);
}

struct arc *
newarc (struct sketch *s)
{
  struct arc *a;

  a = (struct arc *) malloc (sizeof (struct arc));
  a->next = s->arcs;
  s->arcs = a;
  a->depths = 0;
  a->cusps = a->depthsdim = a->dvalues = -1;
  a->endpoints = -1;
  a->tag = ++s->arccount;  /* sono numerati a partire da 1 */
  a->transparent = 0;
  a->regionleft = a->regionright = 0;
  a->link_number = F_UNDEF;
  return (a);
}

/*
 * elimina SOLO un tratto di bordo
 * restituisce il tratto precedente, che
 * ora definisce un circuito senza "b"
 */

struct border *
removeborder (struct border *b)
{
  struct border *bp;

  // if (debug) printf ("removeborder...\n");
  assert (b->next != b);

  ensurecanremoveborder (b);

  bp = prevborder (b);

  bp->next = b->next;
  assert (b->info == 0);
  free (b);
  return (bp);
}

void
ensurecanremoveborder (struct border *b)
{
  assert (b->next != b);
  //region = b->region;
  if (findinborder (b->border->sponda, b))
  {
    if (b->border->sponda == b) b->border->sponda = b->next;
    return;
  }
  assert (0);
}

void
removearc (struct arc *arc, struct sketch *sketch)
{
  int tag;
  struct arc *a, *prevarc;

  assert (arc->regionleft == 0);
  assert (arc->regionright == 0);
  if (arc->depths) free (arc->depths);

  prevarc = 0;
  tag = arc->tag;
  for (a = sketch->arcs; a; a = a->next)
  {
    if (a->next == arc) prevarc = a;
    if (a->tag > tag) a->tag--;
  }
  if (sketch->arcs == arc)
  {
    sketch->arcs = arc->next;
    free (arc);
    return;
  }
  if (prevarc)
  {
    prevarc->next = arc->next;
    free (arc);
    return;
  }
  fprintf (stderr, "Error: cannot find arc to remove\n");
}

/*
 * remove a region from the list of regions.
 * the region border list must already be dealt
 * with (and the pointer must already be zero)
 * if allocated, the "strati" vector is also
 * freed.
 * In the end the data structure of the region is
 * freed
 */

void
removeregion (struct region *region, struct sketch *sketch)
{
  int tag;
  struct region *r, *prevregion;

  assert (region->border == 0);

  if (region->strati) free (region->strati);
  prevregion = 0;
  tag = region->tag;
  for (r = sketch->regions; r; r = r->next)
  {
    if (r->next == region) prevregion = r;
    if (r->tag > tag) r->tag--;
  }
  if (sketch->regions == region)
  {
    sketch->regions = region->next;
    free (region);
    return;
  }
  if (prevregion)
  {
    prevregion->next = region->next;
    free (region);
    return;
  }
  fprintf (stderr, "Error: cannot find region to remove\n");
}

void
freesketch (struct sketch *sketch)
{
  if (sketch->regions) freeregion (sketch->regions);
  if (sketch->arcs) freearc (sketch->arcs);
  free (sketch);
}

void
freeregion (struct region *region)
{
  if (region->border) freeborderlist (region->border);
  if (region->next) freeregion (region->next);
  if (region->strati) free (region->strati);
  free (region);
}

void
freeborderlist (struct borderlist *bl)
{
  if (bl->sponda) freeborder (bl->sponda);
  if (bl->next) freeborderlist (bl->next);
  free (bl);
}

/*
 * a border circular list is removed from memory
 * The pointer from the arcs (regionleft/regionright)
 * that point into the border are set to zero.
 */

void
freeborder (struct border *border)
{
  struct border *nextb;
  /* l'idea e' di spezzare il circuito e usare la ricorsione */
  assert (border->next);
  nextb = border->next;
  border->next = 0;
  freeborderdl (nextb);
  return;
}

void
freeborderdl (struct border *border)
{
  struct arc *arcdata;

  if ((arcdata = border->info))
  {
    if (arcdata->regionleft == border) arcdata->regionleft = 0;
    if (arcdata->regionright == border) arcdata->regionright = 0;
  }

  if (border->next) freeborderdl (border->next);
  free (border);
}

void
freearc (struct arc *arc)
{
  if (arc->regionleft || arc->regionright)
  {
    printf ("Warning: clearing arc data with nonempty region pointer\n");
  }
  if (arc->depths) free (arc->depths);
  if (arc->next) freearc (arc->next);
  free (arc);
}
