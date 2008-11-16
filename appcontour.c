#include <assert.h>
#include "contour.h"

extern int debug;
extern int quiet;

/* prototipi */
int checkdnodecons (struct border *b, int dd[4]);
int getdatnode (struct border *b);
struct border *reverse_border (struct border *);
int add_s1 (struct sketch *s, struct region *r, 
            int stratum, int ori);
/* fine prototipi */

/*
 * find region data from tag
 */

struct region *
findregion (struct sketch *sketch, int tag)
{
  struct region *r;

  for (r = sketch->regions; r; r = r->next)
  {
    if (r->tag == tag) return (r);
  }
  return (0);
}

/*
 * find arc data from tag
 */

struct arc *
findarc (struct sketch *sketch, int tag)
{
  struct arc *a;

  for (a = sketch->arcs; a; a = a->next)
  {
    if (a->tag == tag) return (a);
  }
  return (0);
}

/*
 * show various information about the apparent contour
 */

void
showinfo (struct sketch *sketch)
{
  int numarcs, numsmallarcs, nums1s, numcusps, numregions, numcrossings;
  int numarcsend1, numholes, twiceohmotoinvariant, poscusps;
  int numlcomponents;
  int i, d, dmin;
  double ohmotoinvariant;
  struct arc *arc;
  struct region *r;
  struct borderlist *bl;

  numarcs = numsmallarcs = numcusps = nums1s = numarcsend1 = 0;
  poscusps = 0;
  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    numarcs++;
    if (arc->endpoints == 0) nums1s++;
    if (arc->endpoints == 1) numarcsend1++;
    numcusps += arc->cusps;
    numsmallarcs += arc->cusps + 1;
    if (arc->endpoints == 0 && arc->cusps > 0) numsmallarcs--;
    if (sketch->huffman_labelling && arc->cusps > 0)
    {
      for (i = 0; i < arc->cusps; i++)
      {
        dmin = arc->depths[i];
        if ((d = arc->depths[i+1]) < dmin) dmin = d;
        if ((dmin % 2) == 0) poscusps++;
      }
    }
  }

  numcrossings = (numarcs - nums1s)/2;

  numregions = numholes = 0;
  for (r = sketch->regions; r; r = r->next)
  {
    numregions++;
    for (bl = r->border->next; bl; bl = bl->next) numholes++;
  }

  twiceohmotoinvariant = compute_ohmoto (sketch);
  ohmotoinvariant = twiceohmotoinvariant/2.0;

  numlcomponents = count_link_components (sketch);

  if (! quiet)
    printf ("This is an apparent contour %s Huffman labelling\n", 
      sketch->huffman_labelling?("with"):("without"));
    else printf ("Huffman labelling:  %d\n", 
      sketch->huffman_labelling?1:0);

  if (! quiet) printf ("\nProperties of the 3D surface:\n");
  if (sketch->huffman_labelling)
    printf ("Connected comp.:    %d\n", count_connected_components (sketch));
  printf ("Total Euler ch.:    %d\n", euler_characteristic (sketch));

  if (! quiet) printf ("\nAicardi-Ohmoto invariants:\n");
  printf ("Cusps:              %d\n", numcusps);
  if (sketch->huffman_labelling)
    printf ("Positive cusps:     %d\n", poscusps);
  printf ("Crossings:          %d\n", numcrossings);
  printf ("Bennequin:          %.1lf\n", ohmotoinvariant);

  if (! quiet) printf ("\nProperties of the 2D apparent contour:\n");
  printf ("Arcs:               %d\n", numsmallarcs);
  printf ("Extended arcs:      %d\n", numarcs);
  printf ("Link components:    %d\n", numlcomponents);
  printf ("Loops:              %d\n", nums1s);
  printf ("Nodes (cusps+cross):%d\n", numcusps + numcrossings);
  printf ("Regions:            %d\n", numregions);
  printf ("Connected comp.:    %d\n", numholes);

}

/*
 * ad an S1 to the apparent contour in the given region
 * and with the given d value.
 *
 * if ori == 1  this corresponds to adding a sphere
 * if ori == -1 this corresponds to punching a hole in two
 * strata and glueing them at the cut
 */

int
add_s1 (struct sketch *s, struct region *r, int dval, int ori)
{
  struct region *newr;
  struct arc *newa;
  struct borderlist *newbl1, *newbl2;
  struct border *newbp1, *newbp2;

  assert (s && r && dval >= 0);
  assert (ori == 1 || ori == -1);
  if (ori > 0 && dval > r->f) return (0);
  if (ori < 0 && dval > r->f - 2) return (0);

  /* devo aggiungere un s^1 senza cuspidi */
  newa = newarc (s);
  newa->depths = (int *) malloc (sizeof (int));
  newa->depthsdim = 1;
  newa->cusps = 0;
  newa->dvalues = 1;
  newa->depths[0] = dval;
  newa->endpoints = 0;

  newr = newregion (s);
  newbl1 = newborderlist (newr);
  newbp1 = newborder (newbl1);
  newbl1->sponda = newbp1;
  newbp1->next = newbp1;
  newbl2 = newborderlist (r);
  newbp2 = newborder (newbl2);
  newbl2->sponda = newbp2;
  newbp2->next = newbp2;
  newbp1->info = newbp2->info = newa;
  if (ori > 0) {
    newa->regionleft = newbp1;
    newa->regionright = newbp2;
  } else {
    newa->regionleft = newbp2;
    newa->regionright = newbp1;
  }
  newbp1->orientation = ori;
  newbp2->orientation = -ori;

  if (debug) printsketch (s);
  return (1);
}

/*
 * switch front-back (end up with a specular 3D surface)
 */

int
frontback (struct sketch *s)
{
  struct arc *arc;
  int f_low, i;

  for (arc = s->arcs; arc; arc = arc->next)
  {
    f_low = arc->regionright->border->region->f;
    assert (f_low != F_UNDEF);
    for (i = 0; i <= arc->cusps; i++)
    {
      arc->depths[i] = f_low - arc->depths[i];
    }
  }
  return (1);
}

/*
 * switch left-right (end up with a specular 3D surface)
 */

int
leftright (struct sketch *s)
{
  struct arc *arc;
  struct region *r;
  struct borderlist *bl;
  int ndvmezzi, i, dtemp;

  for (arc = s->arcs; arc; arc = arc->next)
  {
    ndvmezzi = (arc->cusps + 1)/2;
    for (i = 0; i < ndvmezzi; i++)
    {
      dtemp = arc->depths[i];
      arc->depths[i] = arc->depths[arc->cusps - i];
      arc->depths[arc->cusps - i] = dtemp;
    }
  }
  for (r = s->regions; r; r = r->next)
  {
    for (bl = r->border; bl; bl = bl->next)
    {
      if (bl->sponda == 0) continue;
      bl->sponda = reverse_border (bl->sponda);
    }
  }
  return (1);
}

/*
 * recursively implemented.  It is important that the data
 * is not moved, created, removed, and only the "next" pointers
 * are changed.
 * Now it is clear that this can be done more efficiently with just a single
 * iteration, but I am too lazy...
 */

struct border *
reverse_border (struct border *b)
{
  struct border *p;
  struct border *n;

  if (b == b->next) return (b);

  p = prevborder (b);
  n = b->next;
  p->next = n;
  b->next = p;
  reverse_border (p);
  n->next = b;

  return (b);
}

/*
 * change the external region.  This corresponds to working on the
 * compactification of R^2 (S^2) and sliding the north pole into
 * the selected region.  In the end the value of f at infinity can
 * be different from zero.
 */

int
changeextregion (struct sketch *s, int tag)
{
  struct region *r, *newext = 0, *oldext = 0;
  struct borderlist *bl;
  int changes;

  for (r = s->regions; r; r = r->next)
  {
    if (r->tag == tag) newext = r;
    if (r->border->sponda == 0) oldext = r;
  }
  if (newext == oldext) return (1);
  assert (newext && oldext);
  bl = extractborderlist (oldext->border);
  bl->region = newext;
  bl->next = newext->border;
  newext->border = bl;
  changes = adjust_isexternalinfo (s);
  return (1);
}

/*
 * compute the euler characteristic
 */

int
euler_characteristic (struct sketch *sketch)
{
  struct arc *arc;
  struct region *r;
  struct borderlist *bl;
  int cusps = 0;
  int arcweight = 0;
  int layers = 0;
  int halfarcweight, factor;

  assert (sketch->regions->border->sponda == 0);
  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    cusps += arc->cusps;
    if (arc->endpoints == 0) continue;   /* questi non contano */
    arcweight += arcmult (arc);
  }

  for (r = sketch->regions; r; r = r->next)
  {
    factor = 2;
    for (bl = r->border; bl; bl = bl->next)
      if (bl->sponda) factor--;          /* regione esterna, escluso il bordo nullo */
                                         /* peso di uno strato: 1-buchi */
    layers += factor * r->f;
  }

  halfarcweight = arcweight/2;
  assert (arcweight == 2*halfarcweight);
  return (layers - cusps - halfarcweight);
}

int
appcontourcheck (struct sketch *sketch, int huffman, int notquiet)
{
  int fail, globfail, i, diff=0;
  struct region *region;
  struct arc *arc;
  int d, dmin = 0, dmax = 0, fmin = 0;
  struct border *failb = 0;
  struct borderlist *hole;
  int dd[4], res = 0;

  if (notquiet)
  {
    printf ("Checking consistency as apparent contour");
    if (huffman) printf (" with huffman labelling");
    printf ("\n");
  }

  globfail = 0;
  if (notquiet) printf ("1. Checking arc orientation across nodes... ");
  fflush (stdout);
  fail = 0;
  for (region = sketch->regions; region; region = region->next)
  {
    for (hole = region->border; hole; hole = hole->next)
    {
      if (hole->sponda && checkorientationborder (hole->sponda) == 0)
      {fail = 1; break;}
    }
    if (fail) break;
  }
  if (fail)
  {
    globfail = 1;
    if (notquiet) printf ("FAILED while checking region %d\n", region->tag);
  } else if (notquiet) printf ("OK\n");

  if (notquiet) printf ("2. Checking positivity of f...              ");
  fflush (stdout);
  fail = 0;
  for (region = sketch->regions; region; region = region->next)
  {
    if (region->f < 0 || region->f == F_UNDEF) {fail = 1; break;}
  }
  if (fail)
  {
    globfail = 1;
    if (notquiet) printf ("FAILED for region %d (f = %d)\n", region->tag, region->f);
  } else if (notquiet) printf ("OK\n");

  if (huffman)
  {
    if (notquiet) printf ("3. d versus f consistency...                ");
    fflush (stdout);

    fail = 0;
    for (arc = sketch->arcs; arc; arc = arc->next)
    {
      assert (arc->depths);
      dmin = dmax = arc->depths[0];
      for (i = 1; i < arc->dvalues; i++)
      {
        d = arc->depths[i];
        if (d < dmin) dmin = d;
        if (d > dmax) dmax = d;
      }
      fmin = arc->regionleft->border->region->f;
      if (arc->regionright->border->region->f < fmin)
        fmin = arc->regionright->border->region->f;
      if (0 <= dmin && dmax <= fmin) continue;
      fail = 1;
      break;
    }
    if (fail)
    {
      globfail = 1;
      if (notquiet) printf ("FAILED for arc %d (dmin = %d, dmax = %d, fmin = %d)\n",
        arc->tag, dmin, dmax, fmin);
    } else if (notquiet) printf ("OK\n");

    if (notquiet) printf ("4. variation of d across cusps...           ");
    fflush (stdout);

    fail = 0;
    for (arc = sketch->arcs; arc; arc = arc->next)
    {
      if (arc->cusps < 0) {fail = 1; break;}
      if (arc->cusps == 0) continue;
      for (i = 0; i < arc->cusps; i++)
      {
        diff = arc->depths[i] - arc->depths[i+1];
        if (diff == 1 || diff == -1) continue;
        fail = 1;
        break;
      }
      if (fail) break;
    }
    if (fail)
    {
      globfail = 1;
      if (notquiet) printf ("FAILED for arc %d (diff = %d)\n", arc->tag, diff);
    } else if (notquiet) printf ("OK\n");

    if (notquiet) printf ("5. values of d on nodes...                  ");
    fflush (stdout);

    fail = 0;
    for (arc = sketch->arcs; arc; arc = arc->next)
    {
      if (arc->endpoints == 0) continue;
      failb = arc->regionleft;
      if ((res = checkdnodecons (arc->regionleft, dd)) == 1)
      {
        if (arc->endpoints == 1) continue;
        failb = arc->regionright;
        if ((res = checkdnodecons (arc->regionright, dd)) == 1) continue;
      }
      fail = 1;
      break;
    }
  }
  if (fail)
  {
    globfail = 1;
    if (notquiet)
    {
      printf ("FAILED ");
      if (res == 2)
      {
        printf ("due to wrong arc orientation\n");
      } else {
        printf ("due to wrong d values\n");
      }
      printf ("          ");
      printf ("a%d (d=%d), ", failb->info->tag, dd[0]);
      failb = failb->next;
      arc = failb->info;
      if (arc->regionleft == failb) failb = arc->regionright;
        else failb = arc->regionleft;
      printf ("a%d (d=%d), ", failb->info->tag, dd[1]);
      failb = failb->next;
      arc = failb->info;
      if (arc->regionleft == failb) failb = arc->regionright;
        else failb = arc->regionleft;
      printf ("a%d (d=%d), ", failb->info->tag, dd[2]);
      failb = failb->next;
      arc = failb->info;
      if (arc->regionleft == failb) failb = arc->regionright;
        else failb = arc->regionleft;
      printf ("a%d (d=%d)\n", failb->info->tag, dd[3]);
    }
  } else if (notquiet) printf ("OK\n");


  if (globfail && notquiet) printf ("This sketch is NOT an apparent contour");
  if (globfail && notquiet && huffman) printf (" with huffman labelling");
  if (globfail && notquiet) printf ("\n");
  return (globfail == 0);
}

int
checkorientationborder (struct border *b)
{
  struct border *bp, *bpn, *bptrans, *bpopp;
  struct arc *arc;

  bp = b;
  do {
    bpn = bp->next;
    arc = bpn->info;
    if (arc->endpoints == 0) continue;
    if (arc->regionleft == bpn) bptrans = arc->regionright;
      else bptrans = arc->regionleft;
    bpopp = bptrans->next;
    if (bp->orientation != bpopp->orientation) return (0);
    bp = bp->next;
  } while (bp != b);

  return (1);
}

int
getdatnode (struct border *b)
{
  struct arc *arc;

  arc = b->info;
  if (b->orientation > 0)
    /* l'arco e' orientato allo stesso modo, buona l'ultima d */
    return (arc->depths[arc->dvalues - 1]);

  /* orientamento opposto, buona la prima d */
  return (arc->depths[0]);
}

/*
 * controllo di consistenza dei valori di d in un
 * nodo come contorno apparente
 */

int
checkdnodecons (struct border *b, int dd[4])
{
  struct arc *arc;
  int i, diff, imin;
  int ori[5];

  dd[0] = getdatnode (b);
  ori[0] = b->orientation;

  for (i = 1; i < 4; i++)
  {
    b = b->next;
    arc = b->info;
    b = gettransborder (b);

    dd[i] = getdatnode (b);
    ori[i] = b->orientation;
  }
  ori[4] = ori[0];   /* cosi' evito di dover calcolare un mod 4 */

  imin = -1;         /* l'indice in cui si trova il valore minimo che sale */
  /* ora ho i 4 valori di 'd' e le orientazioni */
  if (dd[0] == dd[2])
  {
    if (dd[1] < dd[0] || dd[3] < dd[0]) return (0);
    diff = dd[1] - dd[3];
    imin = (diff > 0)?3:1;
    if (diff != 2 && diff != -2) return (0);
  } else {
    if (dd[1] != dd[3]) return (0);
    if (dd[0] < dd[1] || dd[2] < dd[1]) return (0);
    diff = dd[0] - dd[2];
    imin = (diff > 0)?2:0;
    if (diff != 2 && diff != -2) return (0);
  }
  /* ultimo controllo: l'orientazione in imin+1 deve essere 1 */
  if (ori[imin+1] != 1) return (2);
  return (1);
}

/*
 * rimozione di una componente connessa.
 * NOTA: le componenti connesse sono numerate a partire da 0.
 */

/* prototipi locali */
void remove_transparent_arcs (struct sketch *sketch);
int make_transparent (int ccid, struct sketch *sketch);
int join_consecutive_arcs (struct sketch *sketch);
int list_strati (struct sketch *sketch);    /* for debugging purposes */

int
extract_connected_component (int ccid, struct sketch *sketch)
{
  int ccnum, i;

  ccnum = count_connected_components (sketch);
  if (ccid < 0 || ccid >= ccnum) return (1);

  for (i = 0; i < ccnum; i++)
  {
    if (i != ccid)
    {
      if (remove_connected_component (i, sketch) == 0) return (0);
    }
  }
  return (1);
}

int
remove_connected_component (int ccid, struct sketch *sketch)
{

  if (tag_connected_components (sketch) < 0) return (-1);

  if (debug) printf ("removing component %d\n", ccid);
  if (debug) list_strati (sketch);
  /* primo passo, rendo "trasparente" la superficie */
  make_transparent (ccid, sketch);
  remove_transparent_arcs (sketch);
  if (debug) printsketch (sketch);
  join_consecutive_arcs (sketch);
  adjust_isexternalinfo (sketch);
  postprocesssketch (sketch);
  if (debug) list_strati (sketch);
  if (debug) printsketch (sketch);
  return (1);
}

int
list_strati (struct sketch *sketch)
{
  struct region *r;
  int i;

  for (r = sketch->regions; r; r = r->next)
  {
    printf ("region %d: ", r->tag);
    for (i = 0; i < r->f; i++)
    {
      printf ("%d ", r->strati[i]);
    }
    printf ("\n");
  }
  return (1);
}

int
make_transparent (int ccid, struct sketch *sketch)
{
  struct arc *arc;
  struct region *rs, *r;
  int d, i, j, dmin, nstrati;

  /* etichetto gli archi coinvolti come trasparenti
   * e aggiorno il valore di d per gli altri
   */
  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    rs = arc->regionleft->border->region;
    assert (rs->f > arc->regionright->border->region->f);
    d = arc->depths[0];
    if (debug && rs->strati[d] == ccid) printf ("arco %d trasparente\n", arc->tag);
    if (rs->strati[d] == ccid) arc->transparent = 1;
    /* calcolo il numero di strati trasparenti davanti all'arco */
    dmin = BIG_INT;
    for (i = 0; i < arc->dvalues; i++)
    {
      if (arc->depths[i] < dmin) dmin = arc->depths[i];
    }
    nstrati = 0;
    for (i = 0; i < dmin; i++)
    {
      if (rs->strati[i] == ccid) nstrati++;
    }
    for (i = 0; i < arc->cusps + 1; i++) arc->depths[i] -= nstrati;
  }

  /* aggiorno la f e gli strati delle regioni */
  for (r = sketch->regions; r; r = r->next)
  {
    for (i = j = 0; i < r->f; i++)
    {
      if (r->strati[i] != ccid) r->strati[j++] = r->strati[i];
    }
    if (debug && i != j) printf ("n. strati cambia da %d a %d\n", i, j);
    r->f = j;
  }
  return (1);
}

void
remove_transparent_arcs (struct sketch *sketch)
{
  struct arc *arc;
  struct border *bl, *br, *blp, *brp, *bltemp;
  struct borderlist *bll, *blr;
  struct region *rl, *rr;

  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    if (arc->transparent == 0) continue;
    bl = arc->regionleft;
    bll = bl->border;
    rl = bll->region;
    br = arc->regionright;
    blr = br->border;
    rr = blr->region;
    bl->info = br->info = 0;  /* rimuovi il riferimento all'arco */
    assert (rl->f == rr->f);
    if (rl != rr) rl = rr = regionunion (rl, rr, sketch);
    /*
     * resolving bug of 2008.11.13:
     * it might happen that we remain with a node with a loop
     * and no other arcs coming out (e.g. a swallowtail)
     * in this case we simply treat it as an isolated S1
     */
    if (arc->endpoints == 1 && bl->next == bl && br->next == br)
      arc->endpoints = 0;
    switch (arc->endpoints)
    {
      case 0:
      /* e' un S1 isolato */
      bll = extractborderlist (bll);
      freeborderlist (bll);
      blr = extractborderlist (blr);
      freeborderlist (blr);
      break;

      case 1:
      if (bl->next != bl)
      {
        bltemp = bl;
        bl = br;
        br = bltemp;
        bll = bl->border;
      }
      assert (bl->next == bl);
      assert (br->next != br);
      assert (bl->info == 0 && br->info == 0);
      bll = extractborderlist (bll);
      freeborderlist (bll);
      brp = removeborder (br);
      break;

      case 2:
      assert (bl->next != bl);
      assert (br->next != br);
      assert (bl->info == 0 && br->info == 0);
      blp = removeborder (bl);
      if (br->next == br)
      {   /* caso particolare di rimozione di un baffo */
        blr = extractborderlist (blr);
        freeborderlist (blr);
      } else {
        brp = removeborder (br);
        if (blp != brp && blp != br) topo_change (blp, brp);
      }
      break;
    }
    arc->regionleft = arc->regionright = 0;
    removearc (arc, sketch);
  }
  if (debug) printsketch (sketch);
  return;
}

int
join_consecutive_arcs (struct sketch *sketch)
{
  int goon;
  struct arc *arc, *arc1, *arc2;
  struct border *b1, *b2;

  goon = 1;
  while (goon)
  {
    goon = 0;          /* set goon to 1 if whole procedure to be repeated */
    for (arc1 = sketch->arcs; arc1; arc1 = arc1->next)
    {
      /* this arc has to be merged with the following if
       * the next node has exactly two concurring arcs
       * this can be checked as follows:
       */
      if (arc1->endpoints == 0) continue;
      b1 = arc1->regionleft;
      b2 = gettransborder (b1->next);
      if (b2->next->info == arc1) /* it is interesting to observe that this is
                                   * a necessary and sufficient condition
                                   */
      {
        assert (b1->orientation == 1);
        assert (b2->orientation == -1);
        arc2 = b2->info;
        assert (arc2->regionleft->info == arc2);
        if (debug) printf ("calling mergearcs: %d, %d\n", arc1->tag, arc2->tag);
        arc = mergearcs (arc1, arc2, sketch);
        if (b1->next != b1) removeborder (b1->next);
        if (b2->next != b2) removeborder (b2->next);
        arc->regionleft = b1;
        arc->regionright = b2;
        b1->info = arc;
        b2->info = arc;
        if (debug) printsketch (sketch);
        goon = 1;
        if (arc != arc1) {goon = 1; arc1 = arc;}
      }
    }
  }
  return (1);
}

/*
 * funzione per il calcolo del numero di componenti connesse
 * ed inizializzazione del manyfold topologico
 */

/* prototipi */
void countcc_flood (int ccid, struct sketch *sketch);
int countcc_flood_step (struct region *r, int i);

int
count_connected_components (struct sketch *sketch)
{
  struct region *r;
  int i, ccidmax = -1, ccidmin = BIG_INT;

  if (tag_connected_components (sketch) < 0) return (-1);

  for (r = sketch->regions; r; r = r->next)
  {
    for (i = 0; i < r->f; i++)
    {
      if (r->strati[i] > ccidmax) ccidmax = r->strati[i];
      if (r->strati[i] < ccidmin) ccidmin = r->strati[i];
    }
  }
  // free_connected_components (sketch);
  return (ccidmax - ccidmin + 1);
}

int
tag_connected_components (struct sketch *sketch)
{
  struct region *r;
  int i, ccid = 0, found;
  int fneg = 0;

  if (debug) printf ("count_cc, fase 1: inizializzo i vettori\n");

  if (sketch->cc_tagged) return (0);

  for (r = sketch->regions; r; r = r->next)
  {
    if (r->f < 0) fneg++;
    if (r->f > 0)
    {
      r->strati = (int *) malloc (r->f * sizeof(int));
      for (i = 0; i < r->f; i++) r->strati[i] = -1;
    }
  }

  while (1)
  {
    found = 0;
    for (r = sketch->regions; r; r = r->next)
    {
      for (i = 0; i < r->f; i++)
      {
        if (r->strati[i] == -1)
        {
          if (debug) printf ("inizializzo, ccid %d, r %d, i %d\n",
                       ccid, r->tag, i);
          r->strati[i] = ccid;
          found = 1;
          break;
        }
      }
      if (found) break;
    }
    if (found == 0) break;

    /* ho inizializzato una nuova componente connessa */
    countcc_flood (ccid++, sketch);
  }
  sketch->cc_tagged = 1;
  if (fneg) return (-1);
  return (1);
}

int
free_connected_components (struct sketch *sketch)
{
  struct region *r;

  if (debug) printf ("libero lo spazio allocato\n");
  if (sketch->cc_tagged == 0) return (0);
  for (r = sketch->regions; r; r = r->next)
  {
    if (r->f > 0)
    {
      free (r->strati);
      r->strati = 0;
    }
  }
  sketch->cc_tagged = 0;
  return (1);
}

void
countcc_flood (int ccid, struct sketch *sketch)
{
  int i, goon = 1;
  struct region *r;

  while (goon)
  {
    goon = 0;
    for (r = sketch->regions; r; r = r->next)
    {
      for (i = 0; i < r->f; i++)
      {
        if (r->strati[i] != ccid) continue;
        goon |= countcc_flood_step (r, i);
      }
    }
  }
}

int
countcc_flood_step (struct region *r, int i)
{
  int j, d, ii, goon = 0;
  struct borderlist *bl;
  struct border *bp, *btrans;
  struct region *rtrans, *rr;
  struct arc *arc;

  for (bl = r->border; bl; bl = bl->next)
  {
    bp = bl->sponda;
    if (bp == 0) continue;     /* external region.  This is possible if --fi used */
    do {
      arc = bp->info;
      btrans = gettransborder (bp);
      rtrans = btrans->border->region;
      for (j = 0; j < arc->dvalues; j++)
      {
        d = arc->depths[j];
        if (r->f > rtrans->f)
        {
          if (i == d || i == d + 1)
          {    /* rimango nella stessa regione, altro strato */
            ii = d + 1;
            if (i == d + 1) ii = d;
            rr = r;
          } else {
            ii = i;
            if (i > d + 1) ii = i - 2;
            rr = rtrans;
          }
        } else {
          /* verso regione con piu' strati... */
          rr = rtrans;
          ii = i;
          if (i >= d) ii = i + 2;
        }
        if (rr->strati[ii] != r->strati[i])
        {
          goon = 1;
          assert (rr->strati[ii] = -1);
          rr->strati[ii] = r->strati[i];
        }
      }
      bp = bp->next;
    } while (bp != bl->sponda);
  }
  return (goon);
}

int
compute_ohmoto (struct sketch *sketch)
{
  struct arc *a;
  int invariant;

  compute_link_num_arcs (sketch);
  invariant = morse_ohmoto (sketch);
  invariant *= 2;
  for (a = sketch->arcs; a; a = a->next)
  {
    if (a->cusps > 0)
    {
      //fprintf (stderr, "arc %d, cusps %d\n", a->tag, a->cusps);
      invariant += a->cusps;
      invariant -= 2*a->cusps*a->link_number;
    }
  }
  return (invariant);
}

void
compute_link_num_arcs (struct sketch *sketch)
{
  struct arc *a;
  struct region *rl, *rr;
  int finf;

  assert (sketch->regions->border->sponda == 0);
  finf = sketch->regions->f;
  if (finf) fprintf (stderr, "Warning: f is not zero at infinity (%d)\n", finf);
  for (a = sketch->arcs; a; a = a->next)
  {
    rl = a->regionleft->border->region;
    rr = a->regionright->border->region;
    a->link_number = rl->f + rr->f - 2*finf;
    a->link_number /= 2;
    //fprintf (stderr, "arc %d, link number %d\n", a->tag, a->link_number);
  }
}

int
count_link_components (struct sketch *sketch)
{
  int *arcflags;
  struct arc *a, *an;
  struct border *b;
  int i, tag, arcnum;
  int count = 0;

  arcnum = sketch->arccount;
  arcnum++;
  arcflags = (int *) malloc (arcnum * sizeof (int));

  for (i = 0; i < arcnum; i++) arcflags[i] = 0;

  for (a = sketch->arcs; a; a = a->next)
  {
    tag = a->tag;
    assert (tag < arcnum);
    if (arcflags[tag] != 0) continue;
    count++;
    arcflags[tag] = 1;
    if (a->endpoints == 0) continue;  /* this is a loop, no nodes */
                                      /* to cross                 */
    /* devo seguire questo arco attraverso i nodi */
    an = a;
    do {
      b = an->regionleft;
      b = b->next;
      b = gettransborder (b);
      b = b->next;
      an = b->info;
      if (an->regionleft != b)
      {
        fprintf (stderr, "nonconstant link orientation, cannot count links\n");
        return (-1);
      }
      if (an != a)
      {
        if (arcflags[an->tag] != 0)
        {
          fprintf (stderr, "arc already tagged while traversing ");
          fprintf (stderr, "a link component, should not happen\n");
        }
        arcflags[an->tag] = 2;
      }
    } while (an != a);
  }

  free (arcflags);
  return (count);
}
