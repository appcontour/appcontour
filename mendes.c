/*
 * computation of the mendes graph, see paper
 * D. Hacon, C. Mendes de Jesus, and M.C. Romero Fuster,
 * Stable maps from surfaces to the plane with prescribed
 * branching data, Topology and its Applications, 154 (2007) 166-175
 * and example/mendes.morse
 */

#include <assert.h>
#include <float.h>
#include "contour.h"
#include "mendes.h"

struct mendesgraph *
compute_mendes (struct sketch *s)
{
  struct mendesgraph *hg;
  int nummendesarcs, nummendesnodes;
  int **nodedata;
  int *arcdata;
  int *rstrata, *lstrata, dmin, dmax;
  int ia, dplus, dminus, nodeplus, nodeminus;
  struct arc *a;
  struct region *r, *regright, *regleft;
  int k, in, epc, epc2, epc1, epc0n, epc0c, epcdisks, lepc;
  int narcs, ncusps, necusps, nocusps;
  int evencontributes, oddcontributes;
  int nodes_x_4, dfirst;
  struct borderlist *bl;
  extern int debug;

  hg = (struct mendesgraph *) malloc (sizeof (struct mendesgraph));
  hg->sketch = s;
  hg->count_strata = 0;
  hg->x = 0;
  hg->y = 0;
  hg->nodesdata = nodedata = init_mendes_strata (s);
  hg->nummendesnodes = nummendesnodes = tag_mendes_strata (nodedata, s);

  mendes_node_canonify (hg);

  hg->arcsdata = arcdata = (int *) malloc ((s->arccount + 1) * sizeof(int));
  hg->nummendesarcs = nummendesarcs = tag_mendes_arcs (arcdata, s);

  hg->arcincplus = (int *) malloc (nummendesarcs*sizeof(int));
  hg->arcincminus = (int *) malloc (nummendesarcs*sizeof(int));
  hg->nodessign = (int *) malloc (nummendesnodes*sizeof(int));
  hg->nodesgenus = (int *) malloc (nummendesnodes*sizeof(int));

  /*
   * We compute the incidence information of arcs of the
   * Mendes graph (arcincplus and arcincminus);
   * specifically, arcincplus is used if the depth "d" of
   * the corresponding region is an even number, arcinfminus
   * is used otherwise.
   * We also compute the sign of nodes. Since we start with an
   * apparent contour with Huffman labelling, this gives an
   * embedding of the manifold M in R^3, hence it is orientable
   * and the Mendes graph is bipartite. 
   * The orientation of each mendes node is already defined 
   * when filling the arcincplus and arcincminus vectors.
   */

  for (ia = 0; ia < nummendesarcs; ia++)
  {
    for (a = s->arcs; a; a = a->next)
    {
      int d;

      if (arcdata[a->tag] != ia) continue;
      d = a->depths[0];
      r = a->regionleft->border->region;
      dplus = d;
      dminus = d+1;
      if (dplus % 2 != 0)
      {
        dplus = d+1;
        dminus = d;
      }
      rstrata = nodedata[r->tag];
      nodeplus = rstrata[dplus];
      nodeminus = rstrata[dminus];
      hg->arcincplus[ia] = nodeplus;
      hg->arcincminus[ia] = nodeminus;
      hg->nodessign[nodeplus] = +1;
      hg->nodessign[nodeminus] = -1;
      break;
    }
  }

  mendes_arc_canonify (hg);

  /*
   * Now we compute the Euler-Poincare' characteristic
   * of each Mendes node.
   */

  for (in = 0; in < hg->nummendesnodes; in++)
  {
    hg->nodesgenus[in] = -1;    /* means it is unknown */
    epcdisks = epc2 = epc1 = epc0n = epc0c = 0;
                                /* epcdisks = number of disks to be added:
                                 * number of connected components of
                                 * the boundary of the region associated
                                 * to the node;
                                 * epc2, epc1, epc0: cumulative 
                                 * characteristic of elements of dimension
                                 * 2, 1 and 0 (nodes+cusps) respectively
                                 */
    if (debug) printf ("Computing epc for node %d\n", in);
    for (ia = 0; ia < hg->nummendesarcs; ia++)
    {
      if (hg->arcincplus[ia] == in) epcdisks++;
      if (hg->arcincminus[ia] == in) epcdisks++;
    }
    for (r = s->regions; r; r = r->next)
    {
      rstrata = hg->nodesdata[r->tag];
      for (k = 0; k < r->f; k++)
      {
        if (rstrata[k] != in) continue;
        lepc = 2;
        for (bl = r->border; bl; bl = bl->next) lepc -= 1;
        epc2 += lepc;
        if (debug) printf ("region %d, stratum %d, epc = %d\n",
                      r->tag, k, lepc);
      }
    }

    nodes_x_4 = 0;
    for (a = s->arcs; a; a = a->next)
    {
      mendes_compute_cusps_data (a, &dmin, &dmax, &necusps);
      regright = a->regionright->border->region;
      rstrata = hg->nodesdata[regright->tag];
      regleft = a->regionleft->border->region;
      lstrata = hg->nodesdata[regleft->tag];
      /* first we deal with the trivial strata before dmin and
       * after dmax (horizontal connections)
       */
      if (a->endpoints != 0)
      {
        for (k = 0; k < dmin; k++)
        {
          if (rstrata[k] == in) epc1 -= 1;   /* arc contributes negatively */
        }
        for (k = dmax; k < regright->f; k++)
        {
          if (rstrata[k] == in) epc1 -= 1;
        }
        /*
         * contribution from crossings.  For internal arcs (horiz.
         * connection) we account for both starting and ending nodes
         * so we don't have problems with changing d values due to
         * cusps.  In this way all internal crossings are counted four
         * times
         */
        for (k = 0; k < regright->f; k++)
        {
          if (rstrata[k] == in) nodes_x_4 += 2;
        }
        dfirst = a->depths[0];
        if (lstrata[dfirst] == in || lstrata[dfirst+1] == in)
        {
          /* this stratum is a boundary arc for this mendes node */
          /* we count this 3 times in order to reach 4 (together
           * with the arc that enter the region
           */
          nodes_x_4 += 3;
        }
      }
      oddcontributes = evencontributes = 0;
      if (lstrata[dmin] == in) evencontributes = 1;
      if (lstrata[dmin+1] == in) oddcontributes = 1;
      if (evencontributes == 0 && oddcontributes == 0) continue;
      assert (evencontributes == 0 || oddcontributes == 0);
      ncusps = a->cusps;
      epc0c += ncusps;      /* cusps contribute to the number of vertices */
      nocusps = ncusps - necusps;
      narcs = 1;           /* initial contribution to number of arcs */
      if (a->endpoints == 0) narcs = 0;
      narcs += ncusps;     /* contribution from boundary arcs */
      if (ncusps > 0)      /* contribution from internal arcs */
      {
        if (a->endpoints > 0) narcs ++;
        if (evencontributes) narcs += necusps;
          else narcs += nocusps;
      }
      epc1 -= narcs;
    }

    /*
     * accounting for crossings
     */

    assert ((nodes_x_4 % 4) == 0);
    epc0n = nodes_x_4/4;
    //fprintf (stderr, "Computation of epc0n (nodes) missing\n");

    epc = epcdisks + epc2 + epc1 + epc0n + epc0c;
    if (debug) printf ("epcdisks = %d, epc2 = %d, epc1 = %d, ",
                     epcdisks, epc2, epc1);
    if (debug) printf ("epc0n = %d, epc0c = %d, epc = %d\n",
                     epc0n, epc0c, epc);
    if ((epc % 2) != 0)
    {
      fprintf (stderr, "Warning: Euler characteristic of node %d (%d) ",
                       in, epc);
      fprintf (stderr, "should be even!\n");
    }
    hg->nodesgenus[in] = 1 - epc/2;
  }
  return (hg);
}

/*
 * we want to tag each stratum according to
 * the mendes node to which they belong
 * we do this by first tagging all with a negative
 * dummy integer, then 0 is the tag of the first
 * mendes node.
 * in a cicle we look for an untagged stratum and give it
 * a *new* mendes tag, then we loop on all strata and try to
 * locally extend the mendes node by local adjacency
 */

int
tag_mendes_strata (int **data, struct sketch *s)
{
  int k, mendes_tag;
  int *rstrata;
  struct region *r;

  /* reset all tags to -1 */
  for (r = s->regions; r; r = r->next)
  {
    for (k = 0; k < r->f; k++)
    {
      rstrata = data[r->tag];
      rstrata[k] = -1;
    }
  }

  mendes_tag = 0;
  while (single_tag_mendes_strata(mendes_tag, data, s)) mendes_tag++;

  return (mendes_tag);
}

int
single_tag_mendes_strata (int tag, int **data, struct sketch *s)
{
  struct region *r;
  int i, k, c, count, found;
  int *rstrata;

  found = 0;
  for (r = s->regions; r; r = r->next)
  {
    i = r->tag;   /* tags start from zero */
    rstrata = data[i];
    for (k = 0; k < r->f; k++)
    {
      if (rstrata[k] < 0)
      {
        found = 1;
        break;
      }
    }
    if (found) break;
  }

  if (found == 0) return (0);  /* non ho trovato strati non etichettati */
  /* trovata una regione (r) con una strato (k) non etichettata */
  rstrata = data[r->tag];
  rstrata[k] = tag;

  count = 1;
  while ((c = mendes_try_expand_node (tag, data, s)))
  {
    count += c;
  }
  return (count);
}

int
mendes_try_expand_node (int tag, int **data, 
  struct sketch *s)
{
  int count = 0;
  int k;
  struct region *r;
  struct border *bp, *bpstart;
  struct borderlist *bl;
  int *rstrata;

  for (r = s->regions; r; r = r->next)
  {
    rstrata = data[r->tag];
    for (k = 0; k < r->f; k++)
    {
      if (rstrata[k] != tag) continue;
      for (bl = r->border; bl; bl = bl->next)
      {
        bpstart = bl->sponda;
        if (bpstart == 0) continue;
        bp = bpstart;
        do {
          count += local_mendes_try_expand_node (data, bp, k);
        } while (bp = bp->next, bp != bpstart);
      }
    }
  }
  return (count);
}

int
local_mendes_try_expand_node (int **data, 
  struct border *bp, int k)
{
  int count, d, i, dmin, dmax, ori;
  int htag, rtag, stag;
  struct arc *a;
  struct border *btrans;
  int *rdata, *sdata, *r1data, *s1data;

  a = bp->info;
  btrans = gettransborder (bp);
  stag = btrans->border->region->tag;
  rtag = bp->border->region->tag;
  assert (stag != rtag);
  rdata = data[rtag];
  sdata = data[stag];
  htag = rdata[k];
  assert (htag >= 0);
  ori = bp->orientation;
  /* find dmin and dmax */
  dmin = dmax =  a->depths[0];
  for (i = 1; i < a->dvalues; i++)
  {
    if (a->depths[i] < dmin) dmin = a->depths[i];
    if (a->depths[i] > dmax) dmax = a->depths[i];
  }

  if (k < dmin)
  {
    if (sdata[k] < 0)
    {
      sdata[k] = htag;
      return (1);
    }
    assert (sdata[k] == htag);
    return (0);
  }

  if (ori > 0 && k >= dmax + 2)
  {
    if (sdata[k-2] < 0)
    {
      sdata[k-2] = htag;
      return (1);
    }
    assert (sdata[k-2] == htag);
    return (0);
  }

  if (ori < 0 && k >= dmax)
  {
    if (sdata[k+2] < 0)
    {
      sdata[k+2] = htag;
      return (1);
    }
    assert (sdata[k+2] == htag);
    return (0);
  }

  if ((dmin % 2) != (k % 2)) dmin++;
  r1data = rdata;
  s1data = sdata;
  if (ori < 0) {r1data = sdata; s1data = rdata;}
  count = 0;
  for (d = dmin; d < dmax; d += 2)
  {
    if (s1data[d] < 0)
    {
      s1data[d] = htag;
      count++;
    } else assert (s1data[d] == htag);
  }
  for (d = dmin; d < dmax + 2; d += 2)
  {
    if (r1data[d] < 0)
    {
      r1data[d] = htag;
      count++;
    } else assert (r1data[d] == htag);
  }

  return (count);
}

int **
init_mendes_strata (struct sketch *s)
{
  int rnum, i, rtag;
  struct region *r;
  int **data;

  rnum = s->regioncount;
  data = (int **) 
         malloc (rnum*sizeof(int *));

  for (i = 0; i < rnum; i++) data[i] = 0;
  for (r = s->regions; r; r = r->next)
  {
    rtag = r->tag;
    if (r->f > 0) 
      data[rtag] = (int *)
            malloc (r->f*sizeof(int));
  }
  return (data);
}

void
mendes_compute_cusps_data (struct arc *a, 
  int *dmin, int *dmax, int *necusps)
{
  int k, ne, no, dsum, minval;

  *dmin = *dmax = a->depths[0];
  *necusps = 0;
  if (a->cusps == 0) return;
  ne = no = 0;
  for (k = 1; k <= a->cusps; k++)
  {
    if (a->depths[k] < *dmin) *dmin = a->depths[k];
    if (a->depths[k] > *dmax) *dmax = a->depths[k];
    dsum = a->depths[k] + a->depths[k-1];
    minval = (dsum - 1)/2;
    assert (2*minval + 1 == dsum);
    if ((minval % 2) == 0) ne++;
      else no++;
  }
  if ((*dmin % 2) == 0) *necusps = ne;
    else *necusps = no;

  return;
}

void
describe_mendes_nodes (int nummendesnodes, int **data, 
  struct sketch *s)
{
  int htag, k;
  struct region *r;
  int *rstrata;

  for (htag = 0; htag < nummendesnodes; htag++)
  {
    printf ("Node %d:", htag);
    for (r = s->regions; r; r = r->next)
    {
      rstrata = data[r->tag];
      for (k = 0; k < r->f; k++)
      {
        if (rstrata[k] == htag)
        {
          printf (" (%d,%d)", r->tag, k);
        }
      }
    }
    printf (";\n");
  }
}

/*
 * this data structure is an integer vector dimensioned
 * as the number of extended arcs, and will contain the
 * corresponding tag as mendes arc (these are components
 * of the apparent contour).
 */

void
describe_mendes_arcs (int nummendesarcs, int *arcdata, struct sketch *s)
{
  int htag;
  struct arc *a;

  for (htag = 0; htag < nummendesarcs; htag++)
  {
    printf ("Arc %d:", htag);
    for (a = s->arcs; a; a = a->next)
    {
      if (arcdata[a->tag] == htag)
        printf (" %d", a->tag);
    }
    printf (";\n");
  }
}

int
tag_mendes_arcs (int *arcdata, struct sketch *s)
{
  int k, mendes_tag;

  /* reset all tags to -1 */
  for (k = 0; k <= s->arccount; k++) arcdata[k] = -1;

  mendes_tag = 0;
  while (single_tag_mendes_arc(mendes_tag, arcdata, s)) mendes_tag++;

  return (mendes_tag);
}

int
single_tag_mendes_arc (int tag, int *arcdata, struct sketch *s)
{
  struct arc *a;
  int i, c, count, found;

  found = 0;
  for (a = s->arcs; a; a = a->next)
  {
    i = a->tag;   /* tags start from one */
    if (arcdata[i] < 0)
    {
      found = 1;
      break;
    }
  }

  if (found == 0) return (0);  /* non ho trovato archi non etichettati */
  /* trovato un arco (a) non etichettato */
  arcdata[a->tag] = tag;

  count = 1;
  while ((c = mendes_try_expand_arc (tag, arcdata, s))) count += c;
  return (count);
}

int
mendes_try_expand_arc (int tag, int *arcdata, struct sketch *s)
{
  int count = 0;
  struct arc *a, *an;
  struct border *b;

  for (a = s->arcs; a; a = a->next)
  {
    if (arcdata[a->tag] != tag) continue;
    if (a->endpoints == 0) continue;

    an = a;
    do {
      b = an->regionleft;
      b = b->next;
      b = gettransborder (b);
      b = b->next;
      an = b->info;
      assert (an->regionleft == b);
      if (arcdata[an->tag] < 0)
      {
        count++;
        arcdata[an->tag] = tag;
      } else assert (arcdata[an->tag] == tag);
    } while (an != a);
  }

  return (count);
}

/*
 * display mendes graph
 */

void
print_mendes (struct mendesgraph *h)
{
  int ia, in, kigobj, kigobjfa, kigobjflab;
  extern int mendesge, quiet, verbose;
  char *kigpointstyle[] = {"Round", "RoundEmpty" };
  int ptstyleid;

  switch (mendesge)
  {
    case HGE_TEXT:
    if (quiet == 0)
    {
      printf ("#\n# Mendes graph definition:\n");
      printf ("# Mendes graph with %d nodes and %d arcs.\n",
            h->nummendesnodes, h->nummendesarcs);
      printf ("#\n");
    }
    printf ("mendes {\n");
    if (quiet == 0)
    {
      printf ("# Mendes nodes information follows:\n");
      printf ("# Node id: sign genus;\n");
    }
    for (in = 0; in < h->nummendesnodes; in++)
    {
      printf ("Node %d:%c %d;\n",
        in, (h->nodessign[in] < 0)?'-':'+', h->nodesgenus[in]);
    }
    if (quiet == 0)
    {
      printf ("#\n# Mendes arcs information follows:\n");
      printf ("# Arc id: pos.-node-id neg.-node-id;\n");
    }
    for (ia = 0; ia < h->nummendesarcs; ia++)
    {
      printf ("Arc %d: %d %d;\n",
        ia, h->arcincplus[ia], h->arcincminus[ia]);
    }
    printf ("}\n");
    if (verbose)
    {
      if (quiet == 0)
      {
        printf ("\n#\n");
        printf ("# Relation between Mendes graph and Apparent contour\n");
        printf ("#\n");
      }
      printf ("mendes_contour {\n");
      if (quiet == 0)
        printf ("# Node id: (region,stratum) ... ;\n");
      describe_mendes_nodes (h->nummendesnodes, h->nodesdata, h->sketch);
      if (quiet == 0)
        printf ("# Arc id: extended_arc ...;\n");
      describe_mendes_arcs (h->nummendesarcs, h->arcsdata, h->sketch);
      printf ("}\n");
    }
    break;

    case HGE_KIG:
    mendes_xy_alloc (h);
    mendes_xy_compute (h);
//fprintf (stderr, "graph energy = %lf\n", mendes_energy (h));
    printf ("<!DOCTYPE KigDocument>\n");
    printf ("<KigDocument axes=\"1\" grid=\"1\" CompatibilityVersion=\"0.7.0\" Version=\"0.9.1\" >\n");
    printf (" <CoordinateSystem>Euclidean</CoordinateSystem>\n");
    printf (" <Hierarchy>\n");
    kigobj = 1;
    for (in = 0; in < h->nummendesnodes; in++)
    {
      printf ("  <Data type=\"double\" id=\"%d\">%lf</Data>\n",
        kigobj++, h->x[in]);
      printf ("  <Data type=\"double\" id=\"%d\">%lf</Data>\n",
        kigobj++, h->y[in]);
      printf ("  <Object type=\"FixedPoint\" id=\"%d\">\n", kigobj++);
      printf ("   <Parent id=\"%d\"/>\n", kigobj-3);
      printf ("   <Parent id=\"%d\"/>\n", kigobj-2);
      printf ("  </Object>\n");
    }
    /*
     * in this way the object id as kig object is given by
     * 3*(in + 1)
     */
    kigobjflab = kigobj;
    for (in = 0; in < h->nummendesnodes; in++)
    {
      printf ("  <Data type=\"string\" id=\"%d\">N%d</Data>\n",
        kigobj++, in);
      printf ("  <Data type=\"int\" id=\"%d\">0</Data>\n", kigobj++);
      printf ("  <Data type=\"double\" id=\"%d\">0</Data>\n", kigobj++);
      printf ("  <Data type=\"double\" id=\"%d\">0</Data>\n", kigobj++);
      printf ("  <Object type=\"RelativePoint\" id=\"%d\">\n", kigobj++);
      printf ("   <Parent id=\"%d\"/>\n", kigobj - 3);
      printf ("   <Parent id=\"%d\"/>\n", kigobj - 2);
      printf ("   <Parent id=\"%d\"/>\n", 3*(in+1));
      printf ("  </Object>\n");
      printf ("  <Object type=\"Label\" id=\"%d\">\n", kigobj++);
      printf ("   <Parent id=\"%d\"/>\n", kigobj - 5);
      printf ("   <Parent id=\"%d\"/>\n", kigobj - 2);
      printf ("   <Parent id=\"%d\"/>\n", kigobj - 6);
      printf ("  </Object>\n");
    }
    kigobjfa = kigobj;
    for (ia = 0; ia < h->nummendesarcs; ia++)
    {
      printf ("  <Object type=\"SegmentAB\" id=\"%d\">\n", kigobj++);
      in = h->arcincplus[ia];
      printf ("   <Parent id=\"%d\"/>\n", 3*(in+1));
      in = h->arcincminus[ia];
      printf ("   <Parent id=\"%d\"/>\n", 3*(in+1));
      printf ("  </Object>\n");
    }
    printf (" </Hierarchy>\n");
    printf (" <View>\n");
    for (in = 0; in < h->nummendesnodes; in++)
    {
      ptstyleid = 0;
      if (h->nodessign[in] < 0) ptstyleid = 1;
      printf ("  <Draw width=\"-1\" point-style=\"%s\" namecalcer=\"none\" style=\"SolidLine\" shown=\"true\" color=\"#0000ff\" object=\"%d\"/>\n",
        kigpointstyle[ptstyleid], kigobjflab + 6*in + 5);
      printf ("  <Draw width=\"-1\" point-style=\"%s\" namecalcer=\"%d\" style=\"SolidLine\" shown=\"true\" color=\"#0000ff\" object=\"%d\"/>\n",
        kigpointstyle[ptstyleid], 6*in + kigobjflab, 3*in + 3);
    }
    for (ia = 0; ia < h->nummendesarcs; ia++)
    {
      printf ("  <Draw width=\"-1\" point-style=\"Round\" namecalcer=\"none\" style=\"SolidLine\" shown=\"true\" color=\"#0000ff\" object=\"%d\"/>\n",
        kigobjfa + ia);
    }
    printf (" </View>\n");
    printf ("</KigDocument>\n");
    break;

    case HGE_PYKIG:
    printf ("from random import Random\n");
    printf ("g = Random()\n");
    for (in = 0; in < h->nummendesnodes; in++)
    {
      printf ("N%d = Point (g.random(), g.random(), name=\"N%d\"",
        in, in);
      if (h->nodessign[in] < 0)
        printf (", pointstyle=\"RoundEmpty\")\n");
      else
        printf (")\n");
    }
    for (ia = 0; ia < h->nummendesarcs; ia++)
    {
      printf ("a%d = Segment(N%d, N%d)\n", ia, 
        h->arcincplus[ia], h->arcincminus[ia]);
    }
    break;

    default:
    printf ("Invalid Mendes graphic engine.\n");
    exit (1);
    break;
  }  
}

/*
 * reorder all nodes based on their ordering
 */

void
mendes_node_canonify (struct mendesgraph *h)
{
  int *nodesvec, *nodesvecinv;
  int i, s, oldtag;
  int *rstrata;
  struct region *r;

  assert (h->arcsdata == 0);

  nodesvec = (int *) malloc (h->nummendesnodes*sizeof (int));
  for (i = 0; i < h->nummendesnodes; i++) nodesvec[i] = i;

  l_mendes_reorder_n (nodesvec, h->nummendesnodes, h);

  nodesvecinv = (int *) malloc (h->nummendesnodes*sizeof (int));

  for (i = 0; i < h->nummendesnodes; i++) nodesvecinv[i] = -1;
  for (i = 0; i < h->nummendesnodes; i++) nodesvecinv[nodesvec[i]] = i;
  for (i = 0; i < h->nummendesnodes; i++) assert (nodesvecinv[i] >= 0);

  /* now the actual renumbering */
  free (h->count_strata);
  h->count_strata = 0;

  for (r = h->sketch->regions; r; r = r->next)
  {
    rstrata = h->nodesdata[r->tag];
    for (s = 0; s < r->f; s++)
    {
      oldtag = rstrata[s];
      rstrata[s] = nodesvecinv[oldtag];
    }
  }
  free (nodesvec);
  free (nodesvecinv);
  return;
}

void
l_mendes_reorder_n (int *start, int num, struct mendesgraph *h)
{
  /*
   * divide and conquer reordering strategy
   */
  int numhalf, i, j, k;
  int *vcopy;

  if (num <= 1) return;

  numhalf = num/2;

  l_mendes_reorder_n (start, numhalf, h);
  l_mendes_reorder_n (start+numhalf, num - numhalf, h);

  /* now we have to merge the two ordered halves */
  vcopy = (int *) malloc (num*sizeof(int));
  for (i = 0; i < num; i++) vcopy[i] = start[i];

  for (i = 0, j = numhalf, k = 0; k < num; k++)
  {
    /* fill up start[k] by taking the minimum between vcopy[i] and vcopy[j] */
    if (i < numhalf)
    {
      if (j < num)
      {
        if (h_node_compare (vcopy[i], vcopy[j], h) <= 0)
            start[k] = vcopy[i++];
          else
            start[k] = vcopy[j++];
      } else start[k] = vcopy[i++];
    } else start[k] = vcopy[j++];
  }

  free (vcopy);
  return;
}

/*
 * the first ordering key is the number of couples (region,stratum)
 * compose each mendes node, otherwise the ordering is lessicographical
 */

int
h_node_compare (int tag1, int tag2, struct mendesgraph *h)
{
  int c1, c2;
  int *rstrata, s;
  struct region *r;

  c1 = h_count_strata (tag1, h);
  c2 = h_count_strata (tag2, h);

  if (c1 != c2) return ((c1 < c2)?(-1):1);

  for (r = h->sketch->regions; r; r = r->next)
  {
    rstrata = h->nodesdata[r->tag];
    for (s = 0; s < r->f; s++)
    {
      if (rstrata[s] == tag1 && rstrata[s] != tag2) return (-1);
      if (rstrata[s] == tag2 && rstrata[s] != tag1) return (1);
    }
  }
  return (0);
}

/*
 * count in how many couples (region,strata) a mendes node is
 * divided
 */

int
h_count_strata (int tag, struct mendesgraph *h)
{
  int *cstrata, in, s, *rstrata;
  struct region *r;

  if ((cstrata = h->count_strata) == 0)
  {
    cstrata = h->count_strata = (int *) malloc (h->nummendesnodes*sizeof (int));
    for (in = 0; in < h->nummendesnodes; in++) cstrata[in] = 0;
    for (r = h->sketch->regions; r; r = r->next)
    {
      rstrata = h->nodesdata[r->tag];
      for (s = 0; s < r->f; s++)
      {
        cstrata[rstrata[s]]++;
      }
    }
  }

  return (cstrata[tag]);
}

/*
 * reorder all arcs based on their ordering
 */

void
mendes_arc_canonify (struct mendesgraph *h)
{
  int *arcsvec, *arcsvecinv;
  int *oldarcincplus, *oldarcincminus;
  int i, oldtag;
  struct arc *a;

  arcsvec = (int *) malloc (h->nummendesarcs*sizeof (int));
  for (i = 0; i < h->nummendesarcs; i++) arcsvec[i] = i;

  l_mendes_reorder_a (arcsvec, h->nummendesarcs, h);

  arcsvecinv = (int *) malloc (h->nummendesarcs*sizeof (int));

  for (i = 0; i < h->nummendesarcs; i++) arcsvecinv[i] = -1;
  for (i = 0; i < h->nummendesarcs; i++) arcsvecinv[arcsvec[i]] = i;
  for (i = 0; i < h->nummendesarcs; i++) assert (arcsvecinv[i] >= 0);

  /* now the actual renumbering */

  oldarcincplus = h->arcincplus;
  oldarcincminus = h->arcincminus;
  h->arcincplus = (int *) malloc(h->nummendesarcs*sizeof(int));
  h->arcincminus = (int *) malloc(h->nummendesarcs*sizeof(int));
  for (i = 0; i < h->nummendesarcs; i++)
  {
    h->arcincplus[i] = oldarcincplus[arcsvec[i]];
    h->arcincminus[i] = oldarcincminus[arcsvec[i]];
  }
  free (oldarcincplus);
  free (oldarcincminus);

  for (a = h->sketch->arcs; a; a = a->next)
  {
    oldtag = h->arcsdata[a->tag];
    h->arcsdata[a->tag] = arcsvecinv[oldtag];
  }
  free (arcsvec);
  free (arcsvecinv);
  return;
}

void
l_mendes_reorder_a (int *start, int num, struct mendesgraph *h)
{
  /*
   * divide and conquer reordering strategy
   */
  int numhalf, i, j, k;
  int *vcopy;

  if (num <= 1) return;

  numhalf = num/2;

  l_mendes_reorder_a (start, numhalf, h);
  l_mendes_reorder_a (start+numhalf, num - numhalf, h);

  /* now we have to merge the two ordered halves */
  vcopy = (int *) malloc (num*sizeof(int));
  for (i = 0; i < num; i++) vcopy[i] = start[i];

  for (i = 0, j = numhalf, k = 0; k < num; k++)
  {
    /* fill up start[k] by taking the minimum between vcopy[i] and vcopy[j] */
    if (i < numhalf)
    {
      if (j < num)
      {
        if (h_arc_compare (vcopy[i], vcopy[j], h) <= 0)
            start[k] = vcopy[i++];
          else
            start[k] = vcopy[j++];
      } else start[k] = vcopy[i++];
    } else start[k] = vcopy[j++];
  }

  free (vcopy);
  return;
}

/*
 * arcs are ordered in by comparing first the positive node
 * then the negative one
 */

int
h_arc_compare (int tag1, int tag2, struct mendesgraph *h)
{
  if (h->arcincplus[tag1] != h->arcincplus[tag2])
    return ((h->arcincplus[tag1] < h->arcincplus[tag2])?(-1):1);

  if (h->arcincminus[tag1] != h->arcincminus[tag2])
    return ((h->arcincminus[tag1] < h->arcincminus[tag2])?(-1):1);

  return (0);
}

void
mendes_xy_alloc (struct mendesgraph *h)
{
  if (h->x) free (h->x);
  if (h->y) free (h->y);
  h->x = (double *) malloc (h->nummendesnodes*sizeof(double));
  h->y = (double *) malloc (h->nummendesnodes*sizeof(double));
}

/*
 * this is only a very rought starting point...
 */

void
mendes_xy_compute (struct mendesgraph *h)
{
  uint seed, optseed;
  int i;
  double e, eopt;

  optseed = 0;
  eopt = DBL_MAX;
  for (i = 0; i < 1000; i++)
  {
    seed = random ();
    srandom (seed);
    mendes_xy_randomize (h);
    e = mendes_energy (h);
    if (e < eopt)
    {
      eopt = e;
      optseed = seed;
    }
  }

  srandom (optseed);
  mendes_xy_randomize (h);
  e = mendes_energy (h);
//fprintf (stderr, "opt energy: %lf\n", e);
  return;
}

void
mendes_xy_randomize (struct mendesgraph *h)
{
  int i;

  assert (h->x && h->y);

  for (i = 0; i < h->nummendesnodes; i++)
  {
    h->x[i] = (double)random()/(double)RAND_MAX;
    h->y[i] = (double)random()/(double)RAND_MAX;
  }
}

double
mendes_energy (struct mendesgraph *h)
{
  double e1, e2;
  int in1, in2, ia;
  double dx, dy;

  e1 = e2 = 0.0;
  /*
   * e1 is the repulsive energy among all nodes
   */

  for (in1 = 0; in1 < h->nummendesnodes; in1++)
  {
    for (in2 = in1+1; in2 < h->nummendesnodes; in2++)
    {
      dx = h->x[in1] - h->x[in2];
      dy = h->y[in1] - h->y[in2];
      e1 += 1.0/(dx*dx + dy*dy);
    }
  }

  /*
   * e2 is the attractive energy of arcs
   */

  for (ia = 0; ia < h->nummendesarcs; ia++)
  {
    in1 = h->arcincplus[ia];
    in2 = h->arcincminus[ia];
    dx = h->x[in1] - h->x[in2];
    dy = h->y[in1] - h->y[in2];
    e2 += dx*dx + dy*dy;
  }

  return (e1*e2);
}
