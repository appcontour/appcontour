/*
 * computation of the fundamental group of the interior of
 * the surface
 */

#include <assert.h>
#include "contour.h"
#include "fundamental.h"

/* local prototypes */
void fundamental_printarcs (struct ccomplex *cc);
void fundamental_printnodes (struct ccomplex *cc);

struct ccomplex *
compute_fundamental (struct sketch *s, int fg_type)
{
  extern int finfinity;
  extern int debug;
  struct ccomplex *cc;
  int ccnum, euler, surfeuler, realeuler;

  if (finfinity != 0) fprintf (stderr, "Value of f at infinity (%d) must be zero\n", finfinity);
  assert (finfinity == 0);
  computefvalue (s, s->regions, 0 /* should be finfinity */);

  cc = (struct ccomplex *) malloc (sizeof (struct ccomplex));
  cc->type = fg_type;
  cc->sketch = s;

  cc->nodenum = fundamental_countnodes (s);
  cc->arcnum = fundamental_countarcs (cc->sketch, cc->type);
  cc->facenum = fundamental_countfaces (cc->sketch, cc->type);
  euler = cc->nodenum - cc->arcnum + cc->facenum;
  if (debug)
  {
    printf ("Computing fundamental group for the ");
    switch (fg_type)
    {
      case FG_SURFACE: printf ("surface");
        break;
      case FG_INTERNAL: printf ("internal body");
        break;
      case FG_EXTERNAL: printf ("external body");
        break;
      default: printf ("(invalid choice: %d)", fg_type);
        break;
    }
    printf (".\n");
  }
  if (debug) printf ("Euler characteristic: %d = %d nodes - %d arcs + %d faces.\n",
             euler, cc->nodenum, cc->arcnum, cc->facenum);
  surfeuler = euler_characteristic (s);
  assert ((surfeuler % 2) == 0);
  switch (fg_type)
  {
    case FG_SURFACE:
      realeuler = surfeuler;
      break;
    case FG_INTERNAL:
      realeuler = surfeuler/2;
      break;
    case FG_EXTERNAL:
      realeuler = surfeuler/2 + 1;
      break;
    default:
      realeuler = -9999;
      break;
  }
  if (euler != realeuler) fprintf (stderr, 
     "WARNING: computed euler caracteristic (%d) differs from expected value (%d).\n",
     euler, realeuler);

  cc->nodes = (struct ccomplexnode *) malloc (cc->nodenum * sizeof (struct ccomplexnode));
  cc->arcs = (struct ccomplexarc *) malloc (cc->arcnum * sizeof (struct ccomplexarc));
  if (debug) printf ("Creating nodes\n");
  fundamental_fillnodes (cc);
  if (debug) printf ("Creating arcs\n");
  fundamental_fillarcs (cc);

if (debug) fundamental_printnodes (cc);
if (debug) fundamental_printarcs (cc);
  if (debug) printf ("Constructing spanning tree\n");
  ccnum = find_spanning_tree (cc);
  if (debug) printf ("Found %d connected components\n", ccnum);


  assert (ccnum >= 1);

  assert (ccnum == 1);  // only one connected component allowed for now

  // ora: calcolo generatori (per ogni cc), poi bisogna aggiungere le relazioni
  printf ("Not implemented!\n");
  return (cc);
}

/*
 * find a spanning tree for the 1D cell complex.  As a byproduct
 * we also have a separation of the various connected components
 * of the 1D skeleton, for each connected component a base node is
 * selected.
 */

/* local prototypes */

int cc_find_new_base_node (struct ccomplex *cc);

/* function definition */

int
find_spanning_tree (struct ccomplex *cc)
{
  struct ccomplexcc *cccc;
  struct ccomplexcc *lastcc = 0;
  struct ccomplexnode *nodes = cc->nodes;
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexarc *arc;
  int i, bn, again, n1, n2;

  /* initializations */
  cc->cc = 0;
  cc->ccnum = 0;
  for (i = 0; i < cc->nodenum; i++) nodes[i].cc = 0;
  for (i = 0; i < cc->arcnum; i++) arcs[i].isinspanningtree = 0;

  while ((bn = cc_find_new_base_node (cc)) >= 0)
  {
    cccc = (struct ccomplexcc *) malloc (sizeof (struct ccomplexcc));
    cccc->tag = cc->ccnum;
    cccc->basenode = bn;
    cc->ccnum++;
    if (lastcc == 0)
      cc->cc = cccc;
     else
      lastcc->next = cccc;
    lastcc = cccc;
    nodes[bn].cc = cccc;
    again = 1;
    while (again)   // repeated loop on arcs to expand spanning tree
    {
      again = 0;
      for (i = 0; i < cc->arcnum; i++)
      {
        arc = arcs + i;
        n1 = arc->enda;
        n2 = arc->endb;

        if (nodes[n1].cc == nodes[n2].cc) continue;
        again = 1;
        arc->isinspanningtree = 1;
        if (nodes[n1].cc)
        {
          assert (nodes[n2].cc == 0);
          nodes[n2].cc = nodes[n1].cc;
        } else {
          assert (nodes[n1].cc == 0);
          nodes[n1].cc = nodes[n2].cc;
        }
      }
    }
  }
  return (cc->ccnum);
}

int
cc_find_new_base_node (struct ccomplex *cc)
{
  int i;

  for (i = 0; i < cc->nodenum; i++)
  {
    if (cc->nodes[i].cc == 0) return (i);
  }
  return (-1);
}

/*
 * fill the vector containing all the nodes of the complex
 */

void
fundamental_fillnodes (struct ccomplex *cc)
{
  struct arc *a, *ane, *ase;
  struct border *bord;
  struct ccomplexnode *vecpt;
  struct ccomplexnode *vec = cc->nodes;
  struct sketch *s = cc->sketch;
  int fleft, strata, stratum;
  int vecdim = cc->nodenum;
  int dne, dse, i;

  vecpt = vec;
  for (a = s->arcs; a; a = a->next)
  {
    if (a->cusps > 0)
    {
      for (i = 1; i <= a->cusps; i++)
      {
        assert (vecpt < vec + vecdim);
        vecpt->type = CC_NODETYPE_CUSP;
        vecpt->stratum = a->depths[i];
        vecpt->ne = vecpt->se = a;
        vecpt->cusp = i;
        vecpt++;
      }
    }
    if (a->endpoints == 0)
    {
      fleft = a->regionleft->border->region->f;
      strata = fleft -1;
      for (stratum = 0; stratum < strata; stratum++)
      {
        assert (vecpt < vec + vecdim);
        vecpt->type = CC_NODETYPE_VIRTUALCUT;
        if (stratum == a->depths[0]) vecpt->type = CC_NODETYPE_VIRTUALFOLD;
        vecpt->stratum = stratum;
        vecpt->ne = vecpt->se = a;
        vecpt++;
      }
    } else {
      bord = a->regionright;
      assert (bord->orientation < 0);
      bord = bord->next;
      if (bord->orientation < 0) continue;
      /* archi in posizione canonica */
      ane = a;
      ase = bord->info;
      strata = ane->regionright->border->region->f;
      dne = ane->depths[0];
      dse = ase->depths[0];
      if (dne >= dse + 2)
        dne--;
       else
        dse++;
      for (stratum = 0; stratum < strata; stratum++)
      {
        assert (vecpt < vec + vecdim);
        vecpt->type = CC_NODETYPE_CUT;
        if (stratum == dse || stratum == dne) vecpt->type = CC_NODETYPE_FOLD;
        vecpt->stratum = stratum;
        vecpt->ne = ane;
        vecpt->se = ase;
        vecpt++;
      }
    }
  }
  assert (vecpt == vec + vecdim);
}

/*
 * count the number of nodes in the structure of cells
 * (we remove the 3D cells by deformation retraction)
 */

int
fundamental_countnodes (struct sketch *s)
{
  int virtualnodes = 0;
  int normalnodes = 0;
  // int nodesonfolds = 0;
  // int nodesonstrata = 0;
  int cuspnodes = 0;
  int fleft, fright;
  int strata;
  struct arc *a;
  // struct region *r;
  // struct borderlist *bl;
  // struct border *b, *bp;

  for (a = s->arcs; a; a = a->next)
  {
    if (a->cusps > 0)
    {
      cuspnodes += a->cusps;
    }
    fleft = a->regionleft->border->region->f;
    fright = a->regionright->border->region->f;
    assert (fleft - fright == 2);
    strata = fleft -1;
    if (a->endpoints != 0)
    {
      normalnodes += strata;
    } else {
      virtualnodes += strata;
    }
  }

  //for (r = s->regions; r; r = r->next)
  //{
  //  strata = r->f;
  //  if (strata == 0) continue;
  //  for (bl = r->border; bl; bl = bl->next)
  //  {
  //    bp = b = bl->sponda;
  //    do
  //    {
  //      if (bp->orientation < 0 && bp->next->orientation < 0) nodesonstrata += strata;
  //      bp = bp->next;
  //    } while (bp != b);
  //  }
  //}

  assert ((normalnodes % 2) == 0);  // normal nodes are counted twice!
  normalnodes /= 2;

  return (virtualnodes + normalnodes + cuspnodes);
}

/*
 * fill the vector containing all the arcs of the complex
 */

/*
 * local prototypes
 */

int stratum_start (struct arc *a, int stratum);
int stratum_end (struct arc *a, int stratum);
int stratum_varcend (struct border *bord, int stratum);

void
fundamental_fillarcs (struct ccomplex *cc)
{
  struct arc *a, *anext;
  int stratum, strata;
  int n1, n2, i, n, na, nb;
  int sectiona, sectionb;
  struct ccomplexarc *vecpt;
  struct ccomplexnode *node;
  struct sketch *s = cc->sketch;
  struct border *bord;
  struct region *r;
  struct borderlist *bl;
  int *cuspnodes = 0;

  vecpt = cc->arcs;

  for (a = s->arcs; a; a = a->next)
  {
    strata = a->regionright->border->region->f + 1;
    if (a->cusps > 0) cuspnodes = (int *) malloc ((a->cusps+1)*sizeof(int));
    for (i = 1; i <= a->cusps; i++)
    { /*
       * cerco i nodi corrispondenti alle cuspidi, memorizzati
       * nel vettore cuspnodes a partire dall'indice 1
       */
      for (n = 0; n < cc->nodenum; n++)
      {
        node = cc->nodes + n;
        if (node->type == CC_NODETYPE_CUSP && node->ne == a)
        {
          assert (node->cusp <= a->cusps);
          cuspnodes[node->cusp] = n;
        }
      }
    }

    for (stratum = 0; stratum < strata; stratum++)
    { /*
       * devo trovare i due nodi estremi dell'arco
       */
      na = nb = -1;
      if (a->endpoints == 0)
      {
        na = nb = fund_findnode (cc, a, stratum);
      } else {
        bord = a->regionleft->next;
        bord = gettransborder (bord);
        anext = bord->next->info;
        na = fund_findnode (cc, a, stratum_start (a, stratum));
        nb = fund_findnode (cc, anext, stratum_end (a, stratum));
      }
      assert (na >= 0);
//printf ("a: %d, anext: %d, stratum %d\n", a->tag, anext->tag, stratum);
      assert (nb >= 0);
      sectiona = sectionb = 0;
      while (sectiona <= a->cusps)
      {
        while (sectionb <= a->cusps &&
          ((a->depths[sectiona] == stratum) ==
           (a->depths[sectionb] == stratum)) ) sectionb++;
        /* now sectiona and sectionb-1 indicate an arc */
        
        assert (vecpt < cc->arcs + cc->arcnum);
        vecpt->type = CC_ARCTYPE_CUT;
        if (stratum == a->depths[sectiona]) vecpt->type = CC_ARCTYPE_FOLD;

        n1 = na;
        if (sectiona > 0) n1 = cuspnodes[sectiona];
        n2 = nb;
        if (sectionb <= a->cusps) n2 = cuspnodes[sectionb];
        assert (n1 >= 0);
        assert (n2 >= 0);
        vecpt->enda = n1;
        vecpt->endb = n2;
        vecpt++;
        sectiona = sectionb;
      }
    }
    if (cuspnodes) {free (cuspnodes); cuspnodes = 0;}
  }

  for (r = s->regions; r; r = r->next)
  {
    strata = r->f;
    assert ((strata % 2) == 0);
    if (r->border->sponda == 0)
    {
      assert (strata == 0);
      if (cc->type == FG_EXTERNAL)
      {
        /*
         * in this case we need to create virtual arcs connecting
         * the connected components to each other, they are
         * one less than the number of connected components
         */
        assert (r->border->next);
        n1 = fund_findnode (cc, r->border->next->sponda->info, 0);
        assert (n1 >= 0);
        for (bl = r->border->next->next; bl; bl = bl->next)
        {
          n2 = fund_findnode (cc, bl->sponda->info, 0);
          assert (n2 >= 0);
          assert (vecpt < cc->arcs + cc->arcnum);
          vecpt->type = CC_ARCTYPE_VIRTUAL;
          vecpt->enda = n1;
          vecpt->endb = n2;
          vecpt++;
        }
      }
    } else {
      for (stratum = 0; stratum < strata; stratum++)
      {
        if ((stratum % 2) == 1 && cc->type == FG_INTERNAL) continue;
        if (stratum > 0 && (stratum % 2) == 0 && cc->type == FG_EXTERNAL) continue;
        n1 = fund_findnode (cc, r->border->sponda->info,
                            stratum_varcend (r->border->sponda, stratum));
        assert (n1 >= 0);
        for (bl = r->border->next; bl; bl = bl->next)
        {
          n2 = fund_findnode (cc, bl->sponda->info,
                              stratum_varcend (bl->sponda, stratum));
          assert (n2 >= 0);
          assert (vecpt < cc->arcs + cc->arcnum);
          vecpt->type = CC_ARCTYPE_VIRTUAL;
          vecpt->enda = n1;
          vecpt->endb = n2;
          vecpt++;
        }
      }
    }
  }

  assert (vecpt == cc->arcs + cc->arcnum);
}

int stratum_start (struct arc *a, int stratum)
{
  int d2tilde, d2delta;
  struct arc *ta;
  struct border *bord;

  if (a->endpoints == 0) return (stratum);

  bord = a->regionright->next;
  ta = bord->info;
  d2delta = get_d_increase_across_node (ta, -bord->orientation);
  assert ( d2delta == 0 || d2delta == 2);
  d2delta /= 2;
  if (bord->orientation > 0)
  {
    d2tilde = ta->depths[0] + d2delta;
    if (stratum > d2tilde) stratum--;
  } else {
    d2tilde = ta->depths[ta->cusps] + d2delta;
    if (stratum >= d2tilde) stratum++;
  }
  return (stratum);
}

int stratum_end (struct arc *a, int stratum)
{
  int d2tilde, d2delta;
  struct arc *ta;
  struct border *bord;

  if (a->endpoints == 0) return (stratum);

  bord = a->regionleft->next;
  ta = bord->info;
  d2delta = get_d_increase_across_node (ta, -bord->orientation);
  assert ( d2delta == 0 || d2delta == -2);
  d2delta /= 2;
  if (bord->orientation > 0)
  {
    d2tilde = ta->depths[0] + d2delta;
    if (stratum > d2tilde) stratum--;
  } else {
    d2tilde = ta->depths[ta->cusps] + d2delta;
    if (stratum >= d2tilde) stratum++;
  }
  return (stratum);
}

/*
 * this function computes the stratum number for one of the endpoints
 * of a virtual arc added to make a region simply connected.
 * the stratum in input is that of the virtual arc
 */

int
stratum_varcend (struct border *bord, int stratum)
{
  struct arc *a;
  int d;

  a = bord->info;
  d = a->depths[0];

  if (bord->orientation > 0 && stratum > d) stratum--;
  if (bord->orientation < 0 && stratum >= d) stratum++;
  /* this is the stratum of the corresponding arc */

  if (a->endpoints == 0) return (stratum);

  return (stratum_start (a, stratum));
}

/*
 * =================================
 */

int
fund_findnode (struct ccomplex *cc, struct arc *a, int stratum)
{
  int n;
  struct ccomplexnode *node;

  for (n = 0; n < cc->nodenum; n++)
  {
    node = cc->nodes + n;
    if ((node->ne == a || node->se == a) && node->stratum == stratum) return (n);
  }
  fprintf (stderr, "Warning: cannot find arc endpoint in the cell complex for arc %d, stratum %d\n",
    a->tag, stratum);
  return (-1);
}

/*
 * count the number of arcs (1D cells) in the structure of cells
 * (we remove the 3D cells by deformation retraction, also we
 * remove virtual arcs on top of virtual walls, which are themselves
 * removed)
 */

int
fundamental_countarcs (struct sketch *s, int fg_type)
{
  int foldarcs = 0;
  int cutarcs = 0;
  int virtualarcs = 0;
  int cusparcs = 0;
  int columns = 0;
  int vcolumns = 0;
  int acnodes = 0;
  struct arc *a;
  struct region *r;
  struct borderlist *bl;
  int strata, fright, d;

  for (a = s->arcs; a; a = a->next)
  {
    cutarcs += a->regionright->border->region->f;
    foldarcs++;
    cusparcs += 2*a->cusps;
    if (fg_type != FG_SURFACE)
    {
      /* count columns */
      fright = a->regionright->border->region->f;
      assert ((fright % 2) == 0);
      fright /= 2;
      d = (a->depths[0] % 2);
      if (fg_type == FG_EXTERNAL) d = -d;
      if (a->endpoints == 0)
      {
        vcolumns += fright + d;
      } else {
        columns += fright + 2*d;
        acnodes++;
      }
    }
  }

  if (fg_type != FG_SURFACE)
  {
    assert ((acnodes % 2) == 0);
    acnodes /= 2;
    if (fg_type == FG_INTERNAL)
      columns -= acnodes;
     else
      columns += acnodes;
    assert ((columns % 2) == 0);
    columns /= 2;
  }

  for (r = s->regions; r; r = r->next)
  {
    strata = r->f;
    assert ((strata % 2) == 0);
    if (r->border->sponda == 0)
    {
      assert (strata == 0);
      if (fg_type == FG_EXTERNAL)
      {
        /*
         * in this case we need to add virtual arcs connecting
         * the connected components to each other, they are
         * one less than the number of connected components
         */
        assert (r->border->next);
        for (bl = r->border->next->next; bl; bl = bl->next) virtualarcs++;
      }
    } else {
      if (fg_type != FG_SURFACE) strata /= 2;
      if (fg_type == FG_EXTERNAL) strata++;
      for (bl = r->border->next; bl; bl = bl->next) virtualarcs += strata;
    }
  }

  return (foldarcs + cutarcs + virtualarcs + cusparcs + columns + vcolumns);
}

/*
 * count the number of faces (2D cells) in the structure of cells
 * (we remove the 3D cells by deformation retraction, also we
 * remove virtual vertical walls together with their top arc
 * in case of fg_internal/fg_external)
 */

int
fundamental_countfaces (struct sketch *s, int fg_type)
{
  int orizfaces = 0;
  int walls = 0;
  int cuspwalls = 0;
  struct arc *a;
  struct region *r;
  int strata, fright, i, dmod2;
  /*
   * orizontal faces (no virtual walls because we deformation-retract them)
   */

  for (r = s->regions; r; r = r->next)
  {
    strata = r->f;
    assert ((strata % 2) == 0);
    if (r->border->sponda == 0)
    {
      assert (strata == 0);
      continue;
    }
    switch (fg_type)
    {
      case FG_SURFACE:
        orizfaces += strata;
      break;

      case FG_INTERNAL:
        orizfaces += strata/2;
      break;

      case FG_EXTERNAL:
        orizfaces += strata/2 + 1;
      break;

      default:
        assert (0);
      break;
    }
  }

  /*
   * vertical walls (not for FG_SURFACE)
   */

  if (fg_type != FG_SURFACE)
  {
    for (a = s->arcs; a; a = a->next)
    {
      fright = a->regionright->border->region->f;
      assert ((fright % 2) == 0);
      fright /= 2;
      if ((a->depths[0] % 2) == 1)
      {
        if (fg_type == FG_INTERNAL)
          fright += 1;
         else
          fright -= 1;
      }
      for (i = 1; i <= a->cusps; i++)
      {
        dmod2 = a->depths[i] % 2;
        if (dmod2 == 1 && fg_type == FG_INTERNAL) cuspwalls++;
        if (dmod2 == 0 && fg_type == FG_EXTERNAL) cuspwalls++;
      }
      walls += fright;
    }
  }
  return (orizfaces + walls + cuspwalls);
}

/*
 * ============
 */

void
fundamental_printnodes (struct ccomplex *cc)
{
  int n;
  struct ccomplexnode *node;

  for (n = 0; n < cc->nodenum; n++)
  {
    node = cc->nodes + n;
    if (node->type == CC_NODETYPE_CUSP) 
      printf ("node %d of CUSP type, ne-arc %d, cusp %d\n", n, node->ne->tag, node->cusp);
     else
      printf ("node %d of type %d, ne-arc %d, stratum %d\n", n, node->type, node->ne->tag, node->stratum);
  }
}

void
fundamental_printarcs (struct ccomplex *cc)
{
  int n;
  struct ccomplexarc *arc;

  for (n = 0; n < cc->arcnum; n++)
  {
    arc = cc->arcs + n;
    printf ("arc %d[%d,%d], type %d\n", n, arc->enda, arc->endb, arc->type);
  }
}
