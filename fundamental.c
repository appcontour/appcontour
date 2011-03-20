/*
 * computation of the fundamental group of the interior of
 * the surface
 */

#include <assert.h>
#include "contour.h"
#include "fundamental.h"

struct ccomplex *
compute_fundamental (struct sketch *s, int fg_type)
{
  extern int finfinity;
  extern int debug;
  struct ccomplex *cc;
  int ccnum;

  if (finfinity != 0) fprintf (stderr, "Value of f at infinity (%d) must be zero\n", finfinity);
  assert (finfinity == 0);
  computefvalue (s, s->regions, 0 /* should be finfinity */);

  cc = (struct ccomplex *) malloc (sizeof (struct ccomplex));
  cc->type = fg_type;
  cc->sketch = s;

  cc->nodenum = fundamental_countnodes (s);
  cc->arcnum = fundamental_countarcs (cc->sketch, cc->type);
  cc->facenum = fundamental_countfaces (cc->sketch, cc->type);
  if (debug) printf ("Complex nodes: %d\n", cc->nodenum);
  if (debug) printf ("Complex arcs: %d\n", cc->arcnum);
  if (debug) printf ("Complex faces: %d\n", cc->facenum);
  if (debug) printf ("Euler characteristic: %d\n", cc->nodenum + cc->facenum - cc->arcnum);
  cc->nodes = (struct ccomplexnode *) malloc (cc->nodenum * sizeof (struct ccomplexnode));
  cc->arcs = (struct ccomplexarc *) malloc (cc->arcnum * sizeof (struct ccomplexarc));
  fundamental_fillnodes (cc);
  fundamental_fillarcs (cc);

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
  struct arc *a;
  struct ccomplexnode *vecpt;
  struct ccomplexnode *vec = cc->nodes;
  struct sketch *s = cc->sketch;
  int fleft, strata, stratum;
  int vecdim = cc->nodenum;

  vecpt = vec;
  for (a = s->arcs; a; a = a->next)
  {
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
        vecpt->under = vecpt->over = a;
        vecpt++;
      }
    } else {
      assert (0);
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
      // fprintf (stderr, "Cusps are not allowed for now (arc %d: %d cusps)\n", a->tag, a->cusps);
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

  return (virtualnodes + normalnodes);
}

/*
 * fill the vector containing all the arcs of the complex
 */

void
fundamental_fillarcs (struct ccomplex *cc)
{
  struct arc *a;
  int stratum, strata;
  int n1, n2;
  struct ccomplexarc *vecpt;
  struct sketch *s = cc->sketch;

  vecpt = cc->arcs;

  for (a = s->arcs; a; a = a->next)
  {
    strata = a->regionright->border->region->f + 1;
    for (stratum = 0; stratum < strata; stratum++)
    {
      assert (vecpt < cc->arcs + cc->arcnum);
      vecpt->type = CC_ARCTYPE_CUT;
      if (stratum == a->depths[0]) vecpt->type = CC_ARCTYPE_FOLD;
      if (a->endpoints == 0)
      {
        printf ("devo cercare l'arco chiuso %d, strato %d\n", a->tag, stratum);
        n1 = n2 = 0;  // TODO!
        vecpt->enda = n1;
        vecpt->endb = n2;
      } else {
        printf ("devo cercare l'arco %d, strato %d\n", a->tag, stratum);
assert (0);
        n1 = n2 = 0;  // TODO!
        vecpt->enda = n1;
        vecpt->endb = n2;
      }
      vecpt++;
    }
  }

  assert (vecpt == cc->arcs + cc->arcnum);
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
  struct arc *a;
  struct region *r;
  struct borderlist *bl;
  int strata;

  for (a = s->arcs; a; a = a->next)
  {
    assert (a->cusps == 0);   /* no cusps allowed for now */
    cutarcs += a->regionright->border->region->f;
    foldarcs++;
  }

  for (r = s->regions; r; r = r->next)
  {
printf ("region: %d\n", r->tag);
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
printf ("virtualarcs 1: %d\n", virtualarcs);
        assert (r->border->next);
        for (bl = r->border->next->next; bl; bl = bl->next) virtualarcs++;
printf ("virtualarcs 2: %d\n", virtualarcs);
      }
    } else {
      if (fg_type != FG_SURFACE) strata /= 2;
      if (fg_type == FG_EXTERNAL) strata++;
printf ("virtualarcs 3: %d\n", virtualarcs);
      for (bl = r->border->next; bl; bl = bl->next) virtualarcs += strata;
printf ("virtualarcs 4: %d\n", virtualarcs);
    }
  }

printf ("foldarcs = %d, cutarcs = %d, virtualarcs = %d\n", foldarcs, cutarcs, virtualarcs);
  return (foldarcs + cutarcs + virtualarcs);
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
  struct arc *a;
  struct region *r;
  int strata, fright;
printf ("fg_type = %d\n", fg_type);
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
printf ("fright = %d, huffman=%d\n", fright, a->depths[0]);
      assert ((fright % 2) == 0);
      fright /= 2;
      assert (a->cusps == 0);
      if ((a->depths[0] % 2) == 1)
      {
        if (fg_type == FG_INTERNAL)
          fright += 1;
         else
          fright -= 1;
      }
      walls += fright;
    }
  }
printf ("orizfaces = %d, walls = %d\n", orizfaces, walls);
  return (orizfaces + walls);
}
