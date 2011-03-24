/*
 * computation of the fundamental group of the interior of
 * the surface
 */

#include <assert.h>
#include "contour.h"
#include "fundamental.h"

extern int debug;

void
compute_fundamental (struct ccomplex *cc)
{
  int ccnum, gennum, n;
  struct ccomplexarc *arc;
  struct ccomplexnode *node;
  struct ccomplexcc *cccc;

  if (debug) printf ("Constructing spanning tree\n");
  ccnum = find_spanning_tree (cc);
  if (debug) printf ("Found %d connected components\n", ccnum);

  assert (ccnum >= 1);

  for (cccc = cc->cc; cccc; cccc = cccc->next)
  {
    if (ccnum > 1) printf ("Connected component %d:\n", cccc->tag);
    gennum = 0;
    for (n = 0; n < cc->arcdim; n++)
    {
      arc = cc->arcs + n;
      if (arc->type == CC_REMOVED) continue;
      if (arc->isinspanningtree) continue;
      node = cc->nodes + arc->enda;
      if (node->cc != cccc) continue;
      gennum++;
    }
    printf ("Found %d generators\n", gennum);
  }

  // ora: calcolo generatori (per ogni cc), poi bisogna aggiungere le relazioni
  printf ("Not implemented!\n");
}

/*
 * procedures for the semplification of the cell complex
 */

int
complex_collapse (struct ccomplex *cc)
{
  int goon = 1;
  int count = 0;

  while (goon)
  {
    goon = 0;
    goon += complex_collapse_faces (cc);
    goon += complex_collapse_arcs (cc);
    count += goon;
  }
  return (count);
}

int
complex_collapse_faces (struct ccomplex *cc)
{
  struct ccomplexface *faces = cc->faces;
  struct ccomplexarc *arcs = cc->arcs;
  int goon = 1;
  int count = 0;
  int n, i, *ivec;
  int narc;

  while (goon)
  {
    goon = 0;
    for (n = 0; n < cc->facedim; n++)
    {
      if (faces[n].type == CC_REMOVED) continue;
      for (i = 0; i < faces[n].facebordernum; i++)
      {
        ivec = faces[n].faceborder;
        narc = onarc2narc (ivec[i]);
        assert (arcs[narc].refcount >= 1);
        if (arcs[narc].refcount > 1) continue;
        /* can collapse face n with arc narc */
        goon = 1;
        count++;
        complex_remove_face (cc, n);
        complex_remove_arc (cc, narc);
        break;
      }
    }
  }
  if (debug) printf ("collapsed %d faces\n", count);
  return (count);
}

int
complex_collapse_arcs (struct ccomplex *cc)
{
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexnode *nodes = cc->nodes;
  int goon = 1;
  int count = 0;
  int n, nnode, rnode;

  while (goon)
  {
    goon = 0;
    for (n = 0; n < cc->arcdim; n++)
    {
      if (arcs[n].type == CC_REMOVED) continue;
      if (arcs[n].refcount > 0) continue;
      rnode = -1;
      nnode = arcs[n].enda;
      if (nodes[nnode].refcount == 1) rnode = nnode;
      nnode = arcs[n].endb;
      if (nodes[nnode].refcount == 1) rnode = nnode;
      if (rnode >= 0)
      {
        goon = 1;
        count++;
        complex_remove_arc (cc, n);
        complex_remove_node (cc, rnode);
        continue;
      }
    }
  }
  if (debug) printf ("collapsed %d arcs\n", count);
  return (count);
}

/*
 * melt together pairs of different faced that share
 * a common arc (border of no other face)
 */

int complex_facemelt (struct ccomplex *cc)
{
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexface *faces = cc->faces;
  struct ccomplexarc *arc;
  struct ccomplexface *face, *face1, *face2;
  int count = 0;
  int goon = 1;
  int n, m, i, nface1, nface2, i1, i2;
  int *ivec, *ivec1, *ivec2;

  while (goon)
  {
    goon = 0;
    for (n = 0; n < cc->arcdim; n++)
    {
      arc = arcs + n;
      if (arc->type == CC_REMOVED) continue;
      if (arc->refcount != 2) continue;
      nface1 = nface2 = -1;
      for (m = 0; m < cc->facedim; m++)
      {
        face = faces + m;
        if (face->type == CC_REMOVED) continue;
        ivec = face->faceborder;
        for (i = 0; i < face->facebordernum; i++)
        {
          if (onarc2narc(ivec[i]) == n)
          {
            if (nface1 < 0)
            {
              nface1 = m;
              i1 = i;
            } else {
              assert (nface2 < 0);
              nface2 = m;
              i2 = i;
            }
          }
        }
      }
      assert (nface1 >= 0 && nface2 >= 0);
      if (nface1 == nface2) continue;
      /*
       * We can melt the two faces together
       */
      face1 = faces + nface1;
      face2 = faces + nface2;
      ivec1 = face1->faceborder;
      ivec2 = face2->faceborder;
      if (ivec1[i1] < 0) cc_revert_face (cc, nface1);
      if (ivec2[i2] > 0) cc_revert_face (cc, nface2);
      goon = 1;
      count++;
      complex_do_melt_faces (cc, nface1, nface2, m);
      if (debug) cellcomplex_checkconsistency (cc);
    }
  }
  return (count);
}

/*
 * we know that arc is positively oriented in nface1 and
 * negatively oriented in nface2
 */

void
complex_do_melt_faces (struct ccomplex *cc, int nface1, int nface2, int narc)
{
  struct ccomplexface *face1, *face2;
  int *ivec1, *ivec2, *newivec;
  int found1 = 0;
  int found2 = 0;
  int i, j, k, j2;

  face1 = cc->faces + nface1;
  ivec1 = face1->faceborder;
  face2 = cc->faces + nface2;
  ivec2 = face2->faceborder;
  newivec = (int *) malloc (face1->facebordernum + face2->facebordernum - 2);

  k = 0;
  for (i = 0; i < face1->facebordernum; i++)
  {
    if (ivec1[i] == narc + 1)
    {
      found1++;
      assert (found1 == 1);
      for (j = 0; j < face2->facebordernum; j++)
      {
        if (ivec2[j] == - narc - 1)
        {
          found2++;
          assert (found2 == 1);
          j2 = j;
          continue;
        }
        if (found2) newivec[k++] = ivec2[j];
      }
      assert (found2);
      for (j = 0; j < j2; j++) newivec[k++] = ivec2[j];
    } else {
      newivec[k++] = ivec1[i];
    }
  }
  complex_remove_arc (cc, narc);
  face1->facebordernum += face2->facebordernum;
  free (ivec1);
  face1->faceborder = newivec;
  complex_remove_face (cc, nface2);
}

void
complex_remove_face (struct ccomplex *cc, int nface)
{
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexface *faces = cc->faces;
  int i, *ivec, narc;

  ivec = faces[nface].faceborder;

  for (i = 0; i < faces[nface].facebordernum; i++)
  {
    narc = onarc2narc (ivec[i]);
    arcs[narc].refcount--;
    assert (arcs[narc].refcount >= 0);
  }
  free (ivec);

  faces[nface].type = CC_REMOVED;
  cc->facenum--;
}

void
complex_remove_arc (struct ccomplex *cc, int narc)
{
  struct ccomplexnode *nodes = cc->nodes;
  struct ccomplexarc *arcs = cc->arcs;
  int nnode;

  assert (arcs[narc].refcount == 0);  // Cannot remove part of a face border!
  nnode = arcs[narc].enda;
  nodes[nnode].refcount--;
  assert (nodes[nnode].refcount >= 0);
  nnode = arcs[narc].endb;
  nodes[nnode].refcount--;
  assert (nodes[nnode].refcount >= 0);

  arcs[narc].type = CC_REMOVED;
  cc->arcnum--;
}

void
complex_remove_node (struct ccomplex *cc, int nnode)
{
  struct ccomplexnode *nodes = cc->nodes;

  assert (nodes[nnode].refcount == 0);  // Cannot remove an endpoint of some arc!
  nodes[nnode].type = CC_REMOVED;
  cc->nodenum--;
}

void
complex_countreferences (struct ccomplex *cc)
{
  int n;
  struct ccomplexnode *nodes = cc->nodes;
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexface *faces = cc->faces;
  int nnode, i, *ivec, narc;

  for (n = 0; n < cc->nodedim; n++)
  {
    if (nodes[n].type == CC_REMOVED) continue;
    nodes[n].refcount = 0;
  }

  for (n = 0; n < cc->arcdim; n++)
  {
    if (arcs[n].type == CC_REMOVED) continue;
    arcs[n].refcount = 0;
    nnode = arcs[n].enda;
    nodes[nnode].refcount++;
    nnode = arcs[n].endb;
    nodes[nnode].refcount++;
  }

  for (n = 0; n < cc->facedim; n++)
  {
    if (faces[n].type == CC_REMOVED) continue;
    for (i = 0; i < faces[n].facebordernum; i++)
    {
      ivec = faces[n].faceborder;
      narc = onarc2narc (ivec[i]);
      arcs[narc].refcount++;
    }
  }
}

/*
 * procedures for the construction of the cell complex
 */

struct ccomplex *
compute_cellcomplex (struct sketch *s, int fg_type)
{
  extern int finfinity;
  struct ccomplex *cc;
  int euler, surfeuler, realeuler;
  int status;

  if (finfinity != 0) fprintf (stderr, "Value of f at infinity (%d) must be zero\n", finfinity);
  assert (finfinity == 0);
  computefvalue (s, s->regions, 0 /* should be finfinity */);

  cc = (struct ccomplex *) malloc (sizeof (struct ccomplex));
  cc->type = fg_type;
  cc->sketch = s;

  cc->nodenum = cc->nodedim = fundamental_countnodes (s);
  cc->arcnum = cc->arcdim = fundamental_countarcs (cc->sketch, cc->type);
  cc->facenum = cc->facedim = fundamental_countfaces (cc->sketch, cc->type);
  euler = cc->nodenum - cc->arcnum + cc->facenum;
  if (debug)
  {
    printf ("Computing cell complex for the ");
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

  cc->nodes = (struct ccomplexnode *) malloc (cc->nodedim * sizeof (struct ccomplexnode));
  cc->arcs = (struct ccomplexarc *) malloc (cc->arcdim * sizeof (struct ccomplexarc));
  cc->faces = (struct ccomplexface *) malloc (cc->facedim * sizeof (struct ccomplexface));
  if (debug) printf ("Creating nodes\n");
  fundamental_fillnodes (cc);
  if (debug) printf ("Creating arcs\n");
  fundamental_fillarcs (cc);
  if (debug) printf ("Creating faces\n");
  fundamental_fillfaces (cc);

  complex_countreferences (cc);

  if (debug)
  {
    cellcomplex_print (cc, 2);
    status = cellcomplex_checkconsistency (cc);
    assert (status);
    printf ("Connected components: %d\n", find_spanning_tree (cc));
  }
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
  for (i = 0; i < cc->nodedim; i++) nodes[i].cc = 0;
  for (i = 0; i < cc->arcdim; i++) arcs[i].isinspanningtree = 0;

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
      for (i = 0; i < cc->arcdim; i++)
      {
        arc = arcs + i;
        if (arc->type == CC_REMOVED) continue;
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

  for (i = 0; i < cc->nodedim; i++)
  {
    if (cc->nodes[i].type == CC_REMOVED) continue;
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
  int vecdim = cc->nodedim;
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
  struct arc *a, *anext, *ase;
  int stratum, strata;
  int n1, n2, i, n, na, nb;
  int d, dne, dse, parity;
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
      for (n = 0; n < cc->nodedim; n++)
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
      assert (nb >= 0);
      sectiona = sectionb = 0;
      while (sectiona <= a->cusps)
      {
        while (sectionb <= a->cusps &&
          ((a->depths[sectiona] == stratum) ==
           (a->depths[sectionb] == stratum)) ) sectionb++;
        /* now sectiona and sectionb-1 indicate an arc */
        
        assert (vecpt < cc->arcs + cc->arcdim);
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
        vecpt->arc = a;
        vecpt->stratum = stratum;
        vecpt->cusp1 = sectiona;
        vecpt->cusp2 = sectionb;
        if (sectionb > a->cusps) vecpt->cusp2 = 0;
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
          assert (vecpt < cc->arcs + cc->arcdim);
          vecpt->type = CC_ARCTYPE_VIRTUAL;
          vecpt->enda = n1;
          vecpt->endb = n2;
          vecpt->stratum = 0;
          vecpt->bl = bl;
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
          assert (vecpt < cc->arcs + cc->arcdim);
          vecpt->type = CC_ARCTYPE_VIRTUAL;
          vecpt->enda = n1;
          vecpt->endb = n2;
          vecpt->stratum = stratum;
          vecpt->bl = bl;
          vecpt++;
        }
      }
    }
  }

  /*
   * we now generate all the columns
   */

  if (cc->type != FG_SURFACE)
  {
    for (a = s->arcs; a; a = a->next)
    {
      strata = a->regionright->border->region->f;
      /*
       * parity is the strata parity of the base of
       * columns, it will change when crossing nodes
       */
      parity = 0;
      if (cc->type == FG_EXTERNAL) parity = 1;
      if (a->endpoints == 0)
      {
        strata++;
        /* creating virtual columns */
        d = a->depths[0];
        for (stratum = 0; stratum < strata - 1; stratum++)
        {
          if (stratum == d) parity = 1 - parity;
          if ((stratum % 2) != parity) continue;
          na = fund_findnode (cc, a, stratum);
          nb = fund_findnode (cc, a, stratum+1);
          assert (na >= 0 && nb >= 0);
          assert (vecpt < cc->arcs + cc->arcdim);
          vecpt->enda = na;
          vecpt->endb = nb;
          vecpt->type = CC_ARCTYPE_VCOLUMN;
          vecpt->stratum = stratum;  /* of the column base */
          vecpt++;
        }
      } else {
        bord = a->regionright;
        assert (bord->orientation < 0);
        bord = bord->next;
        if (bord->orientation < 0) continue;
        /* canonically oriented node */
        ase = bord->info;
        dne = a->depths[0];
        dse = ase->depths[0];
        if (dne >= dse + 2)
          dne--;
         else
          dse++;
        for (stratum = 0; stratum < strata - 1; stratum++)
        {
          if (stratum == dne || stratum == dse) parity = 1 - parity;
          if ((stratum % 2) != parity) continue;
          na = fund_findnode (cc, a, stratum);
          nb = fund_findnode (cc, a, stratum+1);
          assert (na >= 0 && nb >= 0);
          assert (vecpt < cc->arcs + cc->arcdim);
          vecpt->enda = na;
          vecpt->endb = nb;
          vecpt->type = CC_ARCTYPE_COLUMN;
          vecpt->stratum = stratum;  /* of the column base */
          vecpt++;
        }
      }
    }
  }

  assert (vecpt == cc->arcs + cc->arcdim);
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
  if (bord->orientation < 0) d = a->depths[a->cusps];

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

  for (n = 0; n < cc->nodedim; n++)
  {
    node = cc->nodes + n;
    if (node->type < CC_NODETYPE_FOLD || node->type > CC_NODETYPE_VIRTUALCUT) continue;
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
 * compute faces
 */

/* local prototypes */
void fund_fillivec (struct ccomplex *cc, struct ccomplexface *face,
                    struct arc *a, int stratum, int cusp1, int cusp2);
int fund_findarc (struct ccomplex *cc, struct arc *a, int stratum, int cusp1, int cusp2);
int fund_findvarc (struct ccomplex *cc, struct borderlist *bl, int stratum);
int fund_findcolumn (struct ccomplex *cc, int nnode, int stratum);

void
fundamental_fillfaces (struct ccomplex *cc)
{
  int stratum, strata, arcnum, i, d, dcusp;
  struct ccomplexface *vecpt;
  struct region *r;
  struct borderlist *bl;
  struct border *b, *bp;
  struct arc *a;
  int ori, startcusp, astratum, na, *ivec;
  int parity, sectiona, sectionb, cusp1, cusp2;

  struct sketch *s = cc->sketch;

  vecpt = cc->faces;

  for (r = s->regions; r; r = r->next)
  {
    strata = r->f;
    if (r->border->sponda == 0) continue;
    if (strata == 0 && cc->type == FG_EXTERNAL) strata = 1;
    /* c'e' uno strato anche nei buchi */
    for (stratum = 0; stratum < strata; stratum++)
    {
      if ((stratum % 2) == 1 && cc->type == FG_INTERNAL) continue;
      if (stratum > 0 && (stratum % 2) == 0 && cc->type == FG_EXTERNAL) continue;
      assert (vecpt < cc->faces + cc->facedim);
      vecpt->type = CC_FACETYPE_HORIZONTAL;
      vecpt->stratum = stratum;
      /* count the number of arcs along the boundary */
      arcnum = 0;
      for (bl = r->border; bl; bl = bl->next)
      {
        b = bl->sponda;
        if (b->orientation < 0) b = b->next;
        /* if the base arc is negatively oriented, then the
         * base node is "after" this arc
         */
        bp = b;
        do
        {
          arcnum++;
          a = bp->info;
          for (i = 0; i < a->cusps; i++)
          {
            dcusp = a->depths[i];
            if (a->depths[i+1] < dcusp) dcusp = a->depths[i+1];
            if (bp->orientation > 0)
            {
              if (stratum == dcusp || stratum == dcusp+1 || stratum == dcusp+2) arcnum++;
            } else {
              if (stratum == dcusp) arcnum++;
            }
          }
          bp = bp->next;
        } while (bp != b);
        if (bl != r->border) arcnum += 2;
      }
      vecpt->facebordernum = arcnum;
      ivec = vecpt->faceborder = (int *) malloc (arcnum * sizeof (int));
      /* now actually create the border list */
      arcnum = 0;
      for (bl = r->border; bl; bl = bl->next)
      {
        if (bl != r->border)
        {
          na = fund_findvarc (cc, bl, stratum);
          ivec[arcnum++] = na+1;   // orientato positivamente
        }
        b = bl->sponda;
        if (b->orientation < 0)
        {
          b = b->next;
        }
        /* if the base arc is negatively oriented, then the
         * base node is "after" this arc
         */
        bp = b;

        do
        {
          a = bp->info;
          ori = bp->orientation;
          if (ori > 0)
          {
            startcusp = 0;
            astratum = stratum;
            if (stratum > a->depths[0]) astratum--;
            for (i = 1; i <= a->cusps; i++)
            {
              dcusp = a->depths[i-1];
              if (a->depths[i] < dcusp) dcusp = a->depths[i];
              if (stratum == dcusp || stratum == dcusp+1 || stratum == dcusp+2)
              {
                na = fund_findarc (cc, a, astratum, startcusp, i);
                assert (na >= 0);
                ivec[arcnum++] = na+1;
                astratum = stratum;
                if (stratum > a->depths[i]) astratum--;
                startcusp = i;
              }
            }
            na = fund_findarc (cc, a, astratum, startcusp, 0);
            ivec[arcnum++] = na+1;
          } else {
            startcusp = 0;
            astratum = stratum;
            if (stratum >= a->depths[a->cusps]) astratum++;
            for (i = a->cusps - 1; i >= 0; i--)
            {
              dcusp = a->depths[i+1];
              if (a->depths[i] < dcusp) dcusp = a->depths[i];
              if (stratum == dcusp)
              {
                assert (r->f > 0);     // be sure we cannot be here even for FG_EXTERNAL
                na = fund_findarc (cc, a, astratum, i+1, startcusp);
                assert (na >= 0);
                ivec[arcnum++] = -na-1;
                astratum = stratum;
                if (stratum >= a->depths[i]) astratum++;
                startcusp = i+1;
              }
            }
            if (cc->type == FG_EXTERNAL && r->f == 0) astratum = 0;  // this is a special case!
            na = fund_findarc (cc, a, astratum, 0, startcusp);
            assert (na >= 0);
            ivec[arcnum++] = -na-1;
          }
          bp = bp->next;
        } while (bp != b);
        if (bl != r->border)
        {
          na = fund_findvarc (cc, bl, stratum);
          ivec[arcnum++] = -na-1;   // orientato negativamente
        }
      }
      if ((stratum % 2) == 1) cc_revert_face (cc, vecpt - cc->faces);
      vecpt++;
    }
  }

  /*
   * now create the vertical walls
   * only for FG_INTERNAL and FG_EXTERNAL
   */
  if (cc->type != FG_SURFACE)
  {
    for (a = s->arcs; a; a = a->next)
    {
      strata = a->regionleft->border->region->f;  // strata with fold lines counting twice
      d = a->depths[0];
      parity = 0;
      if (cc->type == FG_EXTERNAL) parity = 1;
      for (stratum = 0; stratum < strata - 1; stratum++)
      {
        if ((stratum % 2) != parity) continue;
        sectionb = 0;
        for (sectiona = 0; sectiona <= a->cusps; sectiona++)
        {
          if (a->depths[sectiona] == stratum) continue;
          for (sectionb = sectiona + 1;
               a->depths[sectionb] != stratum && sectionb <= a->cusps;
               sectionb++);
          /* now sectiona-(sectionb-1) is a wall range */
          cusp1 = sectiona;
          cusp2 = sectionb;
          sectiona = sectionb;
          if (cusp2 > a->cusps) cusp2 = 0;
          astratum = stratum;
          if (stratum > a->depths[cusp1]) astratum--;
          /* devo creare un muro verticale */
          fund_fillivec (cc, vecpt, a, astratum, cusp1, cusp2);
          vecpt++;
        }
      }
    }
  }

  assert (vecpt == cc->faces + cc->facedim);
}

/* create integer vector with boundary of a vertical face */

void
fund_fillivec (struct ccomplex *cc, struct ccomplexface *face,
               struct arc *a, int stratum, int cusp1, int cusp2)
{
  int na, na1, na2;
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexarc *arc1, *arc2;
  struct ccomplexnode *nodes = cc->nodes;
  struct ccomplexnode *n1, *n2;
  int arcnum, *ivec;
  int i, *buffer1, *buffer2;
  int ib1, ib2, nnode;
  int cuspstart;

  buffer1 = (int *) malloc ((a->cusps + 1)*sizeof(int));
  buffer2 = (int *) malloc ((a->cusps + 1)*sizeof(int));

  ib1 = 0;
  cuspstart = cusp1;
  do {
    na1 = fund_findarc (cc, a, stratum, cuspstart, -1);  // find arc starting at cusp1
    assert (na1 >= 0);
    buffer1[ib1++] = na1;
    arc1 = arcs + na1;
    cuspstart = arc1->cusp2;
  } while (cuspstart != cusp2);

  ib2 = 0;
  cuspstart = cusp1;
  do {
    na2 = fund_findarc (cc, a, stratum + 1, cuspstart, -1);  // find arc starting at cusp1
    assert (na2 >= 0);
    buffer2[ib2++] = na2;
    arc2 = arcs + na2;
    cuspstart = arc2->cusp2;
  } while (cuspstart != cusp2);

  arcnum = ib1 + ib2;
  arc1 = arcs + buffer1[0];
  arc2 = arcs + buffer2[0];
  n1 = nodes + arc1->enda;
  n2 = nodes + arc2->enda;
  assert (n2->stratum >= n1->stratum);
  arcnum += n2->stratum - n1->stratum;

  arc1 = arcs + buffer1[ib1-1];
  arc2 = arcs + buffer2[ib2-1];
  n1 = nodes + arc1->endb;
  n2 = nodes + arc2->endb;
  assert (n2->stratum >= n1->stratum);
  arcnum += n2->stratum - n1->stratum;

  face->type = CC_FACETYPE_WALL;
  face->stratum = stratum;
  face->facebordernum = arcnum;
  ivec = face->faceborder = (int *) malloc (arcnum * sizeof (int));
  arcnum = 0;
  for (i = 0; i < ib1; i++)
    ivec[arcnum++] = buffer1[i] + 1;

  arc1 = arcs + buffer1[ib1-1];
  arc2 = arcs + buffer2[ib2-1];
  nnode = arc1->endb;
  n1 = nodes + nnode;
  n2 = nodes + arc2->endb;
  for (i = n1->stratum; i < n2->stratum; i++)
  {
    na = fund_findcolumn (cc, nnode, i);
    assert (na >= 0);
    ivec[arcnum++] = na + 1;
    nnode = arcs[na].endb;
  }
 
  for (i = ib2 - 1; i >= 0; i--)
    ivec[arcnum++] = - buffer2[i] - 1;

  arc1 = arcs + buffer1[0];
  arc2 = arcs + buffer2[0];
  nnode = arc2->enda;
  n1 = nodes + arc1->enda;
  n2 = nodes + nnode;
  for (i = n2->stratum - 1; i >= n1->stratum; i--)
  {
    na = fund_findcolumn (cc, nnode, i);
    assert (na >= 0);
    ivec[arcnum++] = - na - 1;
    nnode = arcs[na].enda;
  }

  free (buffer1);
  free (buffer2);
}

int
fund_findarc (struct ccomplex *cc, struct arc *a, int stratum, int cusp1, int cusp2)
{
  int n;
  struct ccomplexarc *arc;

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type != CC_ARCTYPE_CUT && arc->type != CC_ARCTYPE_FOLD) continue;
    if (arc->arc != a || arc->stratum != stratum) continue;
    if (arc->cusp1 == cusp1 && arc->cusp2 == cusp2) return (n);
    if (cusp2 < 0 && arc->cusp1 == cusp1) return (n);
    if (cusp1 < 0 && arc->cusp2 == cusp2) return (n);
  }
  fprintf (stderr, "Warning: cannot find face border for arc %d, stratum %d, ",
    a->tag, stratum);
  fprintf (stderr, "cusp1 %d, cusp2 %d\n", cusp1, cusp2);
  return (-1);
}

int
fund_findvarc (struct ccomplex *cc, struct borderlist *bl, int stratum)
{
  int n;
  struct ccomplexarc *arc;

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type != CC_ARCTYPE_VIRTUAL) continue;
    if (arc->bl == bl && arc->stratum == stratum) return (n);
  }
  fprintf (stderr, "Warning: cannot find face border for varc in region %d, stratum %d\n",
    bl->region->tag, stratum);
  return (-1);
}

int
fund_findcolumn (struct ccomplex *cc, int nnode, int stratum)
{
  int n;
  struct ccomplexarc *arc;

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type != CC_ARCTYPE_COLUMN && arc->type != CC_ARCTYPE_VCOLUMN) continue;
    if (arc->stratum != stratum) continue;
    if (arc->enda == nnode || arc->endb == nnode) return (n);
  }
  fprintf (stderr, "Warning: cannot find vertical column at node %d, stratum %d\n",
    nnode, stratum);
  return (-1);
}

/*
 * functions for ccomplex manipulation
 */

void
cc_revert_face (struct ccomplex *cc, int nface)
{
  struct ccomplexface *face = cc->faces;
  int *newvec, *oldvec;
  int i, size;

  face += nface;
  size = face->facebordernum;
  newvec = (int *) malloc (size * sizeof (int));
  oldvec = face->faceborder;

  for (i = 0; i < size; i++)
    newvec[size - i - 1] = -oldvec[i];

  free (oldvec);
  face->faceborder = newvec;
}

int
onarc2narc (int onarc)
{
  assert (onarc != 0);
  if (onarc < 0) onarc = -onarc;
  return (onarc - 1);
}

/*
 * functions for printing ccomplex content
 */

void
cellcomplex_printnodes (struct ccomplex *cc, int verbose)
{
  int n;
  struct ccomplexnode *node;

  for (n = 0; n < cc->nodedim; n++)
  {
    node = cc->nodes + n;
    if (node->type == CC_REMOVED) continue;
    printf ("node %d", n);
    if (verbose) printf (" ref %d", node->refcount);
    if (verbose >= 2)
    {
      if (node->type == CC_NODETYPE_CUSP) 
        printf (" of CUSP type, ne-arc %d, cusp %d", node->ne->tag, node->cusp);
       else
        printf (" of type %d, ne-arc %d, stratum %d", node->type, node->ne->tag, node->stratum);
    }
    printf ("\n");
  }
}

void
cellcomplex_printarcs (struct ccomplex *cc, int verbose)
{
  int n;
  struct ccomplexarc *arc;

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type == CC_REMOVED) continue;
    printf ("arc %d[%d,%d]", n, arc->enda, arc->endb);
    if (verbose) printf (" ref %d", arc->refcount);
    if (verbose >= 2)
    {
      printf (", stratum %d, type %d", arc->stratum, arc->type);
      if (arc->type == CC_ARCTYPE_CUT || arc->type == CC_ARCTYPE_FOLD)
        printf (" (%d-%d), tag %d", arc->cusp1, arc->cusp2, arc->arc->tag);
    }
    printf ("\n");
  }
}

void
cellcomplex_printfaces (struct ccomplex *cc, int verbose)
{
  int n, i, na;
  int *ivec;
  struct ccomplexface *face;

  for (n = 0; n < cc->facedim; n++)
  {
    face = cc->faces + n;
    if (face->type == CC_REMOVED) continue;
    ivec = face->faceborder;
    printf ("face %d [", n);
    for (i = 0; i < face->facebordernum; i++)
    {
      na = onarc2narc(ivec[i]);
      printf ("%c", (ivec[i]>0)?'+':'-');
      printf ("%d ", na);
    }
    printf ("]");
    if (verbose >= 2) printf (", stratum %d", face->stratum);
    printf ("\n");
  }
}

void
cellcomplex_print (struct ccomplex *cc, int verbose)
{
  cellcomplex_printnodes (cc, verbose);
  cellcomplex_printarcs (cc, verbose);
  cellcomplex_printfaces (cc, verbose);
}

/*
 * controlli di consistenza
 */

int
cellcomplex_checkconsistency (struct ccomplex *cc)
{
  struct ccomplexnode *nodes = cc->nodes;
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexface *faces = cc->faces;
  int count, n, nnode, i, inext, *ivec;
  int a1, a2, nodeto, nodefrom;

  for (count = 0, n = 0; n < cc->nodedim; n++)
    if (nodes[n].type != CC_REMOVED) count++;
  assert (count == cc->nodenum);

  for (count = 0, n = 0; n < cc->arcdim; n++)
  {
    if (arcs[n].type == CC_REMOVED) continue;
    count++;
    nnode = arcs[n].enda;
    assert (nnode >= 0 && nnode < cc->nodedim && nodes[nnode].type != CC_REMOVED);
    nnode = arcs[n].endb;
    assert (nnode >= 0 && nnode < cc->nodedim && nodes[nnode].type != CC_REMOVED);
  }
  assert (count == cc->arcnum);

  for (count = 0, n = 0; n < cc->facedim; n++)
  {
    if (faces[n].type == CC_REMOVED) continue;
    count++;
    assert (faces[n].facebordernum > 0);
    ivec = faces[n].faceborder;
    for (i = 0; i < faces[n].facebordernum; i++)
    {
      inext = i + 1;
      if (inext >= faces[n].facebordernum) inext = 0;
      a1 = onarc2narc (ivec[i]);
      assert (a1 >= 0 && a1 < cc->arcdim && arcs[a1].type != CC_REMOVED);
      a2 = onarc2narc (ivec[inext]);
      assert (a2 >= 0 && a2 < cc->arcdim && arcs[a2].type != CC_REMOVED);
      nodeto = arcs[a1].endb;
      if (ivec[i] < 0) nodeto = arcs[a1].enda;
      nodefrom = arcs[a2].enda;
      if (ivec[inext] < 0) nodefrom = arcs[a2].endb;
      if (nodeto != nodefrom)
      {
        printf ("Consistency error in face %d for contiguous arcs %c%d %c%d\n",
          n, (ivec[i]>0)?'+':'-', a1, (ivec[inext]>0)?'+':'-', a2);
        return (0);
      }
    }
  }
  assert (count == cc->facenum);

  return (1);
}

