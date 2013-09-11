/*
 * Canonification procedure for region description.  This is based
 * on a suggestion of Giovanni Paolini (my son).
 * Each connected component of the apparent contour has a dual graph
 * having a node for each region and an arc for each (extended) arc.
 * This dual graph is itself planar, with a base node (external region).
 *
 * the connected components of the apparent contour form a hierarchy
 * where each component is directly contained in a (holed) region of
 * a "parent" connected component, with the exception of the most external
 * component(s).
 *
 * Suppose for the moment that the apparent contour is connected (i.e.
 * all regions are simply connected (with the exception of the external
 * region).  Consider the dual graph.
 * If we select one of the arcs with the base vertex of the dual graph as
 * a vertex, taking advantage from the planarity we have a natural way
 * to construct a spanning tree based on the "depth-first-search" algorithm.
 *
 * following a (dual) arc from a (dual) node R1 to the dual node R2 (this
 * corresponds to going from region R1 to region R2 traversing a given arc
 * of the apparent contour) gives a natural way to select one of the arcs
 * of R2 (call it entry-point) and the planarity leads to a total ordering
 * of the arcs of R2.
 * Of course the "entry-point" depends on the starting region R1. The DFS
 * (depth-first-search) gives a way to naturally reach all nodes and hence to
 * define all the entry-points:  We enter a region starting from the external
 * region and following an arc, this defines the "entry-point" of that region.
 * Then from that region we recursively follow the arcs in the order indicated
 * by the entry-point, omitting those leading to already visited regions.
 *
 * This whole procedure only depends on the selection of one of the arcs
 * at the external region, which is arbitrary.
 *
 * We now have a standardized region description for each choice of an "external
 * arc" (arcs connected to the external region) and we select the one that
 * minimize the description with respect to some lexicografical ordering.
 * Such algorithm has a computational complexity O(n^2).  There are more efficient
 * algorithms available, however for our porposes this seems sufficient.
 *
 * Considering the general case of a nonconnected apparent contour is not
 * a real problem: we can canonify the most internal components first, then we have
 * a way to order (lexicographically) all the "holes" of a region which in turn
 * leads to a way to compare hierarchies of apparent contours.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "contour.h"
#include "giovecanonify.h"

extern int debug;
static int regioncount;
static short int *rmarks;
static short int *rmarks2;

/*
 * We need a way to mark visited regions. Instead of adding
 * a field to the region structure we prefer to use a pair
 * of integer vectors (rmarks and rmarks2)
 *
 * giovecanonify is the main entry point for the procedure
 */

void
giovecanonify (struct sketch *s)
{
  struct region *extregion, *r;
  struct borderlist *bl;
  struct arc *arc;
  int i, tag;

  if (s->arcs == 0) return;

  if (debug) printf ("canonify arcs...\n");
  if (s->isempty) return;
  for (arc = s->arcs; arc; arc = arc->next) canonifyarc (arc);

  extregion = s->regions;
  bl = extregion->border;
  
  assert (bl->isexternal && bl->sponda == 0);

  tag = 0; for (r = extregion; r; r = r->next) r->tag = tag++;
  s->regioncount = tag;  // sometimes this is not correct, e.g. contour removehole ...
  regioncount = s->regioncount;
  rmarks = (short int *) malloc (regioncount*sizeof(short int));
  rmarks2 = (short int *) malloc (regioncount*sizeof(short int));
  for (i = 0; i < regioncount; i++) rmarks[i] = rmarks2[i] = 0;

  bl->next = giovecanonifyblist (bl->next);

  for (i = 0; i < regioncount; i++)
  {
    assert (rmarks[i] == 0);
    assert (rmarks2[i] == 0);
  }
  free (rmarks);
  free (rmarks2);

  //giovepostcanonify (s);
  return;
}

void
giovepostcanonify (struct sketch *s)
{
  int tag;
  struct arc *arc;
  struct region *r;

  giove_sort_regions (s);
  /* rinumero le regioni */
  tag = 0; for (r = s->regions; r; r = r->next) r->tag = tag++;
  sortarcs (s);
  s->arcs = sortequivarcs (s->arcs);
  /* rinumero gli archi */
  tag = 1; for (arc = s->arcs; arc; arc = arc->next) arc->tag = tag++;
}

/*
 * canonification of a set of holes (of the same region). In particular
 * an apparent contour can be regarded as a set of one or more holes in
 * the external region.  Holes is a completely unordered set, so that
 * ordering them in a natural way is an important part of the procedure.
 * This is done lexicographically once we have canonified each hole.
 */

struct borderlist *
giovecanonifyblist (struct borderlist *bl)
{
  struct borderlist *blrest;

  if (bl == 0) return(0);

  blrest = giovecanonifyblist (bl->next);
  bl->next = 0;

  bl->sponda = giovecanonifyhole (bl->sponda);
  return (gioveinsertholeinblist (bl, blrest));
}

/*
 * This is the canonification procedure for a single hole.  It uses
 * the main "DFS" procedure "giove_normalize" and "giove_renormalize"
 * (see below) for each possible choice of
 * an optimal selection of indistinguishable external arcs, canonify
 * with respect to each of them and finally select the lexicographically
 * optimal choice.
 */

struct border *
giovecanonifyhole (struct border *entrypoint)
{
  struct list_of_borders *cell, *optimalborders;
  int res;

  optimalborders = extract_optimal_borders (entrypoint);

  assert (optimalborders);
  entrypoint = optimalborders->border;
  giove_normalize (entrypoint);  /* canonify with this given entrypoint */

  for (cell = optimalborders->next; cell; cell = cell->next)
  {
    res = giove_compare_holes (cell->border, entrypoint);
    if (res < 0) /* found a better entrypoint */
    {
      entrypoint = cell->border;
      giove_renormalize (entrypoint);  /* do not recurse inside other connected components */
    }
  }

  free_list_of_borders (optimalborders);
  return (entrypoint);
}

/*
 * lexicographical comparison among two different normalizations (depending to
 * the choice of the external arc).  It uses the dfs procedure and does not
 * require the region description to be normalized.
 * This is important, because otherwise we would need to duplicate the data
 * structure.
 */

int
giove_compare_holes (struct border *b1, struct border *b2)
{
  int res, savedmark1, savedmark2;
  struct border *cross1, *cross2;
  struct region *extregion1, *extregion2;

  /* we perform a DFS on the two graphs using the giove walking strategy */
  extregion1 = b1->border->region;
  extregion2 = b2->border->region;

  savedmark1 = rmarks[extregion1->tag];
  savedmark2 = rmarks2[extregion2->tag];

  rmarks[extregion1->tag] = 1;
  rmarks2[extregion2->tag] = 1;

  cross1 = crossriver (b1);
  cross2 = crossriver (b2);

  res = giove_compare_dfs (cross1, cross2);
  reset_marks (b1, rmarks);
  reset_marks (b2, rmarks2);
  rmarks[extregion1->tag] = savedmark1;
  rmarks2[extregion2->tag] = savedmark2;
  return (res);
}

/*
 * DFS comparison procedure.  Note that we need two sets of independent
 * markers since we are working on the same graph with different normalization
 * the comparison often end long before the two graphs are completely traversed,
 * leading to a proper subset of marked regions.  This should not impact negatively
 * on the terminal traversal made to clean up the markers.
 */

int
giove_compare_dfs (struct border *b1, struct border *b2)
{
  int res, tag1, tag2;
  struct borderlist *bl1, *bl2;
  struct region *region1, *region2;
  struct border *bb1, *bb2, *cr1, *cr2;

  region1 = b1->border->region;
  region2 = b2->border->region;

  assert (rmarks[region1->tag] == 0);
  rmarks[region1->tag] = 1;
  assert (rmarks2[region2->tag] == 0);
  rmarks2[region2->tag] = 1;

  assert (region1->border == b1->border);
  assert (region2->border == b2->border);

  /* if there are holes we must compare them */
  for (bl1 = region1->border->next, bl2 = region2->border->next;
       bl1;
       bl1 = bl1->next, bl2 = bl2->next)
  {
    if (bl2 == 0) return (1);  /* region 1 has more holes */
    res = giove_compare_holes (bl1->sponda, bl2->sponda);
    if (res != 0) return (res);
  }
  if (bl2) return (-1);  /* region 1 has less holes */

  bb1 = b1; bb2 = b2;
  do
  {
    if (bb1->orientation != bb2->orientation)
      return ((bb1->orientation < bb2->orientation)?(-1):1);
    res = arccmp (bb1->info, bb2->info, 0);
    if (res != 0) return (res);

    cr1 = crossriver (bb1);
    cr2 = crossriver (bb2);

    tag1 = cr1->border->region->tag;
    tag2 = cr2->border->region->tag;
    if (rmarks[tag1] != rmarks2[tag2])
    {
      if (rmarks[tag1] == 0) return (-1);
      assert (rmarks2[tag2] == 0);
      return (1);
    }
    if (rmarks[tag1] == 0)
    {
      res = giove_compare_dfs (cr1, cr2);
      if (res != 0) return (res);
    }

    bb1 = bb1->next;
    bb2 = bb2->next;
  } while (bb1 != b1 && bb2 != b2);
  if (bb2 != b2) return (-1);  /* region 1 has shorter boundary */
  if (bb1 != b1) return (1);  /* region 1 has longer boundary */

  return (0);
}

/*
 * giove_normalize and giove_renormalize differ in that the first one
 * implies also a recursive canonification of the holes of traversed
 * regions.  Since this can be done once and for all, all subsequent
 * normalizations do not need to do this recursive canonification
 * (giove_renormalize).  Both case are covered by "giove_normalize_common"
 */

void
giove_normalize (struct border *entrypoint)
{
  giove_normalize_common (entrypoint, 1);
}

void
giove_renormalize (struct border *entrypoint)
{
  giove_normalize_common (entrypoint, 0);
}

/*
 * This procedure actually modifies the region description
 * in that it changes the pointer from a borderlist struct
 * into the circular lists of "border" data, thus defining
 * the "entrypoint" for each region.
 * Once this is done, it is possible to recognize the arcs
 * that are part of the constructed spanning tree, allowing
 * the traversal of the graph without the necessity of a
 * marking strategy for the nodes.
 * This is e.g. used in the reordering of the regions made
 * at the end of the canonification procedure.
 */

void
giove_normalize_common (struct border *entrypoint, int canonify_holes)
{
  int savedmark;
  struct region *extregion;
  struct border *bb, *cr;

  extregion = entrypoint->border->region;
  if (debug) printf ("Normalize hole, extregion %d, arc %d with flag %d\n",
              extregion->tag, entrypoint->info->tag, canonify_holes);
  savedmark = rmarks[extregion->tag];

  rmarks[extregion->tag] = 1;

  bb = entrypoint;
  do
  {
    cr = crossriver (bb);
    giove_normalize_dfs (cr, canonify_holes);
    bb = bb->next;
  } while (bb != entrypoint);

  reset_marks (entrypoint, rmarks);
  rmarks[extregion->tag] = savedmark;
  return;
}

void
giove_normalize_dfs (struct border *b, int canonify_holes)
{
  struct region *region;
  struct border *bb, *cr;

  region = b->border->region;
  if (rmarks[region->tag]) return;
  rmarks[region->tag] = 1;

  if (debug) printf ("DFS: entering region %d through arc %d\n", region->tag, b->info->tag);
  b->border->sponda = b;
  if (canonify_holes) region->border->next = giovecanonifyblist (region->border->next);

  bb = b;
  do
  {
    cr = crossriver (bb);
    giove_normalize_dfs (cr, canonify_holes);
    bb = bb->next;
  } while (bb != b);
}

/*
 * The selection of an optimal "normalized" description
 * if made among the normalizations obtained after selecting
 * an external arc.
 * Any intrinsic method leading to a small set of "indistinguishable"
 * arcs can be used in order to reduce the computational complexity.
 *
 * It would be best to make an "intrinsic" ordering of
 * the arcs (using arccmp), identify the equivalent
 * classes of indistinguishable borders, identify the
 * equivalent classes with smallest cardinality and
 * return the equivalent class that comes first in the
 * arccmp ordering.
 *
 * Right now we just extract all borders...
 */

struct list_of_borders *
extract_optimal_borders (struct border *entrypoint)
{
  struct border *b;
  struct list_of_borders *cell;
  struct list_of_borders *prevcell = 0;

  b = entrypoint;
  do
  {
    cell = (struct list_of_borders *) malloc (sizeof (struct list_of_borders));
    cell->border = b;
    cell->next = prevcell;
    b = b->next;
    prevcell = cell;
  } while (b != entrypoint);

  return (prevcell);
}

/*
 * small procedure to free the list of candidate "entry-points"
 */

void
free_list_of_borders (struct list_of_borders *optimalborders)
{
  if (optimalborders == 0) return;
  free_list_of_borders (optimalborders->next);
  free (optimalborders);
}

/*
 * Ordering of the set of holes of some region makes use of the
 * comparison procedure.  It is done here recursively by inserting
 * each hole in the list of already ordered holes.
 */

struct borderlist *
gioveinsertholeinblist (struct borderlist *blhole, struct borderlist *bl)
{
  int res;
  struct borderlist *blrest;

  assert (blhole && blhole->next == 0);

  if (bl == 0) return (blhole);

  /* this comparison takes advantage of the canonification thus constructed */
  res = giove_compare_holes (blhole->sponda, bl->sponda);

  if (res <= 0)
  {
    blhole->next = bl;
    return (blhole);
  }

  blrest = bl->next;
  bl->next = gioveinsertholeinblist (blhole, blrest);
  return (bl);
}

/*
 * utility to extrace the opposite "sponda" across the arc
 * (same arc, but viewed by the opposite region)
 */

struct border *
crossriver (struct border *b)
{
  if (b->orientation < 0) return (b->info->regionleft);
  assert (b->orientation > 0);
  return (b->info->regionright);
}

/*
 * marks cannot be reset by simply blanking the markers vector, since
 * there could be some other DFS procedure in action.  We simply
 * redo the DFS strategy by reversing the meaning of the markers.
 */

void
reset_marks (struct border *b, short int *marks)
{
  int tag;
  struct border *bb, *cr;

  tag = b->border->region->tag;
  marks[tag] = 0;

  bb = b;
  do
  {
    cr = crossriver (bb);
    reset_marks_dfs (cr, marks);
    bb = bb->next;
  } while (bb != b);
}

void
reset_marks_dfs (struct border *b, short int *marks)
{
  int tag;
  struct border *cr, *bb;

  tag = b->border->region->tag;
  if (marks[tag] == 0) return;
  marks[tag] = 0;

  bb = b;
  do
  {
    cr = crossriver (bb);
    reset_marks_dfs (cr, marks);
    bb = bb->next;
  } while (bb != b);
}

/*
 * The canonification procedure produces a natural way
 * to order the regions (and number them), simply by listing
 * regions in the order they are visited in the DFS traversal
 */

static struct region *prevregion;

void
giove_sort_regions (struct sketch *s)
{
  prevregion = 0;
  giove_relink_regions (s->regions);
}

void
giove_relink_regions (struct region *r)
{
  struct border *b, *bb, *cr;
  struct borderlist *bl;

  r->next = 0;
  if (prevregion) prevregion->next = r;
  prevregion = r;
  assert (r->border);
  /* first we walk through the holes */
  for (bl = r->border->next; bl; bl = bl->next)
  {
    bb = bl->sponda;
    do
    {
      cr = crossriver (bb);
      if (cr->border->region->border->sponda == cr)
        giove_relink_regions (cr->border->region);
      bb = bb->next;
    } while (bb != bl->sponda);
  }
  b = r->border->sponda;
  if (b) /* altrimenti e' la regione esterna, che non ha bordo */
  {
    bb = b;
    do
    {
      cr = crossriver (bb);
      if (cr->border->region->border->sponda == cr)
        giove_relink_regions (cr->border->region);
      bb = bb->next;
    } while (bb != b);
  }
}

