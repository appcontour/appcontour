#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "contour.h"
#include "embedding.h"
#include "fundamental.h"
#include "parser.h"
#include "readdtcode.h"
#include "wirtinger.h"

extern struct global_data globals;

extern int debug;
extern int verbose;
extern int quiet;

int print_sketch_ra3 (int only3, int count, struct embedding *emb, struct sketch *s, int *nodemark);

struct sketch *
embedding2sketch (struct embedding *emb)
{
  int i, j, found;
  struct emb_node *node;
  struct sketch *sketch;
  struct arc *arc;
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct vecofintlist *loiv;
  struct dualembedding *dual;
  struct dual_region *dregion;
  int fourvalentnum, inode;
  int choice, bit, count;
  int *nodemark;
  int arcid, nodenum;

  if (emb->k == 0)
  {
    if (debug) printf ("This is a standard knot/link, converting to gausscode\n");
    loiv = embeddingtoloiv (emb);
    if (debug) printloiv (loiv);
    freeembedding (emb);
    sketch = readgausscodeloiv (loiv);
    // freeloiv (loiv); // loiv is freed by readgausscodeloiv!
    return (sketch);
  }

  for (i = 0; i < emb->k + emb->n; i++)
  {
    node = &emb->nodes[i];
    if (i < emb->k) assert (node->valency == 3);
     else assert (node->valency == 4);
  }

  dual = embedding2dual (emb);

  /*
   * WARNING: this code DOES NOT WORK when the final region description contains
   * closed arcs (homeomorphic to an S^1), this can happen when there is a region
   * in the region description that is bounded solely by tri-valent vertices.
   *
   */
    
  found = 0;
  for (dregion = dual->regions; dregion; dregion = dregion->next)
  {
    fourvalentnum = 0;
    for (j = 0; j < dregion->valency; j++)
    {
      /*
       * check if the corresponding node is 4-valent
       */
      inode = dregion->wedgeij[j]/4;
      node = &(emb->nodes[inode]);
      if (node->valency == 4) fourvalentnum++;
    }
    if (fourvalentnum == 0) found++;
    if (debug) printf ("region %d with %d fourvalent bounding nodes (out of %d)\n", dregion->id, fourvalentnum, dregion->valency);
  }

  if (found > 0)
  {
    printf ("\n# WARNING: presence of Regions with no 4-valent bound vertex!\n");
    printf ("# the code at present cannot cope with this situation.\n");
    printf ("#\n# However in this case the embedding is IH equivalent to a one-sum\n");
    printf ("# so we are not really interested in such a case.\n");
    return (0);
  }

  sketch = newsketch ();
  sketch->huffman_labelling = 1;
  sketch->arcs = 0;
  sketch->arccount = 0;

  /*
   * arcs that bound embedding regions.
   * Endpoints are four-valent nodes
   *
   * -- Numbering --
   * each crossing i of the diagram originates 4 crossings in the apparent contour
   * numbered (i.j) j = 0,1,2,3, in counterclockwise order and is the endpoint
   * of some other arc.
   *
   * We fix ideas deciding that if the departing arc goes up, then the arriving arc
   * arrives from the right.
   *
   * the numbering of these (long) args is then computed as (k being the number of 3-valent
   * nodes):
   *
   *   4*(i - k) + j + 1
   *
   * the last of these arcs will be numbered 4n
   *
   * This choice of the orientation leads to arcs oppositely oriented with respect to the
   * regions (f=0).
   */

  /*
   * arcs bounding the small quadrilateral crossing regions
   * (four per node), listed counter-clockwise.
   * the label value defines the under/over choice and can
   * be driven using the bits in the "choice" argument
   */

  count = 4*emb->n+1;
  choice = emb->choice;
  sketch->arccount += 4*emb->n;
  for (i = 1; i <= emb->n; i++)
  {
    bit = choice % 2;
    choice /= 2;

    arc = newarc (sketch);
    arc->tag = count++;
    if (arc->next)
    {
      sketch->arcs = arc->next; 
      arc->next = 0;
      insert_arc_in_list (arc, sketch->arcs);
    }
    arc->depths = (int *) malloc (sizeof (int));
    arc->depthsdim = 1;
    arc->cusps = 0;
    arc->depths[0] = 2*bit;
    arc->endpoints = 2;

    arc = newarc (sketch);
    arc->tag = count++;
    sketch->arcs = arc->next; 
    arc->next = 0;
    insert_arc_in_list (arc, sketch->arcs);
    arc->depths = (int *) malloc (sizeof (int));
    arc->depthsdim = 1;
    arc->cusps = 0;
    arc->depths[0] = 2*(1 - bit);
    arc->endpoints = 2;

    arc = newarc (sketch);
    arc->tag = count++;
    sketch->arcs = arc->next; 
    arc->next = 0;
    insert_arc_in_list (arc, sketch->arcs);
    arc->depths = (int *) malloc (sizeof (int));
    arc->depthsdim = 1;
    arc->cusps = 0;
    arc->depths[0] = 2*bit;
    arc->endpoints = 2;

    arc = newarc (sketch);
    arc->tag = count++;
    sketch->arcs = arc->next; 
    arc->next = 0;
    insert_arc_in_list (arc, sketch->arcs);
    arc->depths = (int *) malloc (sizeof (int));
    arc->depthsdim = 1;
    arc->cusps = 0;
    arc->depths[0] = 2*(1 - bit);
    arc->endpoints = 2;
  }

  sketch->arccount += 4*emb->n;
  for (i = 4*emb->n; i > 0; i--)
  {
    arc = newarc (sketch);
    arc->tag = i;
    arc->depths = (int *) malloc (sizeof (int));
    arc->depthsdim = 1;
    arc->cusps = 0;
    arc->depths[0] = 0;
    arc->endpoints = 2;
  }

  count = 0;
  /*
   * start with real regions (f = 0)
   */
  for (dregion = dual->regions; dregion; dregion = dregion->next)
  {
    int jr;

    region = newregion_tail (sketch);
    assert (region->border == 0);
    assert (region->tag == count);
    region->f = 0;
    if (count == 0)
    {
      bl = newborderlist (region);
      bl->isexternal = 1;
    }
    bl = newborderlist (region);

    blast = 0;
    for (jr = 0; jr < dregion->valency; jr++)
    {
      int iv, jv;

      iv = dregion->wedgeij[jr] / 4;
      node = &(emb->nodes[iv]);
      if (node->valency <= 3) continue;
      jv = dregion->wedgeij[jr] % 4;
      jv++; if (jv >= 4) jv -= 4;
      arcid = 4*(iv - emb->k) + jv + 1;
      b = newborder (bl);
      if (blast == 0)
      {
        bl->sponda = b;
      } else {
        b->next = bl->sponda;
        blast->next = b;
      }
      blast = b;
      b->orientation = -1;

      for (arc = sketch->arcs; arc; arc = arc->next)
      {
        if (arc->tag == arcid)
        {
          b->info = arc;
          break;
        }
      }
    }
    count++;
  }

  nodenum = emb->k + emb->n;
  nodemark = (int *) malloc (nodenum*sizeof(int));
  for (i = 0; i < nodenum; i++) nodemark[i] = 0;

  /*
   * now the elongated regions containing the triple-points
   * (one or two)
   */

  while (print_sketch_ra3 (1, count, emb, sketch, nodemark)) count++;

  while (print_sketch_ra3 (0, count, emb, sketch, nodemark)) count++;

  free (nodemark);

  /*
   * finally regions corresponding to four-valent nodes
   */

  for (i = 0; i < emb->n; i++)
  {
    region = newregion_tail (sketch);
    assert (region->border == 0);
    assert (region->tag == count);
    region->f = 4;
    bl = newborderlist (region);
    blast = 0;
    count++;
    for (j = 0; j < 4; j++)
    {
      b = newborder (bl);
      if (blast == 0)
      {
        bl->sponda = b;
      } else {
        b->next = bl->sponda;
        blast->next = b;
      }
      blast = b;
      b->orientation = 1;

      arcid = 4*emb->n + 4*i + j + 1;
      for (arc = sketch->arcs; arc; arc = arc->next)
      {
        if (arc->tag == arcid)
        {
          b->info = arc;
          break;
        }
      }
    }
  }

  if (verbose) printf ("Done converting embedding with: k = %d, n = %d, choice = %d into an apparent contour\n", emb->k, emb->n, emb->choice);

  freedualembedding (dual);
  postprocesssketch (sketch);
  return (sketch);
}

/*
 * elongated regions containing triple nodes
 */

int
print_sketch_ra3 (int only3, int regionnum, struct embedding *emb, struct sketch *s, int *nodemark)
{
  int vi, vin, vinn, vj, vjn;
  int found, val, arcid;
  struct emb_node *node, *nextnode, *nextnextnode;
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *arc;

  /*
   * now the elongated regions containint the triple-points
   * (one or two)
   */
  found = 0;
  for (vi = 0; vi < emb->k + emb->n; vi++)
  {
    node = &(emb->nodes[vi]);
    if (node->valency != 4) continue;
    for (vj = 0; vj < 4; vj++)
    {
      if ((nodemark[vi] & (1<<vj)) != 0) continue;
      vin = node->ping[vj];
      nextnode = &(emb->nodes[vin]);
      if (only3 && nextnode->valency == 4) continue;
      found++;
      break;
    }
    if (found) break;
  }

  if (! found) return (0);

  region = newregion_tail (s);
  assert (region->border == 0);
  assert (region->tag == regionnum);
  region->f = 2;
  bl = newborderlist (region);
  blast = 0;
  while ((nodemark[vi] & (1<<vj)) == 0)
  {
    nodemark[vi] |= 1 << vj;

    vin = node->ping[vj];
    nextnode = &(emb->nodes[vin]);
    vjn = node->pong[vj];
    val = nextnode->valency;
    vjn = (vjn + 1) % val;
    while (nextnode->valency <= 3)
    {
      nodemark[vin] |= 1 << vjn;
      vinn = nextnode->ping[vjn];
      nextnextnode = &(emb->nodes[vinn]);
      vjn = nextnode->pong[vjn];  // make a step forward for vjn
      val = nextnextnode->valency;
      vjn = (vjn + 1) % val;
      vin = vinn;                     // make a step forward for vin and nextnode
      nextnode = nextnextnode;
    }

    val = nextnode->valency;
    assert (val == 4);
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b; 
    b->orientation = 1;

    arcid = 4*(vi - emb->k) + vj + 1;
    for (arc = s->arcs; arc; arc = arc->next)
    {
      if (arc->tag == arcid)
      {
        b->info = arc;
        break;
      }
    }

    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;  
    b->orientation = -1;

    arcid = 4*(vin - emb->k) + vjn + 1 + 4*emb->n;
    for (arc = s->arcs; arc; arc = arc->next)
    {
      if (arc->tag == arcid)
      {
        b->info = arc;
        break;
      }
    }

    vjn = (vjn + 3) % 4;
    vi = vin;                      // make a step forward for node, vi, vj
    node = nextnode;
    vj = vjn;
  }

  return (1);
}

/*
 *
 */

void
printembedding (struct embedding *emb)
{
  struct emb_node *node;
  int i, k, s, ii, kk;

  start_comment ();
  printf ("Warning: noncanonical\n");

  printf ("embedding:%d {", emb->choice);
  for (i = 0; i < emb->k + emb->n; i++)
  {
    if (i > 0) printf (",");
    printf (" %d: (", i);
    node = &emb->nodes[i];
    for (s = 0; s < node->valency; s++)
    {
      if (s > 0) printf (", ");
      printf ("%d", node->ping[s]);
    }
    printf (")");
  }
  printf (" }\n");

  if (verbose && emb->k > 0)
  {
    if (emb->k == 2)
    {
      if (emb->connections[0]/3 +
          emb->connections[1]/3 +
          emb->connections[2]/3 == 3)
        printf ("type: theta-curve\n");
       else
        printf ("type: handcuff\n");
    } else {
      printf ("Strands connection matrix:\n");
      for (i = 0; i < emb->k; i++)
      {
        for (k = 0; k < 3; k++)
        {
          ii = emb->connections[3*i + k]/3;
          kk = emb->connections[3*i + k] - 3*ii;
          // printf ("ii = %d kk = %d\n", ii, kk);
          printf ("%d.%d ", ii, kk);
        }
        printf ("\n");
      }
    }
  }
  return;
}

/*
 *
 */

void
freedualembedding (struct dualembedding *dual)
{
  assert (dual->regions);
  free (dual->wedgeij);
  freedualregions (dual->regions);
  free (dual);
  return;
}

void
freedualregions (struct dual_region *region)
{
  if (region == 0) return;
  //id = region->id;
  if (region->next) freedualregions (region->next);
  assert (region->ping);
  free (region->ping);
  // free (region->wedgeij);  // this is a portion of a large contiguous vector
  free (region);
  return;
}

/*
 *
 *
 */

struct dualembedding *
embedding2dual (struct embedding *emb)
{
  int i, j, jj, val, nodenum, numwedges, count, totcount, found;
  int regionsize, region_id;
  int inode, iinode, jwedgeminus, jwedgeplus;
  int *wedgemark, *r_wedgeij_all, *r_wedgeij;
  struct emb_node *node, *thisnode, *nextnode;
  int thisnodej, nextnodei, nextnodej;
  struct dualembedding *dual;
  struct dual_region *region, *prevregion, *adjregion;

  dual = (struct dualembedding *) malloc (sizeof (struct dualembedding));

  assert (emb->k % 2 == 0);
  nodenum = dual->v = emb->n + emb->k;
  dual->e = 2*emb->n + 3*emb->k/2;
  numwedges = dual->numwedges = 4*emb->n + 3*emb->k;
  dual->numregions = dual->e - dual->v + 2;
  dual->regions = 0;
  wedgemark = (int *) malloc (4*nodenum*sizeof(int));

  for (i = 0; i < nodenum; i++)
  {
    node = &emb->nodes[i];
    for (j = 0; j < node->valency; j++)
    {
      wedgemark[4*i+j] = 0;  /* will be >0 for wedges of build regions */
    }
  }

  /*
   * now loop through all "wedges" (pair of consecutive arcs for each node
   * a wedge is indexed using the first of the two arcs (counterclockwise ordering)
   */

  totcount = 0;
  region_id = 0;
  r_wedgeij_all = (int *) malloc (numwedges * sizeof(int));
  dual->wedgeij = r_wedgeij_all;
  prevregion = 0;
  for (i = 0; i < nodenum; i++)
  {
    node = &emb->nodes[i];
    for (j = 0; j < node->valency; j++)
    {
      if (wedgemark[4*i + j]) continue;
      /* wedge of new region! */
      r_wedgeij = &r_wedgeij_all[totcount];
      /*
       * walk around the region
       */
      count = 0;
      thisnode = node;
      thisnodej = j;
      r_wedgeij[count] = 4*i + j;
      while ((nextnodei = thisnode->ping[thisnodej]) != i)
      {
        count++;
        nextnode = &emb->nodes[nextnodei];
        val = nextnode->valency;
        nextnodej = thisnode->pong[thisnodej];
        nextnodej = (nextnodej + val - 1) % val;

        wedgemark[4*nextnodei + nextnodej] = 1;
        r_wedgeij[count] = 4*nextnodei + nextnodej;
        thisnode = nextnode;
        thisnodej = nextnodej;
      }
      count++;

      regionsize = sizeof (struct dual_region) + count*sizeof (int);
      region = (struct dual_region *) malloc (regionsize);
      region->valency = count;
      region->id = region_id;
      region->wedgeij = r_wedgeij;
      region->ping = (struct dual_region **) malloc (region->valency * sizeof (struct dual_region *));
      region->next = 0;
      if (prevregion)
      {
        prevregion->next = region;
      } else {
        dual->regions = region;
      }
      //dual->regions = region;
      totcount += count;
      region_id++;
      prevregion = region;
      assert (region_id <= dual->numregions);
    }
  }
  assert (totcount == numwedges);

  /*
   * fill 'ping' vector (pointers to adjregion)
   */

  for (region = dual->regions; region; i++, region = region->next)
  {
    for (j = 0; j < region->valency; j++)
    {
      inode = region->wedgeij[j]/4;
      node = &emb->nodes[inode];
      jwedgeminus = region->wedgeij[j] % 4;
      jwedgeplus = (jwedgeminus + 1) % node->valency;
      /* search for the neighbouring region */
      found = 0;
      for (adjregion = dual->regions; adjregion; adjregion = adjregion->next)
      {
        if (adjregion == region) continue;
        for (jj = 0; jj < adjregion->valency; jj++)
        {
          iinode = adjregion->wedgeij[jj]/4;
          if (iinode != inode) continue;
          if (adjregion->wedgeij[jj]%4 != jwedgeplus) continue;
          found = 1;
          break;
        }
        if (found) break;
      }
      assert (adjregion);
      region->ping[j] = adjregion;
    }
  }

  free (wedgemark);
  return (dual);
}

void
printdualembedding (struct dualembedding *dual, struct embedding *emb)
{
  struct dual_region *region, *adjregion;
  int jj;

  if (verbose) printf ("Dual of given embedding (%d regions):\n", dual->numregions);

  printf ("dualembedding {");
  for (region = dual->regions; region; region = region->next)
  {
    if (region != dual->regions) printf (", ");
    printf ("%d:(", region->id);
    for (jj = 0; jj < region->valency; jj++)
    {
      adjregion = region->ping[jj];
      if (jj > 0) printf (", ");
      printf ("%d", adjregion->id);
      if (verbose) printf ("<%d.%d>", region->wedgeij[jj]/4, region->wedgeij[jj]%4);
    }
    printf (")");
  }
  printf ("}\n");

  if (verbose) print_dual_type (dual, emb);
}

#define MAXTRI 1024

void  
print_dual_type (struct dualembedding *dual, struct embedding *emb)
{   
  int i,j, in, inode, k;
  int numregions;
  struct dual_region *region;
  struct emb_node *node;
  int *vec;
  int counttri, largest, largestj;

  printf ("dual_type: ");

  numregions = dual->numregions;
  vec = (int *) malloc (dual->numregions * sizeof(int));

  for (region = dual->regions, i = 0; region; region = region->next, i++)
  {
    vec[i] = MAXTRI*region->valency;
    counttri = 0;
    for (in = 0; in < region->valency; in++)
    {
      inode = region->wedgeij[in]/4;
      node = &emb->nodes[inode];
      if (node->valency == 3) counttri++;
    }
    assert (counttri < MAXTRI);
    vec[i] += counttri;
  }

  //dosortvec (sortvec, regionsnum);

  for (i = 0; i < numregions; i++)
  {
    if (i > 0) printf (", ");
    largest = -1;
    for (j = 0; j < numregions; j++)
    {
      if (vec[j] > largest)
      {
        largest = vec[j];
        largestj = j;
      }
    }
    counttri = largest % MAXTRI;
    printf ("%d", (largest-counttri)/MAXTRI);
    for (k = 0; k < counttri; k++) printf ("*");
    vec[largestj] = -1;
  }
  printf ("\n");
  free (vec);
}

/*
 * at first we have one generator per short arc.  Numbered as follows
 * 1. Those exiting from a trivalent node are numbered from 0 to open_arcs, in the order
 *    we find them
 * 2. Each crossing has two exiting arcs: even and odd; numbering is
 *    2*(id-k) + open_arcs       [even]
 *    2*(id-k) + open_arcs + 1   [odd]
 *
 * Relators:
 * at a trivalent node: abc if all arcs are starting here.
 *                      otherwise the respective inverse
 * at an overpass: g_in = g_out
 * at an underpass: bdAC or dbCA where a,b are the generators on the overpass, c,d the generators
 *                  of the underpass; the first relator if the overpass sees the underpass
 *                  oriented from right to left
 */

#define NODE_IS_START 1
#define NODE_IS_ARRIVAL 2

struct presentation *
wirtingerfromembedding (struct embedding *emb)
{
  struct presentation *p;
  struct vecofintlist *gaussloiv;
  struct vecofintlist *dtloiv;

  if (emb->k == 0)
  {
    if (debug) printf ("This is a standard knot/link, converting to gausscode\n");
    gaussloiv = embeddingtoloiv (emb);
    if (debug) printloiv (gaussloiv);
    freeembedding (emb);
    dtloiv = gausscode2dtcode (gaussloiv);
    if (debug) printloiv (dtloiv);
    freeloiv (gaussloiv);
    p = wirtingerfromloiv (dtloiv);
    if (debug) print_presentation (p);
    freeloiv (dtloiv);
    return (p);
  }

  p = wirtingerfromembeddingraw (emb);

  assert (emb->numhcomponents >= 1);
  assert (globals.ccemb != 0);
  if (emb->numhcomponents + emb->numrings == 1)
  {
    assert (emb->numrings == 0);
    emb_meridians_longitudes (emb, p, -1);
  } else {
    if (globals.ccemb > 0)
    {
      emb_meridians_longitudes (emb, p, globals.ccemb);
    } else if (!quiet) printf ("Cannot compute meridians and longitudes for links, consider using option --ccemb\n");
  }

  if (debug) print_presentation (p);
  emb_remove_dup_rules (p);

  if (globals.simplifypresentation) simplify_presentation (p);
  return (p);
}

/*
 * interpret a 'ring' component as an element of the fundamental group, computed as
 * a Wirtinger presentation with some components excluded (see excludecc)
 *
 * if globals.loopcc is zero use a default value if only one ring is present
 *
 * the element is added as selected element
 */

struct presentation *
ccasloop (struct embedding *emb)
{
  struct emb_node *node;
  int i, j, k, starti, startj, inext, jback, jorto;
  int loop;
  struct presentation *p;
  struct presentationrule *rule;

  emb_color4 (emb);

  if (emb->numrings <= 0)
  {
    fprintf (stderr, "No loop components present in the embedding\n");
    return (0);
  }

  if (emb->numhcomponents+emb->numrings <= 1)
  {
    fprintf (stderr, "Connected objects are not allowed\n");
    return (0);
  }

  if (globals.loopcc < 0)
  {
    if (quiet)
    {
      printf ("%d\n", emb->numhcomponents);
      printf ("%d\n", emb->numrings);
    } else {
      printf ("There are %d components with genus larger than one\n", emb->numhcomponents);
      printf ("There are %d components of genus one\n", emb->numrings);
    }
    return (0);
  }

  if (globals.loopcc == 0)
  {
    if (emb->numrings > 1)
    {
      fprintf (stderr, "More than one loop component is present, you must indicate which one to use with option\n");
      fprintf (stderr, "   --loopcc <n>\n");
      fprintf (stderr, "Loop components are numbered after the higher-genus components\n");
      return (0);
    }
    globals.loopcc = emb->numhcomponents + 1;
  }

  loop = globals.loopcc;
  if (loop <= emb->numhcomponents) {fprintf (stderr, "indicated component (%d) is not a loop\n", loop); return (0);}
  if (loop > emb->numhcomponents + emb->numrings) {fprintf (stderr, "indicated component (%d) does not exist\n", loop); return (0);}

  if (debug)
  {
    printf ("numhcomponents: %d, numrings: %d\n", emb->numhcomponents, emb->numrings);
    printf ("k: %d, n: %d\n", emb->k, emb->n);
    for (i = 0; i < emb->k + emb->n; i++)
    {
      node = &emb->nodes[i];
      printf ("node %d, valency: %d, color: %d, colorodd: %d\n", i, node->valency, node->color, node->colorodd);
    }
  }

  p = wirtingerfromembeddingraw (emb);

  //assert (emb->numhcomponents >= 1);
  if (emb->numhcomponents + emb->numrings == 1)
  {
    printf ("More than one component required\n");
    return (0);
  }

  /*
   * now remove selected components
   */

  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    assert (node->valency == 3);
    if (isexcluded (node->color, loop))
    {
      for (j = 0; j < 3; j++)
      {
        rule = (struct presentationrule *) malloc (1*sizeof (int) + sizeof (struct presentationrule));
        rule->length = 1;
        rule->next = p->rules;
        p->rules = rule;
        rule->var[0] = node->generator[j] + 1;
      }
    }
  }

  for (i = emb->k; i < emb->k + emb->n; i++)
  {
    node = &emb->nodes[i];
    if (isexcluded (node->color, loop))
    {
      rule = (struct presentationrule *) malloc (1*sizeof (int) + sizeof (struct presentationrule));
      rule->length = 1;
      rule->next = p->rules;
      p->rules = rule;
      rule->var[0] = node->generator[0] + 1;
    }
    if (isexcluded (node->colorodd, loop))
    {
      rule = (struct presentationrule *) malloc (1*sizeof (int) + sizeof (struct presentationrule));
      rule->length = 1;
      rule->next = p->rules;
      p->rules = rule;
      rule->var[0] = node->generator[1] + 1;
    }
  }

  /* build selected element(s) */
  for (i = emb->k; i < emb->k + emb->n; i++)
  {
    node = &emb->nodes[i];
    if (node->color == loop)
    {
      j = 0;
      break;
    }
    if (node->colorodd == loop)
    {
      j = 1;
      break;
    }
  }

  if (node->direction[j] == NODE_IS_ARRIVAL) j += 2;
  starti = i;
  startj = j;
  rule = (struct presentationrule *) malloc (2*emb->n*sizeof (int) + sizeof (struct presentationrule));
  k = 0;
  while (1)
  {
    if (debug) printf ("Along selected: node %d, direction %d.", i, j);
    if (((j + node->overpassisodd) % 2) == 1)
    {
      if (debug) printf (" Is underpass\n");
      jorto = (j + 1) % 4;
      if (node->direction[jorto] == NODE_IS_ARRIVAL)
      {
        rule->var[k++] = node->generator[jorto] + 1;
        if (debug) printf ("generator +%d\n", node->generator[jorto]);
      } else {
        rule->var[k++] = -(node->generator[jorto] + 1);
        if (debug) printf ("generator -%d\n", node->generator[jorto]);
      }
    } else {
      if (debug) printf (" Is overpass\n");
    }

    inext = node->ping[j];
    jback = node->pong[j];

    i = inext;
    node = &emb->nodes[i];
    j = (jback + 2) % 4;

    if (i == starti && j == startj) break;
  }

  rule->length = k;
  if (debug) printf ("Length of selected element: %d\n", k);
  if (!globals.loopasrelator)
  {
    rule->next = p->elements;
    p->elements = rule;
  }
  if (globals.loopasrelator)
  {
    rule->next = p->rules;
    p->rules = rule;
  }

  if (debug) print_presentation (p);
  emb_remove_dup_rules (p);
  if (globals.simplifypresentation) simplify_presentation (p);
  return (p);
}

/*
 * check if component with given color must be excluded
 */

int
isexcluded (int color, int loop)
{
  int k;

  if (color == loop) return (1);
  for (k = 0; k < globals.numexcluded; k++)
  {
    if (color == globals.exclude[k]) return (1);
  }
  return (0);
}

/*
 * compute raw wirtinger presentation (one generator per short arc)
 */

char
lgen (int gen)
{
  if (gen > 0) return ('a' + gen - 1);
  if (gen < 0) return ('A' - gen - 1);
  return ('?');
}

struct presentation *
wirtingerfromembeddingraw (struct embedding *emb)
{
  int open_arcs, short_arcs;
  struct emb_node *node;
  int i, j;
  struct presentation *p;
  struct presentationrule *rule;
  int sign, a, b, c, d;

  assert ((emb->k % 2) == 0);
  open_arcs = emb->k/2*3;
  short_arcs = open_arcs + emb->n*2;

  assert (emb->orientation);

  p = (struct presentation *) malloc (sizeof (struct presentation));
  p->gennum = short_arcs;
  p->elements = 0;
  p->rules = 0;
  p->characteristic = 0;
  p->espected_deficiency = 1;

  if (debug)
  {
    for (i = 0; i < emb->k + emb->n; i++)
    {
      node = &emb->nodes[i];
      printf ("generators at node %d: ", i);
      for (j = 0; j < node->valency; j++)
      {
        printf (" %c", (node->direction[j]==NODE_IS_START)?lgen(node->generator[j]+1):lgen(-node->generator[j]-1));
      }
      printf ("\n");
    }
  }

  for (i = emb->k; i < emb->k + emb->n; i++)
  {
    node = &emb->nodes[i];
    assert (node->valency == 4);
    rule = (struct presentationrule *) malloc (2*sizeof (int) + sizeof (struct presentationrule));
    rule->length = 2;
    rule->next = p->rules;
    p->rules = rule;
    sign = 1;  // positive crossing
    if (node->overpassisodd)
    {
      rule->var[0] = node->generator[1] + 1;
      rule->var[1] = - (node->generator[3] + 1);
      if (node->direction[1] == NODE_IS_ARRIVAL)
      {
        a = node->generator[1];
        b = node->generator[3];
      } else {
        sign *= -1;
        a = node->generator[3];
        b = node->generator[1];
      }
      if (node->direction[0] == NODE_IS_ARRIVAL)
      {
        sign *= -1;
        c = node->generator[0];
        d = node->generator[2];
      } else {
        c = node->generator[2];
        d = node->generator[0];
      }
    } else {
      rule->var[0] = node->generator[0] + 1;
      rule->var[1] = - (node->generator[2] + 1);
      if (node->direction[0] == NODE_IS_ARRIVAL)
      {
        a = node->generator[0];
        b = node->generator[2];
      } else {
        sign *= -1;
        a = node->generator[2];
        b = node->generator[0];
      }
      if (node->direction[1] == NODE_IS_ARRIVAL)
      {
        c = node->generator[1];
        d = node->generator[3];
      } else {
        sign *= -1;
        c = node->generator[3];
        d = node->generator[1];
      }
    }
    rule = (struct presentationrule *) malloc (4*sizeof (int) + sizeof (struct presentationrule));
    rule->length = 4;
    rule->next = p->rules;
    p->rules = rule;
    if (sign > 0)
    {
      rule->var[0] = b + 1;
      rule->var[1] = d + 1;
      rule->var[2] = -(a + 1);
      rule->var[3] = -(c + 1);
    } else {
      rule->var[0] = d + 1;
      rule->var[1] = b + 1;
      rule->var[2] = -(c + 1);
      rule->var[3] = -(a + 1);
    }
  }

  if (verbose) printf ("RELATORS AT NODES:\n");
  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    assert (node->valency == 3);
    rule = (struct presentationrule *) malloc (3*sizeof (int) + sizeof (struct presentationrule));
    rule->length = 3;
    rule->next = p->rules;
    p->rules = rule;

    if (verbose) printf ("NODE RELATOR ");
    for (j = 0; j < 3; j++)
    {
      if (verbose) printf (" %c%d", (node->direction[j]==NODE_IS_ARRIVAL)?'-':'+', node->generator[j]);
      sign = 1;
      if (node->direction[j]==NODE_IS_ARRIVAL) sign = -1;
      rule->var[j] = sign*(node->generator[j]+1);
    }
    if (verbose) printf ("\n");
  }

  return (p);
}

/*
 * count crossings by color type
 */

void
printcrossingcolors (struct embedding *emb)
{
  int ccnum;
  int i, j, rowcol, colrow;
  int *matrix;
  struct emb_node *node;

  emb_color (emb);
  emb_color4 (emb);

  ccnum = emb->numhcomponents + emb->numrings;
  matrix = (int *) malloc (ccnum*ccnum*sizeof(int));
  for (i = 0; i < ccnum*ccnum; i++) matrix[i] = 0;

  for (i = emb->k; i < emb->k+emb->n; i++)
  {
    node = &emb->nodes[i];
    assert (node->valency == 4);

    rowcol = ccnum*(node->color-1) + node->colorodd - 1;
    assert (rowcol >= 0 && rowcol < ccnum*ccnum);

    matrix[rowcol]++;
  }

  for (i = 0; i < ccnum; i++)
  {
    rowcol = i*ccnum + i;
    if (quiet) printf ("%d-%d: %d\n", i+1, i+1, matrix[rowcol]);
     else printf ("number of selfcrossings of component %d: %d\n", i+1, matrix[rowcol]);
  }

  for (i = 0; i < ccnum; i++)
  {
    for (j = i+1; j < ccnum; j++)
    {
      rowcol = i*ccnum + j;
      colrow = j*ccnum + i;
      if (quiet) printf ("%d-%d: %d\n", i+1, j+1, matrix[rowcol]+matrix[colrow]);
       else printf ("number of mutual crossings between components %d and %d: %d\n", i+1, j+1, matrix[rowcol]+matrix[colrow]);
    }
  }
}

/*
 * compute arc-connectedness of the embedding
 */

int is2connected (struct dualembedding *dual);
int is3connected (struct dualembedding *dual);
int pingpong (struct dual_region *r, int ii);
int dual_left_turn (struct dual_region *r, int ii1, int ii2);

int
embedding_connectedness (struct dualembedding *dual, struct embedding *emb)
{
  if (is2connected (dual)) return (2);
  if (is3connected (dual)) return (3);
  
  return (9999);
}

int
is2connected (struct dualembedding *dual)
{
  int i, j;
  struct dual_region *region;

  if (dual == 0) return (0);

  for (region = dual->regions; region; region = region->next)
  {
    for (i = 0; i < region->valency; i++)
    {
      for (j = i+1; j < region->valency; j++)
      {
        if (region->ping[i]->id == region->ping[j]->id) return (region->id + 1);
      }
    }
  }
  return (0);
}

int
is3connected (struct dualembedding *dual)
{
  int i, ii, jj, kk, ripp, rjpp, rkpp;
  struct dual_region *ri, *rj, *rk, *rfinal;

  if (dual == 0) return (0);

  for (i = 0, ri = dual->regions; ri; i++, ri = ri->next)
  {
    assert (i == ri->id);
    for (ii = 0; ii < ri->valency; ii++)
    {
      ripp = pingpong (ri, ii);
      rj = ri->ping[ii];
      if (rj <= ri) continue;

      for (jj = 0; jj < rj->valency; jj++)
      {
        if (jj == ripp) continue; // do not backtrack
        rjpp = pingpong (rj, jj);
        rk = rj->ping[jj];
        if (rk == ri) continue;
        if (rk <= ri) continue;
        for (kk = 0; kk < rk->valency; kk++)
        {
          if (kk == rjpp) continue; //do not backtrack
          rkpp = pingpong (rk, kk);
          rfinal = rk->ping[kk];
          if (rfinal == ri)
          {
            // printf ("found: %d.%d.%d  %d.%d.%d  %d.%d.%d\n", ri->id, ii, pingpong(ri,ii), rj->id, jj, pingpong(rj,jj), rk->id, kk, pingpong(rk,kk));
            if (dual_left_turn(ri, rkpp, ii) &&
                dual_left_turn(rj, ripp, jj) &&
                dual_left_turn(rk, rjpp, kk) ) continue;
            if (dual_left_turn(ri, ii, rkpp) &&
                dual_left_turn(rj, jj, ripp) &&
                dual_left_turn(rk, kk, rjpp) ) continue;
            return (i + 1);
          }
        }
      }
    }
  }

  return (0);
}

int
pingpong (struct dual_region *r, int ii)
{
  struct dual_region *rnext;
  int jj;

  assert (ii < r->valency);
  rnext = r->ping[ii];
  for (jj = 0; jj < rnext->valency; jj++)
  {
    if (rnext->ping[jj] == r) return (jj);
  }
  return (9999);
}

/*
 * check if i1 is i2+1 mod valency
 */

int
dual_left_turn (struct dual_region *r, int i1, int i2)
{
  int i2plus;

  i2plus = (i2+1) % r->valency;
  if (i2plus == i1) return (1);
  return (0);
}

/*
 *
 */

/*
 * Just read the embedding (possibly comprehensive of "choice" information)
 * accepted syntax is:
 *
 *   embedding[:<choice>] { <node>, <node>, ... };
 *
 * where node has syntax:
 *   node ::= <id>: (<dest0>, <dest1>, <dest2>)
 * or
 *   node ::= <id>: (<dest0>, <dest1>, <dest2>, <dest3>)
 * based on whether this is a trivalent node or a crossing
 *
 * <choice> is zero if not present; it has to be interpreted as a binary number
 * where the least significant digit refers to the first crossing with a value of
 * 0 means that the overcrossing connects even entries (first and second in '(e, o, e, o)', index starting at 0)
 * 1 means that the overcrossing connects odd entries (second and fourth)
 * choice can alternatively be specified using the 'contour' option --choice
 * an hexadecimal value (syntax 0xnn) can be given
 *
 * As an example the following input:
 *
 *   embedding:0x6 { 0: (3, 6, 7), 1: (2, 7, 6), 2: (1, 5, 4), 3: (0, 4, 5), 4: (2, 3, 7, 7), 5: (2, 6, 6, 3), 6: (0, 5, 5, 1), 7: (0, 1, 4, 4) };
 *
 * describes an embedding with 4 trivalent nodes, 4 crossings and a 'choice' value of 0x6 (hexadecimal 6) meaning that
 * for the second and third crossing the overpass connects refers to the nodes in odd position in "<id>: (even, odd, even, odd)
 * [position index starts from zero].
 *
 * it is assumed that inside parentheses the first number (id of arrival node) appears first.
 * it is assumed that trivalent nodes are numbered first
 */

struct embedding *
readembedding (FILE *file)
{
  struct embedding *emb;
  struct emb_node *node, *nodeend, *prevnode, *nodesvec;
  int i, j, tok, count, choice;
  int iend, iendprev, jend;

  //int node_id, node_id2, tnode_id, tnode_pt;

  emb = (struct embedding *) malloc (sizeof (struct embedding));
  emb->k = emb->n = 0;
  emb->choice = -1;
  emb->orientation = 0; //this is controlled by a sign or by commandline argument
  emb->connections = 0;
  emb->nodes = 0;

  tok = gettoken (file);
  if (tok == TOK_COLON)
  {
    tok = gettoken (file);
    assert (tok == ISNUMBER);
    emb->choice = gettokennumber ();
  } else ungettoken (tok);

  if (emb->choice >= 0 && globals.choice >= 0 && !quiet)
  {
    printf ("Value of 'choice' given via '--choice 0x%x' option takes precedence on the value 0x%x indicated in", globals.choice, emb->choice);
    printf (" the embedding:0x%x description\n", emb->choice);
  }
  if (globals.choice >= 0) emb->choice = globals.choice;
  if (emb->choice < 0)
  {
    if (!quiet) printf ("'choice' value not given, assuming zero\n");
    emb->choice = 0;
  }

  tok = gettoken (file);
  if (tok == TOK_MINUS || tok == TOK_PLUS)
  {
    if (tok == TOK_MINUS) emb->orientation = -1;
    if (tok == TOK_PLUS) emb->orientation = 1;
    tok = gettoken (file);
  }
  if (globals.rotation)
  {
    if (emb->orientation && !quiet) printf ("Command line option forces orientation of embedding to be %s\n",
       (globals.rotation > 0)?"counterclockwise":"clockwise");
    emb->orientation = (globals.rotation > 0)?1:(-1);
  }
  if (emb->orientation == 0)
  {
    if (verbose) printf ("Default orientation for embedding is counterclockwise\n");
    emb->orientation = 1;
  }
  assert (tok == TOK_LBRACE);
  node = 0;
  while (1)
  {
    tok = gettoken (file);
    if (tok == TOK_RBRACE) break;
    assert (tok == ISNUMBER);
    prevnode = node;
    node = (struct emb_node *) malloc (sizeof (struct emb_node));
    node->id = gettokennumber ();
    node->valency = 0;
    node->overpassisodd = -1;  // undefined
    node->next = 0;
    tok = gettoken (file);
    assert (tok == TOK_COLON);
    assert (gettoken (file) == TOK_LPAREN);

    for (i = 0; i < 4; i++)
    {
      assert (gettoken (file) == ISNUMBER);
      node->ping[i] = gettokennumber ();
      node->valency++;
      tok = gettoken (file);
      if (tok == TOK_RPAREN) break;
      assert (tok == TOK_COMMA);
    }

    assert (node->valency == 3 || node->valency == 4);
    if (node->valency == 3) emb->k++;
    if (node->valency == 4) emb->n++;
    if (prevnode) prevnode->next = node;
      else emb->nodes = node;

    tok = gettoken (file);
    if (tok == TOK_RBRACE) break;
    assert (tok == TOK_COMMA);
  }

  /* move each node data into a single vector for easier referencing */
  nodesvec = (struct emb_node *) malloc ((emb->k+emb->n)*sizeof (struct emb_node));
  for (node = emb->nodes, count = 0; node; node = node->next, count++)
  {
    assert (node->id == count);
    memcpy (nodesvec+count, node, sizeof (struct emb_node));
    nodesvec[count].next = node;
  }
  for (i = 0; i < emb->k + emb->n; i++)
  {
    node = nodesvec+i;
    free (node->next);
    node->next = 0;
  }

  /* postprocess: add information such as "pong" and overpasses */
  /* sanity check: we require id of crossings to be >= emb->k */

  choice = emb->choice;
  for (i = 0; i < emb->k + emb->n; i++)
  {
    node = nodesvec + i;
    iendprev = -1;
    for (j = 0; j < node->valency; j++)
    {
      iend = node->ping[j];
      if (iend == iendprev)
      {
         node->pong[j] = node->pong[j-1] - 1;
      } else {
        nodeend = nodesvec + iend;
        node->pong[j] = -1;
        for (jend = nodeend->valency - 1; jend >= 0; jend--)
        {
          if (nodeend->ping[jend] != i) continue;
          node->pong[j] = jend;
          break;
        }
        assert (node->pong[j] >= 0);
      }
      iendprev = iend;
      //printf ("computing pong: from node (%d,%d) -> (%d,%d)\n", i, j, iend, node->pong[j]);
    }
    if (node->valency == 4)
    {
      node->overpassisodd = 0;
      if ((choice & 1) != 0) node->overpassisodd = 1;
      choice /= 2;
    }
  }
  emb->nodes = nodesvec;

  emb_orient (emb);

  return (emb);
}

/*
 * color the trivalent nodes according to connected components
 * return with the number of colors
 */

int
emb_color (struct embedding *emb)
{
  int i, ii, iito, jj, color, colored;
  int goon;
  struct emb_node *node, *nodefrom, *nodeto;

  if (emb->k < 2) return (0);

  /* first reset all colors */
  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    node->color = 0;
  }

  color = 0;
  colored = 0;

  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    if (node->color) continue;
    color++;
    colored++;
    node->color = color;
    goon = 1;
    while (goon) /* expand colors */
    {
      goon = 0;
      for (ii = 0; ii < emb->k; ii++)
      {
        nodefrom = &emb->nodes[ii];
        if (nodefrom->color == 0) continue;
        for (jj = 0; jj < 3; jj++)
        {
          iito = emb->connections[3*ii + jj]/3;
          nodeto = &emb->nodes[iito];
          if (nodeto->color != 0) {assert (nodeto->color == nodefrom->color); continue;}
          if (debug) printf ("extending color %d from node %d.%d to node %d\n", nodefrom->color, ii, jj, iito);
          nodeto->color = nodefrom->color;
          goon++;
          colored++;
        }
      }
    }
  }
  assert (colored == emb->k);
  return (color);
}

/*
 * color the four-valent nodes
 */

/* spread coloring function */

int expand_from_4 (struct embedding *emb);

int
expand_from_4 (struct embedding *emb)
{
  int i, jj, jjto;
  int goon;
  struct emb_node *node, *nodeto;

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (i = emb->k; i < emb->k+emb->n; i++)
    {
      node = &emb->nodes[i];
      assert (node->valency == 4);
      if (node->color)
      {
        /* even is colored, try to expand */
        for (jj = 0; jj < 4; jj += 2)
        {
          jjto = node->ping[jj];
          nodeto = &emb->nodes[jjto];
          if ((node->pong[jj] % 2) == 0)
          {
            if (nodeto->color == 0) goon++;
             else assert (nodeto->color == node->color);
            nodeto->color = node->color;
          } else {
            if (nodeto->colorodd == 0) goon++;
             else assert (nodeto->colorodd == node->color);
            nodeto->colorodd = node->color;
          }
        }
      }
      if (node->colorodd)
      {
        /* odd is colored, try to expand */
        for (jj = 1; jj < 4; jj += 2)
        {
          jjto = node->ping[jj];
          nodeto = &emb->nodes[jjto];
          if ((node->pong[jj] % 2) == 0)
          {
            if (nodeto->color == 0) goon++;
             else assert (nodeto->color == node->colorodd);
            nodeto->color = node->colorodd;
          } else {
            if (nodeto->colorodd == 0) goon++;
             else assert (nodeto->colorodd == node->colorodd);
            nodeto->colorodd = node->colorodd;
          }
        }
      }
    }
  }

  return (1);
}

int
emb_color4 (struct embedding *emb)
{
  int i, jj, jjto, istart;
  int iring;
  struct emb_node *node, *nodeto;

  /* first reset all colors */
  for (i = emb->k; i < emb->k+emb->n; i++)
  {
    node = &emb->nodes[i];
    node->color = node->colorodd = 0;
  }

  /* first: expand from trivalent nodes */
  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    assert (node->valency == 3);
    for (jj = 0; jj < 3; jj++)
    {
      jjto = node->ping[jj];
      nodeto = &emb->nodes[jjto];
      if (nodeto->valency != 4) continue;
      if ((node->pong[jj] % 2) == 0) nodeto->color = node->color;
       else nodeto->colorodd = node->color;
    }
  }

  /* second: expand from four-valent nodes */
  expand_from_4 (emb);

  /* now color rings */

  for (iring = 0; iring < emb->numrings; iring++)
  {
    /* find an uncolored node */
    istart = -1;
    for (i = emb->k; i < emb->k+emb->n; i++)
    {
      node = &emb->nodes[i];
      assert (node->valency == 4);
      if (node->color == 0)
      {
        istart = i;
        node->color = emb->numhcomponents + iring + 1;
        break;
      }
      if (node->colorodd == 0)
      {
        istart = i;
        node->colorodd = emb->numhcomponents + iring + 1;
        break;
      }
    }
    assert (istart >= 0);
    expand_from_4 (emb);
  }

  /* check everything is colored */

  for (i = emb->k; i < emb->k+emb->n; i++)
  {
    node = &emb->nodes[i];
    assert (node->color);
    assert (node->colorodd);
  }

  return (emb->numrings);
}

int
emb_meridians_longitudes (struct embedding *emb, struct presentation *p, int ccemb)
{
  int i, ii, iii, kk, iinext, kknext;
  int ikparent;
  int k, goon;
  int *node_flood, *arcs, *underpasses;
  struct emb_node *node, *node2;
  struct presentationrule *rule;
  int llength, llengthpre, llengthpost;
  int u, count;

  if (ccemb >= 0)
  {
    assert (ccemb > 0);
    if (ccemb > emb->numhcomponents)
    {
      return (emb_meridians_longitudes_torus (emb, p, ccemb));
    }
  }

  /*
   * step 1. construct a spanning tree for the spatial graph
   * rooted at the first node with the ccemb color (or 1 if ccemb < 0)
   */

  assert (emb->k > 0);

  node_flood = (int *) malloc (emb->k * sizeof(int));
  arcs = (int *) malloc (3*emb->k * sizeof(int));
  underpasses = (int *) malloc (3*emb->k * sizeof(int));

  for (i = 0; i < emb->k; i++)
  {
    node_flood[i] = -1;
    for (k = 0; k < 3; k++) arcs[3*i+k] = 0;  // mark as NOT spanning
  }

  if (ccemb > 0)
  {
    for (i = 0; i < emb->k; i++)
    {
      node = &emb->nodes[i];
      if (node->color == ccemb)
      {
        node_flood[i] = 0;
        //node_flood[i] = emb->connections[3*i + 0];
        break;
      }
    }
  } else {
    node_flood[0] = 0;
  }

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (i = 0; i < emb->k; i++)
    {
      if (node_flood[i] >= 0) continue;
      for (k = 0; k < 3; k++)
      {
        ii = emb->connections[3*i + k]/3;
        kk = emb->connections[3*i + k] - 3*ii;
        if (node_flood[ii] >= 0)
        {
          // printf ("using: %d.%d -> %d.%d\n", i, k, ii, kk);
          //node_flood[i] = ii;
          node_flood[i] = emb->connections[3*i + k];
          arcs[3*i + k] = arcs[3*ii + kk] = 1;  // mark as SPANNING
          goon = 1;
          break;
        }
      }
    }
  }

  if (debug)
  {
    for (i = 0; i < emb->k; i++)
    {
      printf ("flood[%d] = %d.%d - arcs:", i, node_flood[i]/3, node_flood[i]%3);
      for (k = 0; k < 3; k++) printf ("%d ", arcs[3*i + k]);
      printf ("\n");
    }
  }

  /*
   * we need to know the number of underpasses for each arc
   */

  for (i = 0; i < emb->k; i++)
  {
    for (k = 0; k < 3; k++) underpasses[3*i + k] = 0;
  }

  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    for (k = 0; k < 3; k++)
    {
      if (node->direction[k] == NODE_IS_ARRIVAL) continue;
      /*
       * follow arc all the way to the next trivalent node
       */
      ii = i;
      kk = k;
      node2 = &emb->nodes[ii];
      while (1)
      {
        assert (node2->direction[kk] == NODE_IS_START);
        iinext = node2->ping[kk];
        kknext = node2->pong[kk];
        kknext = (kknext + 2) % 4;
        node2 = &emb->nodes[iinext];
        if (node2->valency == 3) break;
        if ( ((node2->overpassisodd + kknext) % 2) == 1 ) underpasses[3*i + k]++;
        ii = iinext;
        kk = kknext;
      }
      underpasses[emb->connections[3*i + k]] = underpasses[3*i + k];


      //printf ("Arc starting at %d.%d has %d underpasses\n", i, k, underpasses[3*i + k]);
    }
  }

  /*
   * now build meridian and longitude for all non-spanning arcs
   */

  if (ccemb <= 0) ccemb = 1;
  for (i = 0; i < emb->k; i++)
  {
    node = &emb->nodes[i];
    if (node->color != ccemb) continue; // do not consider unwanted components
    for (k = 0; k < 3; k++)
    {
      if (arcs[3*i + k] != 0) continue;
      if (node->direction[k] == NODE_IS_ARRIVAL) continue;

      ii = emb->connections[3*i + k]/3;
      llengthpre = numunderpasses_on_spanning_tree (i, node_flood, underpasses);
      llengthpost = numunderpasses_on_spanning_tree (ii, node_flood, underpasses);
      llength = underpasses[3*i + k] + llengthpre + llengthpost;
      // printf ("building longitude for arc starting at %d.%d (of length %d)\n", i, k, llength);
      /* adding 1 so that we do not have problems in case llength is zero */
      rule = (struct presentationrule *) malloc ((llength+1)*sizeof (int) + sizeof (struct presentationrule));
      rule->length = llength;

      //printf ("===== ARC =======\n");
      count = underpasses_on_arc (3*i+k, &(rule->var[llengthpre]), emb);

      //printf ("===== POST =======\n");
      u = llengthpre + count;
      iii = ii;
      while (iii != 0)
      {
        ikparent = node_flood[iii];
        count = underpasses_on_arc (emb->connections[ikparent], &(rule->var[u]), emb);
        u += count;
        iii = ikparent/3;
      }

      //printf ("===== PRE =======\n");
      iii = i;
      while (iii != 0)
      {
        ikparent = node_flood[iii];
        u = numunderpasses_on_spanning_tree (ikparent/3, node_flood, underpasses);
        //printf ("Moving towards root %d -> %d.%d --- u = %d\n", iii, ikparent/3, ikparent % 3, u);
        underpasses_on_arc (ikparent, &(rule->var[u]), emb);
        iii = ikparent/3;
      }
      rule->next = p->elements;
      p->elements = rule;

      //printf ("building meridian for arc starting at %d.%d\n", i, k);
      rule = (struct presentationrule *) malloc ((2*llengthpre + 1)*sizeof (int) + sizeof (struct presentationrule));
      rule->length = 2*llengthpre + 1;

      //printf ("===== PRE =======\n");
      iii = i;
      while (iii != 0)
      {
        ikparent = node_flood[iii];
        u = numunderpasses_on_spanning_tree (ikparent/3, node_flood, underpasses);
        underpasses_on_arc (ikparent, &(rule->var[u]), emb);
        iii = ikparent/3;
      }

      rule->var[llengthpre] = node->generator[k] + 1;

      //printf ("===== POST =======\n");
      u = llengthpre + 1;
      iii = i;
      while (iii != 0)
      {
        ikparent = node_flood[iii];
        count = underpasses_on_arc (emb->connections[ikparent], &(rule->var[u]), emb);
        u += count;
        iii = ikparent/3;
      }

      rule->next = p->elements;
      p->elements = rule;
    }
  }

  if (ccemb >= 0)
    if (verbose) printf ("Computing meridians and longitudes for a component of a link\n");

  free (node_flood);
  free (arcs);
  free (underpasses);
  return (emb->k/2 + 1);
}

int
emb_meridians_longitudes_torus (struct embedding *emb, struct presentation *p, int ccemb)
{
  struct emb_node *node;
  int i, ii, k, istart, idir, dir, dirplus, newdir;
  int generator, u, sign, rulenpsize, linkingnumber, rulewtsize;
  struct presentationrule *rulenp, *rule, *rulemer, *rulewt;

  if (verbose) printf ("Computing meridian and longitude for a torus-like component of a link\n");

  assert (ccemb >= emb->k);
  emb_color4 (emb);

  istart = -1;
  for (i = emb->k; i < emb->k + emb->n; i++)
  {
    node = &emb->nodes[i];
    if (node->color == ccemb)
    {
      istart = i;
      if (node->direction[0] == NODE_IS_ARRIVAL) idir = 2; else idir = 0;
      break;
    }
    if (node->colorodd == ccemb)
    {
      istart = i;
      if (node->direction[1] == NODE_IS_ARRIVAL) idir = 3; else idir = 1;
      break;
    }
  }

  i = istart;
  node = &emb->nodes[i];
  dir = idir;
  rulenpsize = emb->n + 1;
  rulenp = (struct presentationrule *) malloc (rulenpsize*sizeof (int) + sizeof (struct presentationrule));
  rulemer = (struct presentationrule *) malloc (2*sizeof (int) + sizeof (struct presentationrule));
  rulemer->length = 1;
  generator = node->generator[(dir+2)%4];
  rulemer->var[0] = generator + 1;

  u = 0;
  linkingnumber = 0;
  while (1)
  {
    sign = 1;
    if ((node->overpassisodd + dir) % 2 == 1)  // this is an underpass
    {
      dirplus = (dir + 1) % 4;
      if (node->direction[dirplus] == NODE_IS_START) sign = -1;
      assert (u < rulenpsize);
      rulenp->var[u++] = sign*(node->generator[dirplus] + 1);
      if (node->color == node->colorodd)
      {
        linkingnumber += sign;
      }
    }
    newdir = (node->pong[dir] + 2) % 4;
    i = node->ping[dir];
    node = &emb->nodes[i];
    dir = newdir;
    
    if ((i == istart) && (dir == idir)) break;
  }
  rulenp->length = u;
  assert (istart >= emb->k);

  linkingnumber = -linkingnumber;  // we have to "undo" the effect of the linkingnumber
  sign = 1;
  if (linkingnumber < 0)
  {
    sign = -1;
    linkingnumber = -linkingnumber;
  }
  rule = (struct presentationrule *) malloc ((rulenp->length+linkingnumber)*sizeof (int) + sizeof (struct presentationrule));
  for (i = 0; i < rulenp->length; i++) rule->var[i] = rulenp->var[i];
  for (i = rulenp->length; i < rulenp->length + linkingnumber; i++) rule->var[i] = sign*(generator + 1);  //XXX CHECK, perhaps negated
  rule->length = rulenp->length + linkingnumber;
  free (rulenp);

  if (globals.longitudeasrelator)
  {
    if (globals.twists == 0)
    {
      rule->next = p->rules;
      p->rules = rule;
    } else {
      rule->next = p->elements;
      p->elements = rule;
      rulewtsize = rulemer->length + rule->length*abs(globals.twists);
      rulewt = (struct presentationrule *) malloc (rulewtsize*sizeof (int) + sizeof (struct presentationrule));
      ii = 0;
      for (i = 0; i < rulemer->length; i++)
      {
        rulewt->var[ii++] = rulemer->var[i];
      }
      for (k = 0; k < abs(globals.twists); k++)
      {
        for (i = 0; i < rule->length; i++)
        {
          if (globals.twists > 0)
          {
            rulewt->var[ii++] = rule->var[i];
          } else {
            rulewt->var[ii++] = -rule->var[rule->length - i - 1];
          }
        }
      }
      assert (ii == rulewtsize);
      rulewt->length = rulewtsize;
      rulewt->next = p->rules;
      p->rules = rulewt;
      //assert (globals.twists == 0);
    }
  }
  if (globals.meridianasrelator)
  {
    rulemer->next = p->rules;
    p->rules = rulemer;
  }

  if (! globals.longitudeasrelator)
  { 
    rule->next = p->elements;
    p->elements = rule;
  }
  if (! globals.meridianasrelator)
  {
    rulemer->next = p->elements;
    p->elements = rulemer;
  }
  return (1);
}

/*
 *
 */

int
underpasses_on_arc (int i_and_k, int *var, struct embedding *emb)
{
  int i, k, inext, knext, knextplus;
  int sign;
  int u = 0;
  struct emb_node *node;

  i = i_and_k/3;
  k = i_and_k % 3;
  node = &emb->nodes[i];
  //printf ("  entering underpasses_on_arc, starting node %d.%d\n", i, k);
  while (1)
  {
    inext = node->ping[k];
    knext = node->pong[k];
    knext = (knext + 2) % 4;

    node = &emb->nodes[inext];
    //printf ("  in underpasses_on_arc, arc after crossing %d.%d, u=%d\n", inext, knext, u);
    if (node->valency == 3) break;
    knextplus = (knext + 1) % 4;
    k = knext;
    if ((node->overpassisodd + knext) % 2 == 0) continue;
    sign = 1;
    if (node->direction[knextplus] == NODE_IS_START) sign = -1;
    var[u++] = sign*(node->generator[knextplus] + 1);
  }
  //printf ("  exiting underpasses_on_arc, wrote %d chars\n", u);
  return (u);
}

/*
 * compute the number of underpasses on the spanning tree from node 0 to node i
 */

int
numunderpasses_on_spanning_tree (int i, int *node_flood, int *underpasses)
{
  int parent, k;

  if (i == 0) return (0);

  parent = node_flood[i]/3;
  k = node_flood[i] % 3;

  //printf ("UNDERPASSES_ON_SPANNING_TREE for node %d; parent: %d.%d --- underpasses to parent: %d\n", i, parent, k, underpasses[3*parent + k]);
  return (underpasses[3*parent + k] + numunderpasses_on_spanning_tree (parent, node_flood, underpasses));
}


int
emb_remove_dup_rules (struct presentation *p)
{
  int count = 0;
  struct presentationrule *r, *s;
  int k, g1, g2;

  for (r = p->rules; r; r = r->next)
  {
    if (r->length != 2) continue;
    g1 = r->var[0];
    g2 = -r->var[1];
    if (g1 < 0 || g2 < 0) continue;
    //assert (g1 > 0 && g2 > 0);
    r->length = 0;
    if (g1 == g2) continue;
    /* now substitute g2 -> g1 in all rules and selected elements */
    count++;
    for (s = p->rules; s; s = s->next)
    {
      for (k = 0; k < s->length; k++)
      {
        if (s->var[k] == g2) s->var[k] = g1;
        if (s->var[k] == -g2) s->var[k] = -g1;
      }
    }
    r->length = 1;
    r->var[0] = g2;
    for (s = p->elements; s; s = s->next)
    {
      for (k = 0; k < s->length; k++)
      {
        if (s->var[k] == g2) s->var[k] = g1;
        if (s->var[k] == -g2) s->var[k] = -g1;
      }
    }
  }
  if (debug) print_presentation (p);

  while (sp_eliminatevar (p));
  sp_removeemptyrules (p);
  if (debug) print_presentation (p);

  return (count);
}

/*
 * ==============================================
 * place construction of the orientation in a
 * function by itself
 * ==============================================
 */

int
emb_orient (struct embedding *emb)
{
  int in, nextin, id, thisid, nextid, count;
  struct emb_node *node, *nextnode, *thisnode;
  int open_arcs, short_arcs, shortcount;
  int generator, connection;

  open_arcs = emb->k/2*3;
  short_arcs = open_arcs + emb->n*2;

  shortcount = 0;

  if (debug) printf ("FASE 1: orientazione archi aperti (inizio e fine su due (o uno) nodo trivalente\n");
  if (debug) printf ("ce ne sono %d\n", open_arcs);

  /* unknown direction */
  for (in = 0; in < emb->k + emb->n; in++)
  {
    node = &emb->nodes[in];
    for (id = 0; id < node->valency; id++) node->direction[id] = 0;
  }

  if (emb->k > 0)
  {
    emb->connections = (int *) malloc (3*emb->k*sizeof(int));
    count = 0;
    for (in = 0; in < emb->k; in++)
    {
      node = &emb->nodes[in];
      for (id = 0; id < 3; id++)
      {
        //printf ("Starting from node %d direction %d\n", in, id);
        assert (node->direction[id] != NODE_IS_START);
        if (node->direction[id] != 0) continue;
        assert (count < open_arcs);
        node->direction[id] = NODE_IS_START;
        node->generator[id] = count;
        nextin = node->ping[id];
        nextnode = &emb->nodes[nextin];
        nextnode->generator[node->pong[id]] = count;
        /*
         * now follow the arc untill we reach a trivalent node
         */
        thisnode = node;
        thisid = id;
        while (1)
        {
          nextin = thisnode->ping[thisid];
          nextid = thisnode->pong[thisid];
          shortcount++;
          nextnode = &emb->nodes[nextin];
          if (thisnode->valency == 4)
          {
            generator = 2*(thisnode->id - emb->k) + open_arcs;
            if ((thisid % 2) != 0) generator++;
            thisnode->generator[thisid] = generator;
            nextnode->generator[thisnode->pong[thisid]] = generator;
          }
          if (nextnode->valency == 3)
          {
            assert (nextnode->direction[nextid] == 0);
            nextnode->direction[nextid] = NODE_IS_ARRIVAL;
            break;
          }
          nextnode->direction[nextid] = NODE_IS_ARRIVAL;
          thisnode = nextnode;
          thisid = (nextid + 2) % 4;
          assert (nextnode->direction[thisid] == 0);
          nextnode->direction[thisid] = NODE_IS_START;
        }
        emb->connections[3*in + id] = 3*nextin + nextid;
        emb->connections[3*nextin + nextid] = 3*in + id;
        count++;
      }
    }
    assert (count == open_arcs);
    if (debug)
    {
      printf ("Connections between trivalent nodes:\n");
      for (in = 0; in < emb->k; in++)
      {
        for (id = 0; id < 3; id++)
        {
          connection = emb->connections[3*in + id];
          printf ("%d.%d -> %d.%d\n", in, id, connection/3, connection % 3);
        }
      }
    }
    emb->numhcomponents = emb_color (emb);
  }

  emb->numrings = 0;

  if (short_arcs > shortcount)
  {
    if (debug) printf ("FASE 2: orientazione archi chiusi short arcs: %d, processed: %d\n", short_arcs, shortcount);
    //
    // this look very similar to the 'open arcs' case
    // perhaps we should merge the two cases?
    //
    count = 0;
    for (in = emb->k; in < emb->k + emb->n; in++)
    {
      node = &emb->nodes[in];
      assert (node->valency == 4);
      for (id = 0; id < 4; id++)
      {
        //if (debug) printf ("Starting from node %d direction %d\n", in, id);
        if (node->direction[id] != 0) continue;
        count++;
        node->direction[id] = NODE_IS_START;
        //
        // now follow the arc untill we are back at the starting node
        //
        thisnode = node;
        thisid = id;
        while (1)
        {
          assert (thisnode->valency == 4); 
          nextin = thisnode->ping[thisid];
          nextid = thisnode->pong[thisid];
          shortcount++;
          nextnode = &emb->nodes[nextin];
          assert (nextnode->valency == 4);
          generator = 2*(thisnode->id - emb->k) + open_arcs;
          if ((thisid % 2) != 0) generator++; 
          thisnode->generator[thisid] = generator;
          nextnode->generator[thisnode->pong[thisid]] = generator;
          nextnode->direction[nextid] = NODE_IS_ARRIVAL;
          thisnode = nextnode;
          thisid = (nextid + 2) % 4;
          assert (nextnode->direction[thisid] != NODE_IS_ARRIVAL);
          if (nextnode->direction[thisid] == NODE_IS_START) break;
          nextnode->direction[thisid] = NODE_IS_START;
        }
      }
    }
    assert (short_arcs == shortcount);
    emb->numrings = count;
    if (verbose && emb->numrings) printf ("There %s %d toric component(s)\n", (emb->numrings > 1)?"are":"is", emb->numrings);
  }

  return (1);
}

/*
 * ==============================================
 */

struct vecofintlist *
embeddingtoloiv (struct embedding *emb)
{
  struct emb_node *node, *thisnode, *nextnode;
  int il, in, thisin, nextin, id, thisid, nextid;
  int sign, signature, handedness;
  struct vecofintlist *loiv, *lv;
  int ic, iv, *visited;

  if (emb->k > 0)
  {
    printf ("Fatal: cannot compute dtcode/gausscode of handlebodies knots/links of higher genus\n");
    exit (1001);
  }

  loiv = 0;
  for (ic = 0; ic < emb->numrings; ic++)
  {
    lv = (struct vecofintlist *) malloc ( SIZEOFLOIV (2*emb->n) );
    lv->type = LOIV_ISGAUSSCODE;
    lv->dim = 2*emb->n;
    lv->handedness = (int *) malloc (lv->dim*sizeof(int));
    lv->next = loiv;
    loiv = lv;
  }
  loiv = lv;

  visited = (int *) malloc (2*emb->n * sizeof(int));
  for (iv = 0; iv < 2*emb->n; iv++) visited[iv] = 0;

  il = 0;
  in = 0;
  sign = emb->orientation;

  lv = loiv;
  for (in = emb->k; in < emb->k + emb->n; in++)
  {
    node = &emb->nodes[in];
    for (id = 0; id < 4; id++)
    {
      //if (debug) printf ("Starting from node %d direction %d\n", in, id);
      if (node->direction[id] != NODE_IS_START) continue;
      if (visited[2*in + (id % 2)]) break;
      assert (lv);
      //
      // now follow the arc untill we are back at the starting node
      //
      thisnode = node;
      thisid = id;
      thisin = in;
      while (1)
      {
        nextin = thisnode->ping[thisid];
        nextid = thisnode->pong[thisid];
        nextnode = &emb->nodes[nextin];

        signature = (thisnode->overpassisodd + thisid) % 2;
        signature = 1 - 2*signature;
        signature = sign*signature;   // overall orientation
        lv->vec[il] = signature*(thisin+1);
        handedness = 1;
        if (thisnode->direction[(thisid + 1) % 4] == NODE_IS_START) handedness = -1;
        lv->handedness[il] = handedness;
        visited[2*thisin + (thisid % 2)]++;

        thisin = nextin;
        thisnode = nextnode;
        thisid = (nextid + 2) % 4;
        il++;
        assert (nextnode->direction[thisid] != NODE_IS_ARRIVAL);
        if (nextin == in && (nextid + 2)%4 == id) break;
      }
      lv->len = il;
      lv = lv->next;
      il = 0;
    }
    if (visited[2*in + (id % 2)]) break;
  }

  free (visited);

  return (loiv);
}

/*
 * free memory allocated for an embedding
 */

void
freeembedding (struct embedding *emb)
{
  free (emb->nodes);
  if (emb->connections) free (emb->connections);

  free (emb);

  return;
}

/*
 * print simplifying "Reidemeister" rules that simplify embedding
 * the "typeII" move applies as in this diagram
 *
 *  |
 * ----.
 *  |  |
 *  `--|--
 *     |
 *
 * fork:
 *
 *  |
 * ----.
 *  |  |
 *  *--|--
 * /   |
 *
 * twist:
 *      ___
 *     |   |
 *   ----. |
 *     | | |
 *   ------'
 *     | |
 *
 * loop-flip
 *      ___
 *     |   |
 *   ----. |
 *     | | |
 *     *---'
 *    /  |
 *
 * barrel
 *
 *       |
 *    .--|-.
 * --*   |  *--
 *    `----'
 *       |
 *
 * twist3
 *
 *          /
 *     ,---*
 *     |   |
 *   ----. |
 *     | | |
 *   ------'
 *     | |
 *
 * twist4
 *
 *         | 
 *     ,---|--
 *     |   |
 *   ----. |
 *     | | |
 *   ------'
 *     | |
 * NOTE: applying a loop-flip or barrel move (change both overpasses) does not change the planar
 * embedding.  However the value of 'choice' changes.  We say that the move is 'simplifying'
 * if applying it decreases the 'choice' value.
 */

int check_for_loop_flip (struct embedding *emb, struct dualembedding *dual, struct dual_region *r);
int check_for_twist (struct embedding *emb, struct dualembedding *dual, struct dual_region *r);
int check_for_barrel (struct embedding *emb, struct dualembedding *dual, struct dual_region *r);
int check_for_twist3 (struct embedding *emb, struct dualembedding *dual, struct dual_region *r);
int check_for_twist4 (struct embedding *emb, struct dualembedding *dual, struct dual_region *r);

void
printembrules (struct embedding *emb, struct dualembedding *dual)
{
  struct dual_region *region;
  struct emb_node *node;
  int jj, val, parity, countcrossings;
  int inode;
  int isbetween, between, otherside;

  /*
   * does reverting all overpasses decrese 'choice'?
   * It is sufficient to check the overpass of the last node
   */

  inode = emb->k + emb->n - 1;
  node = &(emb->nodes[inode]);
  if (node->overpassisodd)
  {
    printf ("mirror\n");
    return;  // no more checks in this case
  }

  /*
   * first: search for typeII
   */

  for (region = dual->regions; region; region = region->next)
  {
    //printf ("%d:(", region->id);
    if (check_for_twist3 (emb, dual, region)) printf ("twist3\n");
    if (check_for_twist4 (emb, dual, region)) printf ("twist4\n");
    val = region->valency;
    assert (val >= 2);
    parity = 0;
    countcrossings = 0;
    isbetween = between = 0;
    for (jj = 0; jj < val; jj++)
    {
      inode = region->wedgeij[jj]/4;
      node = &(emb->nodes[inode]);
      if (node->valency != 4)
      {
        assert (node->valency == 3);
        if (isbetween) between++;
        continue;
      }
      countcrossings++;
      isbetween = 1 - isbetween;
      parity += node->overpassisodd + (region->wedgeij[jj]%4);
    }

    //printf ("<%d.%d>\n", region->wedgeij[0]/4, region->wedgeij[0]%4);
    //printf ("<%d.%d>\n", region->wedgeij[1]/4, region->wedgeij[1]%4);
    //printf ("parity: %d\n", parity);
    otherside = val - countcrossings - between;
    //printf ("region %d: %d = %d + %d + %d (val = crossings + between + otherside)\n", region->id, val, countcrossings, between, otherside);
    if (countcrossings != 2) continue;
    assert (otherside >= 0);
    if (otherside > 0 && between > 0) continue;
    if ( countcrossings == 2 && (parity%2) == 1)
    {
      if (val == 2) printf ("typeII\n");
       else printf ("fork\n");
    }
    /*
     * explanation: a region with exactly 2 crossings that are adjacent is
     * either a bigon (with the obvious typeII simplifying move for 2 out of 4 choices of the overpasses)
     * or (possibly after IH moves) that can allow for a "fork" type simplifying move.
     *
     * The parity computation allows to pick the overpasses choices that allow the moves
     */
    if (check_for_loop_flip (emb, dual, region)) printf ("loop-flip\n");
    if (check_for_twist (emb, dual, region)) printf ("twist\n");
    if (check_for_barrel (emb, dual, region)) printf ("barrel\n");
  }

  return;
}

int
check_for_loop_flip (struct embedding *emb, struct dualembedding *dual, struct dual_region *r)
{
  struct dual_region *r01;
  struct emb_node *nodekk, *lastother, *node[2];
  int inodekk, inode[2];
  int jj, kk, parity;
  int ih, il, ilast, ilastother;

  if (r->valency != 2) return (0);
  for (jj = 0; jj < 2; jj++)
  {
    inode[jj] = r->wedgeij[jj]/4;
    node[jj] = &(emb->nodes[inode[jj]]);
  }

  ilast = emb->k + emb->n - 1;
  ih = 0; il = 1;
  if (inode[ih] < inode[il]) { ih = 1; il = 0; }

  if (emb->n >= 3 && inode[ih] == ilast)
  {
    assert (node[ih]->overpassisodd == 0); // otherwise we got "mirror" above
    ilastother = ilast - 1;
    if (ilastother == inode[il]) ilastother--;
    lastother = &(emb->nodes[ilastother]);
    if (lastother->overpassisodd == 0) return (0);
  } else {
    if (node[ih]->overpassisodd == 0) return (0);
    // value of choice cannot be decreased
  }

  parity =  node[0]->overpassisodd + (r->wedgeij[0]%4);
  parity += node[1]->overpassisodd + (r->wedgeij[1]%4);
  if ((parity % 2) == 1) return (0); // typeII not allowed here

  for (jj = 0; jj < 2; jj++)
  {
    r01 = r->ping[jj];
    if (r01->valency == 3)
    {
      for (kk = 0; kk < 3; kk++)
      {
        inodekk = r01->wedgeij[kk]/4;
        nodekk = &(emb->nodes[inodekk]);
        if (nodekk->valency == 3) return (1);
      }
    }
  }

  //printf ("n1: %d, n2: %d\n", r->wedgeij[0]/4, r->wedgeij[1]/4);
  return (0);
}

/*
 * check for a situation like this:
 *      ___
 *     |   |
 *   ----. |
 *     | | |
 *   ------'
 *     | |
 */

int
check_for_twist (struct embedding *emb, struct dualembedding *dual, struct dual_region *r)
{
  struct dual_region *r01;
  struct emb_node *nodekk, *node[2];
  int inodekk, inode[2];
  int jj, kk, parity, par[3];
  int totval;

  if (r->valency != 2) return (0);
  //if (verbose) printf ("Checking for twist on bigon %d\n", r->id);
  for (jj = 0; jj < 2; jj++)
  {
    inode[jj] = r->wedgeij[jj]/4;
    node[jj] = &(emb->nodes[inode[jj]]);
  }

  parity =  node[0]->overpassisodd + (r->wedgeij[0]%4);
  parity += node[1]->overpassisodd + (r->wedgeij[1]%4);
  if ((parity % 2) == 1) return (0); // typeII not allowed here

  for (jj = 0; jj < 2; jj++)
  {
    r01 = r->ping[jj];
    if (r01->valency == 3)
    {
      totval = 0;
      for (kk = 0; kk < 3; kk++)
      {
        inodekk = r01->wedgeij[kk]/4;
        nodekk = &(emb->nodes[inodekk]);
        par[kk] = ( nodekk->overpassisodd + (r01->wedgeij[kk]%4) ) % 2;
        totval += nodekk->valency;
      }
      if (totval == 12)
      {
        if (par[0] != par[1] || par[1] != par[2]) return (1);
      }
    }
  }

  return (0);
}

/*
 * check for a situation like this:
 *       |
 *    .--|-.
 * --*   |  *--
 *    `----'
 *       |
 */

int
check_for_barrel (struct embedding *emb, struct dualembedding *dual, struct dual_region *r)
{
  struct dual_region *radj;
  struct emb_node *nodekk, *lastother, *node[3];
  int inodekk, inode[3];
  int jj, parity, jtri, jtriadj, jj1, jj2;
  int ih, il, ilast, ilastother;

  if (r->valency != 3) return (0);
  jtri = -1;
  for (jj = 0; jj < 3; jj++)
  {
    inode[jj] = r->wedgeij[jj]/4;
    node[jj] = &(emb->nodes[inode[jj]]);
    if (node[jj]->valency == 3)
    {
      if (jtri >= 0) return (0);
      jtri = jj;
    }
  }
  if (jtri < 0) return (0);
  //fprintf (stderr, "nodeid: %d\n", node[jtri]->id);

  jj1 = (jtri + 1)%3;
  jj2 = (jtri + 2)%3;
  radj = r->ping[jj2];
  if (radj->valency != 3) return (0);
  jtriadj = -1;
  for (jj = 0; jj < 3; jj++)
  {
    inodekk = radj->wedgeij[jj]/4;
    nodekk = &(emb->nodes[inodekk]);
    if (nodekk->valency == 3)
    {
      jtriadj = jj;
      //fprintf (stderr, "nodeadjid: %d\n", nodekk->id);
      if (nodekk->id < node[jtri]->id) return (0);  // only one output
      break;
    }
  }
  if (jtriadj < 0) return (0);
  parity =  node[jj1]->overpassisodd + (r->wedgeij[jj1]%4);
  parity += node[jj2]->overpassisodd + (r->wedgeij[jj2]%4);
  parity = parity % 2;
  if (parity != 0) return (0);

  /*
   * we are in a "barrel" situation
   * now we signal it only if the choice value will decrease
   */

  //fprintf (stderr, "region: %d, radj: %d, jj1=%d jj2=%d jtri=%d\n", r->id, radj->id, jj1, jj2, jtri);
  //fprintf (stderr, "parity: %d\n", parity%2);

  ilast = emb->k + emb->n - 1;
  ih = jj1; il = jj2;
  if (inode[ih] < inode[il]) { ih = jj2; il = jj1; }

  if (emb->n >= 3 && inode[ih] == ilast)
  {
    assert (node[ih]->overpassisodd == 0); // otherwise we got "mirror" above
    ilastother = ilast - 1;
    if (ilastother == inode[il]) ilastother--;
    lastother = &(emb->nodes[ilastother]);
    if (lastother->overpassisodd == 0) return (0);
  } else {
    if (node[ih]->overpassisodd == 0) return (0);
    // value of choice cannot be decreased
  }

  return (1);
}

/*
 * check for a situation like this:
 *
 *          /
 *     ,---*
 *     |   |
 *   ----. |
 *     | | |
 *   ------'
 *     | |
 */

int
check_for_twist3 (struct embedding *emb, struct dualembedding *dual, struct dual_region *r)
{
  struct dual_region *radj;
  struct emb_node *nodejj, *node[3];
  int inodejj, inode[3];
  int jj, jj2, jopp, parity[3];
  int totparity;

  if (r->valency != 3) return (0);
  totparity = 0;
  for (jj = 0; jj < 3; jj++)
  {
    inode[jj] = r->wedgeij[jj]/4;
    node[jj] = &(emb->nodes[inode[jj]]);
    parity[jj] = (r->wedgeij[jj] + node[jj]->overpassisodd) % 2;
    if (node[jj]->valency != 4) return (0);
    totparity += parity[jj];
  }
  if (totparity < 1 || totparity > 2) return (0);
  jopp = -1;
  for (jj = 0; jj < 3; jj++)
  {
    if (((totparity + parity[jj]) % 2) == 0)
    {
      assert (jopp < 0);
      jopp = jj;
    }
  }
  assert (jopp >= 0);

  //jj1 = (jopp + 1)%3;
  jj2 = (jopp + 2)%3;
  radj = r->ping[jj2];
  if (radj->valency != 3) return (0);
  //fprintf (stderr, "r: %d, radj: %d\n", r->id, radj->id);

  for (jj = 0; jj < 3; jj++)
  {
    inodejj = radj->wedgeij[jj]/4;
    nodejj = &(emb->nodes[inodejj]);
    if (nodejj->valency == 3) return (1);
  }

  return (0);
}

/*
 * check for a situation like this:
 *
 *         |
 *     ,---|--
 *     |   |
 *   ----. |
 *     | | |
 *   ------'
 *     | |
 */

int
check_for_twist4 (struct embedding *emb, struct dualembedding *dual, struct dual_region *r)
{
  struct dual_region *radj;
  struct emb_node *nodejj, *node[3];
  int inodejj, inode[3];
  int jj, jj2, jopp, parity[3];
  int totparity;

  if (r->valency != 3) return (0);
  totparity = 0;
  for (jj = 0; jj < 3; jj++)
  {
    inode[jj] = r->wedgeij[jj]/4;
    node[jj] = &(emb->nodes[inode[jj]]);
    parity[jj] = (r->wedgeij[jj] + node[jj]->overpassisodd) % 2;
    if (node[jj]->valency != 4) return (0);
    totparity += parity[jj];
  }
  if (totparity < 1 || totparity > 2) return (0);
  jopp = -1;
  for (jj = 0; jj < 3; jj++)
  {
    if (((totparity + parity[jj]) % 2) == 0)
    {
      assert (jopp < 0);
      jopp = jj;
    }
  }
  assert (jopp >= 0);

  //jj1 = (jopp + 1)%3;
  jj2 = (jopp + 2)%3;
  radj = r->ping[jj2];
  if (radj->valency != 3) return (0);
  //fprintf (stderr, "r: %d, radj: %d\n", r->id, radj->id);

  if ((parity[jopp] % 2) == 1) return (0);
  //fprintf (stderr, "r: %d, parityopp = %d\n", r->id, parityopp);
  for (jj = 0; jj < 3; jj++)
  {
    inodejj = radj->wedgeij[jj]/4;
    nodejj = &(emb->nodes[inodejj]);
    if (nodejj->valency == 3) return (0);
    if ((( /* parityopp + */ radj->wedgeij[jj] + nodejj->overpassisodd) % 2) == 1) return (1);
  }

  return (0);
}

/*
 * a region description of a ``thin'' apparent contour can be converted to an embedding!
 */

/*
 * local functions
 */

void ls2e_adjust_cyclic_lists (struct embedding *emb);
struct embedding *ls2e_sanity_check (struct sketch *s);
int ls2e_getadj3node (struct embedding *emb, int i, struct emb_node *node, int j, struct border *transbp);
int ls2e_build3nodenet (struct embedding *emb, int ik);

static struct region **n2r = 0;
static int *r2n = 0;

/* main function */

struct embedding *
trysketch2embedding (struct sketch *s)
{
  struct embedding *emb;
  struct emb_node *node;
  int i, j, k, kcount;
  struct region *r;
  struct borderlist *border;
  struct border *bp, *transbp;

  emb = ls2e_sanity_check (s);
  if (emb == 0) return (0);

  if (verbose) printf ("number of 3-vertices: %d\n", emb->k);
  if (verbose) printf ("number of crossings: %d\n", emb->n);

  n2r = (struct region **) malloc ((emb->k + emb->n)*sizeof(struct region *));
  r2n = (int *) malloc (s->regioncount*sizeof(int));

  i = kcount = 0;
  for (r = s->regions; r; r = r->next)
  {
    if (r->f == 4)
    {
      n2r[i+emb->k] = r;
      r2n[r->tag] = i;
      i++;
    }
    if (r->f == 2)
    {
      r2n[r->tag] = 0;
      bp = r->border->sponda;
      do {
        r2n[r->tag]++;
        bp = bp->next;
        bp = bp->next;
      } while (bp != r->border->sponda);
      if (r2n[r->tag] > 2)
      {
        for (k = 0; k < r2n[r->tag] - 2; k++)
        {
          n2r[kcount] = r;
          kcount++;
        }
      }
    }
  }
  assert (kcount == emb->k);

  /* set adj of 3-nodes to an invalid value */

  for (i = 0; i < emb->k; i++)
  {
    node = &(emb->nodes[i]);
    assert (node->valency == 3);
    for (j = 0; j < 3; j++) node->ping[j] = i;
  }

  for (k = 0; k < emb->k; k++)
  {
    ls2e_build3nodenet (emb, k);
  }

  for (i = 0; i < emb->n; i++)
  {
    node = &(emb->nodes[i + emb->k]);
    border = n2r[i+emb->k]->border;
    bp = border->sponda;
    if (bp->info->depths[0] == 0)
    {
      node->overpassisodd = 1;
      emb->choice += (1 << i);
    }
    for (j = 0; j < 4; j++)
    {
      transbp = gettransborder (bp);
      if (r2n[transbp->border->region->tag] > 2)
      {
        node->ping[j] = ls2e_getadj3node (emb, i, node, j, transbp);
      } else {
        transbp = gettransborder(transbp->next->next);
        assert (transbp->border->region->f == 4);
        node->ping[j] = r2n[transbp->border->region->tag] + emb->k;
      }
      bp = bp->next;
    }
  }

  ls2e_adjust_cyclic_lists (emb);

  //printembedding (emb);
  free (n2r);
  free (r2n);
  return (emb);
}

/*
 * build 3-node net
 */

int
ls2e_build3nodenet (struct embedding *emb, int ik)
{
  int j, rtag, crossing;
  struct region *r;
  struct border *bn, *bntrans;
  struct emb_node *node;

  node = &(emb->nodes[ik]);
  r = n2r[ik];
  rtag = r->tag;

  if (node->ping[0] != ik) return(0);  // work already done for this 3-node

  assert (r2n[rtag] >= 3);

  bn = r->border->sponda;
  bntrans = gettransborder (bn);
  if (bntrans->border->region->f != 4) bn = bn->next;
  assert (gettransborder(bn)->border->region->f == 4);
  bntrans = gettransborder (bn);
  crossing = r2n[bntrans->border->region->tag];
  node->ping[0] = crossing + emb->k;
  bn = bn->next->next;

  for (j = 0; j < r2n[rtag] - 2; j++)
  {
    bntrans = gettransborder (bn);
    crossing = r2n[bntrans->border->region->tag];
    node = &(emb->nodes[ik+j]);
    assert (n2r[ik] == n2r[ik+j]);
    if (j > 0) node->ping[0] = ik + j - 1;
    if (j < r2n[rtag] - 3) node->ping[2] = ik + j + 1;
    node->ping[1] = crossing + emb->k;
    bn = bn->next->next;
  }

  bntrans = gettransborder (bn);
  crossing = r2n[bntrans->border->region->tag];
  node->ping[2] = crossing + emb->k;

  return (1);
}

/*
 * crossing i is adjacent through j to a 3-valent node
 * return the id of such 3-valent node
 */

int
ls2e_getadj3node (struct embedding *emb, int i, struct emb_node *node, int j, struct border *transbp)
{
  int ik, iik, rtag;
  struct emb_node *knode;

  assert (node == &(emb->nodes[i+emb->k]));
  rtag = transbp->border->region->tag;
  assert (r2n[rtag] > 2);

  /* find id of (first) 3-node */
  for (ik = 0; ik < emb->k; ik++)
  {
    if (n2r[ik] == transbp->border->region) break;
  }
  assert (ik < emb->k);

  if (r2n[rtag] == 3) return (ik);

  for (iik = ik; iik < ik + r2n[rtag] - 2; iik++)
  {
    assert (n2r[ik] == n2r[iik]);
    knode = &(emb->nodes[iik]);
    for (j = 0; j < 3; j++)
    {
      if (knode->ping[j] == i + emb->k) return (iik);
    }
  }

  fprintf (stderr, "Something went wrong: i = %d, j = %d, r2n[%d] = %d - ik = %d\n", i+emb->k, j, rtag, r2n[rtag], ik);
  fprintf (stderr, "case of network with multiple 3-nodes (there are %d) not implemented\n", r2n[rtag] - 2);
  return (0);
}

/*
 * move here all sanity checks
 * return a pointer to an embedding if they pass
 */

struct embedding *
ls2e_sanity_check (struct sketch *s)
{
  struct embedding *emb;
  struct emb_node *node;
  int i;
  int count, xcount = 0, kcount = 0, acount = 0, rcount = 0;
  int euler;
  struct region *r;
  struct borderlist *border;
  struct border *bp, *bpstart, *transbp;
  struct arc *arc;

  for (r = s->regions; r; r = r->next)
  {
    border = r->border;
    switch (r->f)
    {
      case 0:
        rcount++;
        break;

      case 2:
        if (border->next != 0)
        {
          if (verbose) printf ("Invalid non-simply-connected arc region\n");
          return (0);
        }
        bpstart = border->sponda;
        transbp = gettransborder (bpstart);
        if (transbp->border->region->f != 4) bpstart = bpstart->next;
        bp = bpstart;
        count = 0;
        do {
          count++;
          transbp = gettransborder (bp);
          if (transbp->border->region->f != 4)
          {
            if (verbose) printf ("Invalid f = %d while walking around an 'arc' region\n", transbp->border->region->f);
            return (0);
          }
          bp = bp->next;
          assert (bp != bpstart);
          transbp = gettransborder (bp);
          if (transbp->border->region->f != 0)
          {
            if (verbose) printf ("Invalid f = %d while walking around an 'arc' region\n", transbp->border->region->f);
            return (0);
          }
          bp = bp->next;
        } while (bp != bpstart);
        assert (count >= 2);
        acount++;
        kcount += (count - 2);
        acount += 2*(count-2);
        break;

      case 4:
        xcount++;
        if (border->next != 0)
        {
          if (verbose) printf ("Invalid non-simply-connected vertex region\n");
          return (0);
        }
        bp = border->sponda;
        count = 0;
        do {
          count++;
          assert (bp->border->region == r);
          if (bp->orientation != 1)
          {
            if (verbose) printf ("Orientation must be 1\n");
            return (0);
          }
          bp = bp->next;
        } while (bp != border->sponda);
        if (count != 4)
        {
          if (verbose) printf ("a 'vertex' region must have 4 sides.  There are %d instead\n", count);
          return (0);
        }
        break;

      default:
        if (verbose) printf ("Invalid apparent contour with a region with f=%d\n", r->f);
        return (0);
    }
  }

  /*
   * sanity check of arcs: no cusps allowed
   */

  for (arc = s->arcs; arc; arc = arc->next)
  {
    if (arc->cusps != 0)
    {
      if (verbose) printf ("No cusps allowed in region description, arc %d has %d cusps\n", arc->tag, arc->cusps);
      return (0);
    }
    assert (arc->endpoints == 1 || arc->endpoints == 2);
    if (arc->endpoints == 1)
    {
      printf ("Found a loop! This is not allowed in an embedding description\n");
      return (0);
    }
  }

  euler = xcount + kcount + rcount - acount;
  if (euler != 2)
  {
    if (verbose) printf ("Euler characteristic is %d instead of 2.  Perhaps this is a split diagram\n", euler);
    return (0);
  }

  /*
   * all sanity checks passed!
   */

  emb = (struct embedding *) malloc (sizeof (struct embedding));
  emb->k = kcount;
  emb->n = xcount;
  emb->choice = 0;
  emb->orientation = 1;
  emb->connections = 0;
  emb->nodes = (struct emb_node *) malloc ((emb->k+emb->n)*sizeof (struct emb_node));

  for (i = 0; i < emb->k + emb->n; i++)
  {
    node = &(emb->nodes[i]);
    node->next = 0;
    if (i+1 < emb->k+emb->n) node->next = &(emb->nodes[i+1]);
    node->id = i;
    node->valency = 4;
    if (i < emb->k) node->valency = 3;
  }

  if (verbose) printf ("number of arcs: %d\n", acount);
  if (verbose) printf ("number of regions: %d\n", rcount);

  return (emb);
}

void
ls2e_adjust_cyclic_lists (struct embedding *emb)
{
  int i, j, minval, maxval, jrot;
  int pingsave[4];
  struct emb_node *node;

  for (i = emb->k; i < emb->k + emb->n; i++)
  {
    node = &(emb->nodes[i]);
    minval = 9999;
    maxval = -1;
    for (j = 0; j < 4; j++)
    {
      if (node->ping[j] < minval) minval = node->ping[j];
      if (node->ping[j] > maxval) maxval = node->ping[j];
      pingsave[j] = node->ping[j];
    }
    if (minval == maxval) continue;
    for (j = 0; j < 4; j++)
    {
      if (node->ping[j] == minval && node->ping[(j+3)%4] != minval) {jrot = j; break;}
    }
    for (j = 0; j < 4; j++)
    {
      node->ping[j] = pingsave[(j + jrot) % 4];
    }
    if ((jrot % 2) == 1)
    {
      node->overpassisodd = !(node->overpassisodd);
      emb->choice ^= (1<<(i-emb->k));
    }
  }
}
