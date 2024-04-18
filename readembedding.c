#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "contour.h"
#include "readembedding.h"
#include "fundamental.h"
#include "parser.h"
#include "readdtcode.h"
#include "wirtinger.h"

extern struct global_data globals;

extern int debug;
extern int verbose;
extern int quiet;

struct sketch *
embedding2sketch (struct embedding *emb)
{
  int i, ir, found;
  struct emb_node *node;
  struct sketch *s;
  struct vecofintlist *loiv;
  struct dualembedding *dual;

  if (emb->k == 0)
  {
    if (debug) printf ("This is a standard knot/link, converting to gausscode\n");
    loiv = embeddingtoloiv (emb);
    if (debug) printloiv (loiv);
    freeembedding (emb);
    s = readgausscodeloiv (loiv);
    // freeloiv (loiv); // loiv is freed by readgausscodeloiv!
    return (s);
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
  for (ir = 0; ir < dual->numregions; ir++)
  {

printf ("region %d\n", ir);
    
//    fourvalentnum = 0;
//    for (j = 0; j < r_adjacencynum[ir]; j++)
//    {
//      /*
//       * check if the corresponding node is 4-valent
//       */
//      inode = r_adjacencyi[ir][j];
//      if (adjacencynum[inode] == 4) fourvalentnum++;
//    }
//    if (fourvalentnum == 0) found++;
  }

  if (found > 0)
  {
    printf ("\n# WARNING: presence of Regions with no 4-valent bound vertex!\n");
    printf ("# the code at present cannot cope with this situation.\n");
    printf ("#\n# However in this case the embedding is IH equivalent to a one-sum\n");
    printf ("# so we are not really interested in such a case.\n");
    printf ("#\n# we simply output an empty contour\n\n");
    printf ("sketch {\n");
    printf ("Region 0 (f = 0): ();\n");
    printf ("}\n");
    return (0);
  }

#if IMPORTEDFROMEMBEDDINGANALYZE

  printf ("sketch {\n");

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
  for (i = 1; i <= 4*n; i++) printf ("Arc %d: [0];\n", i);

  /*
   * arcs bounding the small quadrilateral crossing regions
   * (four per node), listed counter-clockwise.
   * the label value defines the under/over choice and can
   * be driven using the bits in the "choice" argument
   */
  count = 4*n+1;
  for (i = 1; i <= n; i++)
  {
    bit = choice % 2;
    choice /= 2;

    printf ("Arc %d: [%d];\n", count++, 2*bit);
    printf ("Arc %d: [%d];\n", count++, 2*(1 - bit));
    printf ("Arc %d: [%d];\n", count++, 2*bit);
    printf ("Arc %d: [%d];\n", count++, 2*(1 - bit));
  }

  count = print_sketch_rr (0);
  count = print_sketch_ra (count);
  count = print_sketch_rv (count);

  printf ("}\n");

int
print_sketch_rr (int count)
{
  int arcid, notfirst;
  int ir, jr, iv, jv;

  /*
   * start with real regions (f = 0)
   */
  for (ir = 0; ir < numregions; ir++)
  {
    printf ("Region %d (f = 0): (", count);
    if (count++ == 1) printf (") (");
    notfirst = 0;
    for (jr = 0; jr < r_adjacencynum[ir]; jr++)
    {
      iv = r_adjacencyi[ir][jr];
      if (adjacencynum[iv] <= 3) continue;
      if (notfirst++) printf (" ");
      jv = r_adjacencyj[ir][jr];
      jv++; if (jv >= 4) jv -= 4;
      arcid = 4*(iv - k) + jv + 1;
      printf ("-a%d", arcid);
    }
    printf (");\n");
  }
  return (count);
}

/*
 * then sketch regions corresponding to arcs
 */

int
print_sketch_ra (int count)
{
  int i;

  nodemark = (int *) malloc (nodenum*sizeof(int));
  for (i = 0; i < k+n; i++) nodemark[i] = 0;

  /*
   * now the elongated regions containint the triple-points
   * (one or two)
   */
  while (print_sketch_ra3 (1, count, nodemark))
  {
    count++;
  }


  while (print_sketch_ra3 (0, count, nodemark))
  {
//for (i = 0; i < nodenum; i++) printf ("nodemark[%d] = %d\n", i, nodemark[i]);
    count++;
  }

  free (nodemark);

  return (count);
}

int
print_sketch_ra3 (int only3, int regionnum, int *nodemark)
{
  int vi, vj, vinext, vjnext, vinextnext;
  int found, notfirst;

  /*
   * now the elongated regions containint the triple-points
   * (one or two)
   */
  found = 0;
  for (vi = k; vi < k+n; vi++)
  {
    for (vj = 0; vj < 4; vj++)
    {
      if ((nodemark[vi] & (1<<vj)) != 0) continue;
      if (only3 && adjacency[vi][vj] >= k) continue;
      found++;
      break;
    }
    if (found) break;
  }

  if (! found) return (0);

  printf ("Region %d (f = 2): (", regionnum);
  notfirst = 0;
  while ((nodemark[vi] & (1<<vj)) == 0)
  {
    nodemark[vi] |= 1 << vj;
//printf ("\n=== nodemark[%d] |= 1<<%d\n", vi, vj);

    vinext = adjacency[vi][vj];
    vjnext = get_back_index (vi, vj);
    //nodemark[vinext] |= 1 << vjnext;
    vjnext++; if (vjnext >= adjacencynum[vinext]) vjnext -= adjacencynum[vinext];
    //printf ("LOOP: vi.vj = %d.%d -> %d.%d\n", vi, vj, vinext, vjnext);
//    nodemark[vinext] |= 1 << vjnext;
//printf ("\n=== nodemark[%d] |= 1<<%d\n", vinext, vjnext);
    while (adjacencynum[vinext] <= 3)
    {
      nodemark[vinext] |= 1 << vjnext;
//printf ("\n=== nodemark[%d] |= 1<<%d\n", vinext, vjnext);
      vinextnext = adjacency[vinext][vjnext];
      vjnext = get_back_index (vinext, vjnext);
      //nodemark[vinextnext] |= 1 << vjnext;
      vjnext++; if (vjnext >= adjacencynum[vinextnext]) vjnext -= adjacencynum[vinextnext];
//      nodemark[vinextnext] |= 1 << vjnext;
//printf ("=== nodemark[%d] |= 1<<%d\n", vinextnext, vjnext);
      vinext = vinextnext;
      //printf ("  -> %d.%d\n", vinext, vjnext);
    }
    assert (adjacencynum[vinext] == 4);
    if (notfirst++) printf (" ");
    printf ("+a%d ", 4*(vi - k) + vj + 1);
    printf ("-a%d", 4*(vinext - k) + vjnext + 1 + 4*n);

    vjnext--;
    if (vjnext < 0) vjnext += 4;
    vi = vinext;
    vj = vjnext;
  }

  printf (");\n");
  return (1);
}

int
print_sketch_rv (int count)
{
  int i, j;

  /*
   * finally regions corresponding to four-valent nodes
   */

  for (i = 0; i < n; i++)
  {

    printf ("Region %d (f = 4): (", count++);
    for (j = 0; j < 4; j++)
    {
      if (j > 0) printf (" ");
      printf ("+a%d", 4*n + 4*i + j + 1);
    }
    printf (");\n");
  }

  return (count);
}
#endif //IMPORTEDFROMEMBEDDINGANALYZE






  printf ("NOT YET IMPLEMENTED: k = %d, n = %d, choice = %d\n", emb->k, emb->n, emb->choice);

  freedual (dual);
  return (0);
}

/*
 *
 */

void
freedual (struct dualembedding *dual)
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
  int i, j, jj, val, nodenum, numwedges, count, totcount;
  int regionsize, region_id;
  int inode, jwedgeminus, jwedgeplus, found;
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
printf ("region %d: starting from wedge %d.%d\n", region_id, i, j);
      //count = around_region (emb, i, j, &r_wedgeij_all[totcount]);
      r_wedgeij = &r_wedgeij_all[totcount];
      /*
       * walk around the region
       */
      count = 0;
      thisnode = node;
      thisnodej = j;

      while ((nextnodei = thisnode->ping[thisnodej]) != i)
      {
        count++;
        nextnode = &emb->nodes[nextnodei];
        val = nextnode->valency;
        nextnodej = thisnode->pong[thisnodej];
        nextnodej = (nextnodej + val - 1) % val;
printf ("  next wedge: %d.%d\n", nextnodei, nextnodej);

        wedgemark[4*nextnodei + nextnodej] = 1;
        //assert (jwedge >= 0);
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
      jwedgeminus = region->wedgeij[j] % 4;
      jwedgeplus = (jwedgeminus + 1) % 4;
      /* search for the neighbouring region */
      //found = 0;
      for (adjregion = dual->regions; adjregion; adjregion = adjregion->next)
      {
        jj = jwedgeplus;
        if (adjregion->wedgeij[jj]/4 == inode) continue;
        break;
      }
      region->ping[j] = region;
    }
  }

  free (wedgemark);
  return (dual);
}

void
printdual (struct dualembedding *dual)
{
  struct dual_region *region, *adjregion;
  int jj;

  if (verbose) printf ("Dual of given embedding (%d regions):\n", dual->numregions);

  printf ("dualembedding {");
  for (region = dual->regions; region; region = region->next)
  {
    printf ("%d:(", region->id);
    for (jj = 0; jj < region->valency; jj++)
    {
      adjregion = region->ping[jj];
      printf ("%d, ", adjregion->id);
    }
    printf ("), ");
  }
  printf ("}\n");
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
  int open_arcs, short_arcs;
  struct emb_node *node;
  int i, j;
  struct presentation *p;
  struct presentationrule *rule;
  int sign, a, b, c, d;
  struct vecofintlist *gaussloiv;
  struct vecofintlist *dtloiv;

  assert ((emb->k % 2) == 0);
  open_arcs = emb->k/2*3;
  short_arcs = open_arcs + emb->n*2;

  assert (emb->orientation);

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
        printf (" %d", node->generator[j]);
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

  assert (emb->numhcomponents >= 1);
  if (emb->numhcomponents + emb->numrings == 1)
  {
    assert (emb->numrings == 0);
    emb_meridians_longitudes (emb, p);
  } else if (!quiet) printf ("Cannot compute meridians and longitudes for links\n");

  if (debug) print_presentation (p);
  emb_remove_dup_rules (p);
  if (globals.simplifypresentation) simplify_presentation (p);
  return (p);
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

int
emb_meridians_longitudes (struct embedding *emb, struct presentation *p)
{
  int i, ii, iii, kk, iinext, kknext;
  int ikparent;
  int k, goon;
  int *node_flood, *arcs, *underpasses;
  struct emb_node *node, *node2;
  struct presentationrule *rule;
  int llength, llengthpre, llengthpost;
  int u, count;

  /*
   * step 1. construct a spanning tree for the spatial graph
   * rooted at node 0
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

  node_flood[0] = 0;
  goon = 1;
  while (goon)
  {
    goon = 0;
    for (i = 1; i < emb->k; i++)
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

  for (i = 0; i < emb->k; i++)
  {
    for (k = 0; k < 3; k++)
    {
      if (arcs[3*i + k] != 0) continue;
      node = &emb->nodes[i];
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

  free (node_flood);
  free (arcs);
  free (underpasses);
  return (emb->k/2 + 1);
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
    assert (g1 > 0 && g2 > 0);
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


