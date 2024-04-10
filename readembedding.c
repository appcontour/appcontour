#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "contour.h"
#include "readembedding.h"
#include "fundamental.h"
#include "parser.h"

extern struct global_data globals;

struct sketch *
readembedding (FILE *file)
{
  struct embedding *emb;

  emb = readembedding_low (file);

  printf ("NOT YET IMPLEMENTED: k = %d, n = %d, choice = %d\n", emb->k, emb->n, emb->choice);


  return (0);
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

extern int debug;
extern int verbose;
extern int quiet;

struct presentation *
wirtingerfromembedding (struct embedding *emb)
{
  int open_arcs, short_arcs;
  struct emb_node *node, *thisnode, *nextnode;
  int inode, nextinode;
  int i, j, thisi, nexti;
  int thisj, nextj, nextjplus, k;
  int count, shortcount, numrings, numhcomponents;
  struct presentation *p;
  struct presentationrule *rule;
  int sign, generator, a, b, c, d;
  int connection, *connections;

  assert ((emb->k % 2) == 0);
  open_arcs = emb->k/2*3;
  short_arcs = open_arcs + emb->n*2;
  numhcomponents = 0;

  assert (emb->k <= 4); /* for now we restrict to max 4 trivalent nodes */
  shortcount = 0;
  for (inode = 0; inode < emb->k + emb->n; inode++)
  {
    node = &emb->nodes[inode];
    for (i = 0; i < node->valency; i++) node->direction[i] = 0;
  }
  if (debug) printf ("FASE 1: orientazione archi aperti (inizio e fine su due (o uno) nodo trivalente\n");
  if (debug) printf ("ce ne sono %d\n", open_arcs);
  if (emb->k > 0)
  {
    connections = (int *) malloc (3*emb->k*sizeof(int));
    count = 0;
    for (inode = 0; inode < emb->k; inode++)
    {
      node = &emb->nodes[inode];
      for (i = 0; i < 3; i++)
      {
        //printf ("Starting from node %d direction %d\n", inode, i);
        assert (node->direction[i] != NODE_IS_START);
        if (node->direction[i] != 0) continue;
        assert (count < open_arcs);
        node->direction[i] = NODE_IS_START;
        node->generator[i] = count;
        nextinode = node->ping[i];
        nextnode = &emb->nodes[nextinode];
        nextnode->generator[node->pong[i]] = count;
        /*
         * now follow the arc untill we reach a trivalent node
         */
        thisnode = node;
        thisi = i;
        while (1)
        {
          nextinode = thisnode->ping[thisi];
          nexti = thisnode->pong[thisi];
          shortcount++;
          nextnode = &emb->nodes[nextinode];
          if (thisnode->valency == 4)
          {
            generator = 2*(thisnode->id - emb->k) + open_arcs;
            if ((thisi % 2) != 0) generator++;
            thisnode->generator[thisi] = generator;
            nextnode->generator[thisnode->pong[thisi]] = generator;
          }
          if (nextnode->valency == 3)
          {
            assert (nextnode->direction[nexti] == 0);
            nextnode->direction[nexti] = NODE_IS_ARRIVAL;
            break;
          }
          nextnode->direction[nexti] = NODE_IS_ARRIVAL;
          thisnode = nextnode;
          thisi = (nexti + 2) % 4;
          assert (nextnode->direction[thisi] == 0);
          nextnode->direction[thisi] = NODE_IS_START;
        }
        connections[3*inode + i] = 3*nextinode + nexti;
        connections[3*nextinode + nexti] = 3*inode + i;
        count++;
      }
    }
    assert (count == open_arcs);
    if (debug)
    {
      printf ("Connections between trivalent nodes:\n");
      for (inode = 0; inode < emb->k; inode++)
      {
        for (i = 0; i < 3; i++)
        {
          connection = connections[3*inode + i];
          printf ("%d.%d -> %d.%d\n", inode, i, connection/3, connection % 3);
        }
      }
    }
    numhcomponents = emb_color (emb, connections);
  }
  numrings = 0;
  if (short_arcs > shortcount)
  {
    if (debug) printf ("FASE 2: orientazione archi chiusi short arcs: %d, processed: %d\n", short_arcs, shortcount);
    /*
     * this look very similar to the 'open arcs' case
     * perhaps we should merge the two cases?
     */
    count = 0;
    for (inode = emb->k; inode < emb->k + emb->n; inode++)
    {
      node = &emb->nodes[inode];
      assert (node->valency == 4);
      for (i = 0; i < 4; i++)
      {
        //if (debug) printf ("Starting from node %d direction %d\n", inode, i);
        if (node->direction[i] != 0) continue;
        count++;
        node->direction[i] = NODE_IS_START;
        /*
         * now follow the arc untill we are back at the starting node
         */
        thisnode = node;
        thisi = i;
        while (1)
        {
          assert (thisnode->valency == 4);
          nextinode = thisnode->ping[thisi];
          nexti = thisnode->pong[thisi];
          shortcount++;
          nextnode = &emb->nodes[nextinode];
          assert (nextnode->valency == 4);
          generator = 2*(thisnode->id - emb->k) + open_arcs;
          if ((thisi % 2) != 0) generator++;
          thisnode->generator[thisi] = generator;
          nextnode->generator[thisnode->pong[thisi]] = generator;
          nextnode->direction[nexti] = NODE_IS_ARRIVAL;
          thisnode = nextnode;
          thisi = (nexti + 2) % 4;
          assert (nextnode->direction[thisi] != NODE_IS_ARRIVAL);
          if (nextnode->direction[thisi] == NODE_IS_START) break;
          nextnode->direction[thisi] = NODE_IS_START;
        }
      }
    }
    assert (short_arcs == shortcount);
    numrings = count;
    if (verbose && numrings) printf ("There %s %d toric component(s)\n", (numrings > 1)?"are":"is", numrings);
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

  assert (numhcomponents + numrings >= 1);
  if (numhcomponents + numrings == 1)
  {
    printf ("numhcomponents: %d - numrings: %d\n", numhcomponents, numrings);
    if (numrings == 1)
    {
      i = emb->k;
      node = &emb->nodes[i];
      j = (node->overpassisodd)?1:0;
      if (node->direction[j] == NODE_IS_ARRIVAL) j = (j + 2) % 4;
      rule = (struct presentationrule *) malloc (sizeof(int) + sizeof (struct presentationrule));
      rule->length = 1;
      rule->var[0] = node->generator[j] + 1;
      rule->next = 0;
      p->elements = rule;

      rule = (struct presentationrule *) malloc (emb->n*sizeof(int) + sizeof(struct presentationrule));
      thisi = i;
      thisj = j;
      rule->length = emb->n;
      k = 0;
      while (1)
      {
        thisnode = &emb->nodes[thisi];
        nexti = thisnode->ping[thisj];
        nextj = thisnode->pong[thisj];
        nextj = (nextj + 2) % 4;
        nextnode = &emb->nodes[nexti];
        if ( ((nextj + nextnode->overpassisodd) % 2) == 1)
        {
          nextjplus = (nextj + 1) % 4;
          sign = 1;
          if (nextnode->direction[nextjplus] == NODE_IS_START) sign = -1;
          rule->var[k] = sign*(nextnode->generator[nextjplus] + 1);
          if (debug) printf ("generator: %d (1 = a; -1 = A)\n", rule->var[k]);
          k++;
        }
        if (nexti == i && nextj == j) break;
        thisi = nexti;
        thisj = nextj;
      }

      rule->next = 0;

      p->elements->next = rule;
    } else {
      emb_meridians_longitudes (emb, connections, p);
      if (!quiet) printf ("Cannot compute meridians and longitudes for links\n");
    }
  }

  if (emb->k > 0) free (connections);

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
readembedding_low (FILE *file)
{
  struct embedding *emb;
  struct emb_node *node, *nodeend, *prevnode, *nodesvec;
  int i, j, tok, count, choice;
  int iend, iendprev, jend;

  //int node_id, node_id2, tnode_id, tnode_pt;

  emb = (struct embedding *) malloc (sizeof (struct embedding));
  emb->choice = emb->k = emb->n = 0;
  emb->nodes = 0;

  tok = gettoken (file);
  if (tok == TOK_COLON)
  {
    tok = gettoken (file);
    assert (tok == ISNUMBER);
    emb->choice = gettokennumber ();
  } else ungettoken (tok);

  tok = gettoken (file);
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
  //assert (gettoken (file) == TOK_SEMICOLON);
  return (emb);
}

/*
 * color the trivalent nodes according to connected components
 * return with the number of colors
 */

int
emb_color (struct embedding *emb, int *connections)
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
          iito = connections[3*ii + jj]/3;
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
emb_meridians_longitudes (struct embedding *emb, int *connections, struct presentation *p)
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
        ii = connections[3*i + k]/3;
        kk = connections[3*i + k] - 3*ii;
        if (node_flood[ii] >= 0)
        {
          // printf ("using: %d.%d -> %d.%d\n", i, k, ii, kk);
          //node_flood[i] = ii;
          node_flood[i] = connections[3*i + k];
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
      underpasses[connections[3*i + k]] = underpasses[3*i + k];


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

      ii = connections[3*i + k]/3;
      llengthpre = numunderpasses_on_spanning_tree (i, node_flood, underpasses);
      llengthpost = numunderpasses_on_spanning_tree (ii, node_flood, underpasses);
      llength = underpasses[3*i + k] + llengthpre + llengthpost;
printf ("should build longitude for arc starting at %d.%d (of length %d)\n", i, k, llength);
      /* adding 1 so that we do not have problems in case llength is zero */
      rule = (struct presentationrule *) malloc ((llength+1)*sizeof (int) + sizeof (struct presentationrule));
      rule->length = llength;

/* TODO: for now fake the longitude as aaa... */
for (u = 0; u < llength; u++) rule->var[u] = 1;
printf ("===== ARC =======\n");
      count = underpasses_on_arc (3*i+k, &(rule->var[llengthpre]), emb);

printf ("===== POST =======\n");
      u = llengthpre + count;
      iii = ii;
      while (iii != 0)
      {
        ikparent = node_flood[iii];
        count = underpasses_on_arc (connections[ikparent], &(rule->var[u]), emb);
        u += count;
        iii = ikparent/3;
      }

printf ("===== PRE =======\n");
      iii = i;
      while (iii != 0)
      {
        ikparent = node_flood[iii];
        u = numunderpasses_on_spanning_tree (ikparent/3, node_flood, underpasses);
printf ("Moving towards root %d -> %d.%d --- u = %d\n", iii, ikparent/3, ikparent % 3, u);
        underpasses_on_arc (ikparent, &(rule->var[u]), emb);
        iii = ikparent/3;
      }
      rule->next = p->elements;
      p->elements = rule;

printf ("building meridian for arc starting at %d.%d\n", i, k);
      rule = (struct presentationrule *) malloc (sizeof (int) + sizeof (struct presentationrule));
      rule->length = 1;
      rule->var[0] = node->generator[k] + 1;
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
printf ("  entering underpasses_on_arc, starting node %d.%d\n", i, k);
  while (1)
  {
    inext = node->ping[k];
    knext = node->pong[k];
    knext = (knext + 2) % 4;

    node = &emb->nodes[inext];
printf ("  in underpasses_on_arc, arc after crossing %d.%d, u=%d\n", inext, knext, u);
    if (node->valency == 3) break;
    knextplus = (knext + 1) % 4;
    k = knext;
    if ((node->overpassisodd + knext) % 2 == 0) continue;
    sign = 1;
    if (node->direction[knextplus] == NODE_IS_START) sign = -1;
    var[u++] = sign*(node->generator[knextplus] + 1);
  }
printf ("  exiting underpasses_on_arc, wrote %d chars\n", u);
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

printf ("UNDERPASSES_ON_SPANNING_TREE for node %d; parent: %d.%d --- underpasses to parent: %d\n", i, parent, k, underpasses[3*parent + k]);
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
