#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include "contour.h"
#include "readembedding.h"
#include "parser.h"

struct sketch *
readembedding (FILE *file)
{
  struct embedding *emb;

  emb = readembedding_low (file);

  printf ("NOT YET IMPLEMENTED: k = %d, n = %d, choice = %d\n", emb->k, emb->n, emb->choice);



/*
  while (isspace (*linept++));
  linept--;
  assert (*linept++ == '{');
  node_id2 = 0;
  while (1)
  {   
    if (*linept == ',') linept++;
    while (isspace (*linept++));
    linept--;
    if (*linept == '}') break;
    node_id = strtol (linept, &linept, 10);
    assert (node_id == node_id2);
    //printf ("node_id: %d\n", node_id);
    assert (*linept++ == ':');
    while (isspace (*linept++)){}; linept--;
    assert (*linept++ == '(');
    tnode_pt = 0;
    while (1)
    {
      tnode_id = strtol (linept, &linept, 10);
      //printf ("tnode_id: %d (%d,%d)\n", tnode_id, node_id2, tnode_pt);
      adjacency[node_id2][tnode_pt] = tnode_id;
      while (isspace (*linept++)){}; linept--;
      tnode_pt++;
      assert (tnode_pt < MAXADJ);
      if (*linept == ')') break;
      assert (*linept++ == ',');
      while (isspace (*linept++)){}; linept--;
    }
    adjacency[node_id2][tnode_pt] = adjacency[node_id2][0];
    adjacencynum[node_id2] = tnode_pt;
    node_id2++;
    assert (node_id2 < MAXNODENUM);
    linept++;
  }
  nodenum = node_id2;
  assert (*linept++ == '}');
  assert (*linept == 0);

  return (nodenum);
 */






  return (0);
}

struct embedding *
readembedding_low (FILE *file)
{
  struct embedding *emb;
  struct emb_node *node;
  int i, tok;

  //int node_id, node_id2, tnode_id, tnode_pt;

  emb = (struct embedding *) malloc (sizeof (struct embedding));
  emb->choice = emb->k = emb->n = 0;
  emb->nodes = 0;
  emb->crossings = 0;

  tok = gettoken (file);
  if (tok == TOK_COLON)
  {
    tok = gettoken (file);
    assert (tok == ISNUMBER);
    emb->choice = gettokennumber ();
  } else ungettoken (tok);

  tok = gettoken (file);
  assert (tok == TOK_LBRACE);
  while (1)
  {
    tok = gettoken (file);
    if (tok == TOK_RBRACE) break;
    assert (tok == ISNUMBER);
    node = (struct emb_node *) malloc (sizeof (struct emb_node));
    node->id = gettokennumber ();
    node->valency = 0;
    tok = gettoken (file);
    assert (tok == TOK_COLON);
    assert (gettoken (file) == TOK_LPAREN);

    for (i = 0; i < 4; i++)
    {
      assert (gettoken (file) == ISNUMBER);
      node->adj[i] = gettokennumber ();
      node->valency++;
      tok = gettoken (file);
      if (tok == TOK_RPAREN) break;
      assert (tok == TOK_COMMA);
    }

    assert (node->valency == 3 || node->valency == 4);
    if (node->valency == 3)
    {
      emb->k++;
      node->next = emb->nodes;
      emb->nodes = node;
    }
    if (node->valency == 4)
    {
      emb->n++;
      node->next = emb->crossings;
      emb->crossings = node;
    }
    tok = gettoken (file);
    if (tok == TOK_RBRACE) break;
    assert (tok == TOK_COMMA);
  }

  assert (gettoken (file) == TOK_SEMICOLON);

  return (emb);
}
