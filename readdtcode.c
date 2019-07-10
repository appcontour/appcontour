/*
 * CREDITS
 *
 * the function "dt_realize" is taken almost verbatim from the KnotScape
 * package (src/decode_new_DT.c) where it is called "realize".
 * Credits for this routine must go to Jeff Weeks, here is a short
 * note about this routine from Morwen Thistlethwaite:
 >
 > The code was written by Jeff Weeks, and was adapted from pseudocode
 > in:  C.H. Dowker and M.B. Thistlethwaite, On the classification of
 > knots,  C.R. Math. Rep. Acad. Sci. Canada, IV (1982) no.2, 129-131.
 >
 > He wrote it so that he could input the knots into Snappea to work
 > out the canonical cell decompositions of their complements.
 */

#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "contour.h"
#include "parser.h"
#include "readdtcode.h"

#define MAXFILELENGTH 2000

extern int debug;
extern int verbose;
extern int quiet;

/*
 * the code vector is twice as the number of nodes, we allow
 * for each odd/even label
 */

static char *rolfsen_to_dt[] = {
#include "rolfsen_to_dt.h"
};

/*
 * a few preliminary definitions.  See at end for the main
 * function, based on Jim Hoste code
 */

extern int experimental;
static int link_components = 1;
static int *start_component = 0;
static int *ori_component = 0;

//typedef unsigned char   Boolean;

static void dt_realize (int *anInvolution, int *aRealization, int aNumCrossings);
static char examplesfilename[MAXFILELENGTH+1];

void uFatalError (char *function, char *file)
{ printf ("FATAL error in function %s, file %s\n", function, file); exit (3); }

struct sketch *
readdtcode (FILE *file)
{
  struct vecofintlist *loiv, *gaussloiv;
  struct sketch *sketch;

  loiv = readvecofintlist (file, LOIV_ISDTCODE);
  assert (loiv->type == LOIV_ISDTCODE || loiv->type == LOIV_ISRDTCODE);

  if (loiv->next)
  {
    /* in this case we need the gauss code! */
    gaussloiv = dtcode2gausscode (loiv);
    assert (gaussloiv);
    freeloiv (loiv);
    sketch = readgausscodeloiv (gaussloiv);  /* this function frees the gaussloiv structure! */
    return (sketch);
  }

  sketch = realize_dtcode (loiv->len, loiv->vec, loiv->handedness);
  freeloiv (loiv);
  return (sketch);
}

/*
 * revert overpass/underpass in a gauss-code for a given nodeid
 */

void
chg_underpass (struct vecofintlist *loiv, int nodenum)
{
  int i, found = 0;
  struct vecofintlist *lv;

  for (lv = loiv; lv; lv = lv->next)
  {
    for (i = 0; i < lv->len; i++)
    {
      if (abs(lv->vec[i]) == nodenum)
      {
        found++;
        lv->vec[i] *= -1;
      }
    }
  }
  assert (found == 2);
}

struct sketch *
readgausscodeloiv (struct vecofintlist *loiv)
{
  struct sketch *sketch;
  struct vecofintlist *lv, *newloiv;
  int i;
  int *dtcode, *dt_involution, *dt_realization;

  /*
   * read gauss code
   */

  assert (loiv->type == LOIV_ISGAUSSCODE);
  if (loiv->next)
  {
    /* case of link */
    newloiv = gausscode_link_to_knot (loiv);
    if (newloiv == 0) return (0);

    dtcode = (int *) malloc (newloiv->len/2 * sizeof (int));
    gauss2dt_knot (newloiv, dtcode, 0);
    dt_involution = (int *) malloc (newloiv->len * sizeof (int));
    dt_realization = (int *) malloc (newloiv->len * sizeof (int));
    for (i = 0; i < newloiv->len/2; i++)
    {
      dt_involution[2*i] = abs(dtcode[i]) - 1;
      dt_involution[abs(dtcode[i]) - 1] = 2*i;
    }
    dt_realize (dt_involution, dt_realization, newloiv->len/2);

    /*
     * at this point dt_realization[i] contains the handedness
     * of newloin->vec[i]
     */

    inherit_gauss2gauss (newloiv, loiv, dt_realization);
    freeloiv (newloiv);
    free (dt_involution);
    free (dt_realization);

    if (verbose)
    {
      printf ("Realized gauss code: ");
      for (lv = loiv; lv; lv = lv->next)
      {
        printf ("{");
        for (i = 0; i < lv->len; i++) printf ("%d%c ", lv->vec[i], (lv->handedness[i]>0)?'<':'>');
        printf ("}");
      }
      printf ("\n");
    }

    sketch = orientedgauss2sketch (loiv);

    free (dtcode);
    freeloiv (loiv);
    return (sketch);
  }
  assert (loiv->next == 0);  /* we do not do links for the moment */

  gauss2dt_knot (loiv, loiv->vec, loiv->handedness);

  sketch = realize_dtcode (loiv->len/2, loiv->vec, loiv->handedness);
  freeloiv (loiv);
  return (sketch);
}

/*
 * add a few new nodes and make surgeries in order to transform a k-components link
 * into a single knot
 */

struct vecofintlist *
gausscode_link_to_knot (struct vecofintlist *loiv)
{
  struct vecofintlist *lv, *loiv2, *newloiv, *newnewloiv;
  int commonnode, maxint;
  int i, j, js;

  assert (loiv && loiv->next);
  assert (loiv->type == LOIV_ISGAUSSCODE);
  maxint = 0;
  /* we need to find a component that is linked to the first one */
  for (lv = loiv; lv->next; lv = lv->next)
  {
    commonnode = -1;
    for (i = 0; i < loiv->len; i++)
    {
      for (j = 0; j < lv->next->len; j++)
      {
        if (abs(loiv->vec[i]) == abs(lv->next->vec[j]))
        {
          commonnode = abs(loiv->vec[i]);
          break;
        }
      }
      if (commonnode >= 0) break;
    }
    if (commonnode >= 0)
    {
      loiv2 = lv->next;
      lv->next = loiv2->next; /* remove lv->next from the list */
      loiv2->next = loiv->next; /* and place it in second position */
      loiv->next = loiv2;
      break;
    }
  }
  if (commonnode < 0)
  {
    printf ("Cannot find common crossing between the first components and the others\n");
    freeloiv (loiv);
    return (0);
  }
  for (lv = loiv; lv; lv = lv->next)
  {
    for (i = 0; i < lv->len; i++) if (abs(lv->vec[i]) > maxint) maxint = abs(lv->vec[i]);
  }

  /* creating new unique gauss code with one more node */
  newloiv = (struct vecofintlist *) malloc (SIZEOFLOIV(2 + loiv->len + loiv2->len));
  newloiv->type = LOIV_ISGAUSSCODE;
  /* fill the new unique gauss code */
  newloiv->next = 0;
  newloiv->dim = newloiv->len = 2 + loiv->len + loiv2->len;
  newloiv->handedness = 0;
  i = 0;
  newloiv->vec[i++] = maxint + 1;  /* this is the new node */
  /* find node in second component */
  for (js = 0; js < loiv2->len; js++) if (abs(loiv2->vec[js]) == commonnode) break;
  assert (js < loiv2->len);
  for (j = js; j < loiv2->len; j++) newloiv->vec[i++] = loiv2->vec[j];
  for (j = 0; j < js; j++) newloiv->vec[i++] = loiv2->vec[j];
  newloiv->vec[i++] = -(maxint + 1);
  /* find node in first component */
  for (js = 0; js < loiv->len; js++) if (abs(loiv->vec[js]) == commonnode) break;
  assert (js < loiv->len);
  for (j = js; j < loiv->len; j++) newloiv->vec[i++] = loiv->vec[j];
  for (j = 0; j < js; j++) newloiv->vec[i++] = loiv->vec[j];
  assert (i == 2 + loiv->len + loiv2->len);
  if (loiv2->next == 0) return (newloiv); /* finished! */

  newloiv->next = loiv2->next;  /* temporarily join the remaining components */
  newnewloiv = gausscode_link_to_knot (newloiv); /* recursive call */
  loiv2->next = newloiv->next;  /* pointers might have changed */

  newloiv->next = 0;
  free (newloiv);

  return (newnewloiv);
}

/*
 *
 */

void
inherit_gauss2gauss (struct vecofintlist *loiv_knot, struct vecofintlist *loiv_link, int *dt_realization)
{
  struct vecofintlist *lv;
  int i, j, newori;
  int ori = 0;

  for (lv = loiv_link; lv; lv = lv->next)
  {
    if (lv->handedness == 0)
    {
      lv->handedness = (int *) malloc (lv->len * sizeof(int));
      for (i = 0; i < lv->len; i++) lv->handedness[i] = 0;
    }
    for (i = 0; i < lv->len; i++)
    {
      for (j = 0; j < loiv_knot->len; j++)
      if (loiv_knot->vec[j] == lv->vec[i])
      {
        newori = lv->handedness[i]*dt_realization[j];
        if (newori)
        {
          if (ori) assert (ori == newori);
          ori = newori;
        }
        lv->handedness[i] = dt_realization[j];
        break;
      }
    }
  }
  switch (ori)
  {
    case 0:
      if (!quiet)
      {
        start_comment ();
        printf ("Warning: knot/link orientation is undefined by its gausss code!\n");
      }
      break;

    case 1:
      break;

    case -1:
      /* need to revert orientation */
      for (lv = loiv_link; lv; lv = lv->next)
      {
        for (i = 0; i < lv->len; i++) lv->handedness[i] = -lv->handedness[i];
      }
      break;

    default:
      printf ("FATAL: ori: %d must be one of 0, -1, 1\n", ori);
      exit (111);
      break;
  }

  return;
}

/*
 * possibilities:
 * [a b c]
 * {a b c}
 * {[a b c]...}
 * {{a b c}...}
 */

struct vecofintlist *
readvecofintlist (FILE *file, int type)
{
  int tok, tok2;
  struct vecofintlist *loiv;
  int dressed = 1;

  tok = gettoken (file);
  if (tok != TOK_LBRACE && tok != TOK_LBRACKET)
  {
    fprintf (stderr, "Error: left brace or bracket expected\n");
    return (0);
  }
  tok2 = gettoken (file);
  if (tok2 != TOK_LBRACE && tok2 != TOK_LBRACKET) dressed = 0;
  ungettoken (tok2);

  if (dressed && tok != TOK_LBRACE)
  {
    fprintf (stderr, "Error: left brace expected\n");
    return (0);
  }
  if (dressed == 0) ungettoken (tok);

  loiv = readnakedvecofintlist (file, type);

  if (dressed)
  {
    tok = gettoken (file);
    if (tok != TOK_RBRACE)
    {
      printf ("Expected terminating }\n");
      freeloiv (loiv);
      return (0);
    }
  }
  return (loiv);
}

/*
 * vector "handedness" is set to +1 for labels of the dtcode postfixed by '>' and
 * to -1 for label postfixed by '<'
 *
 * possibilities: [a b c] or {a b c}
 *                [a b c][d e f]...
 *                {a b c}{d e f}...
 */

struct vecofintlist *
readnakedvecofintlist (FILE *file, int type)
{
  struct vecofintlist *loiv;
  int i, j, tok;
  int startwithlbracket = 0;

  tok = gettoken (file);
  loiv = (struct vecofintlist *) malloc (SIZEOFLOIV(MAXDTCODELEN));
  /*
   * Gauss-code: {-1 10 -11 16 -5 15 -7 12 -9 13 -3 14} {1 -2 3 -4 5 -6 7 -8 9 -10 11 -12 8 -13 2 -14 4 -15 6 -16}
   * DT-code: [4 6 2]
   */
  loiv->type = type;
  assert (tok == TOK_LBRACKET || tok == TOK_LBRACE);
  if (tok == TOK_LBRACKET) startwithlbracket = 1;
  loiv->dim = MAXDTCODELEN;
  loiv->next = 0;
  loiv->handedness = 0;
  i = 0;
  while ((tok = gettoken (file)) == ISNUMBER || tok == TOK_MINUS)
  {
    if (tok == TOK_MINUS)
    {
      tok = gettoken (file);
      assert (tok == ISNUMBER);
      loiv->vec[i] = - gettokennumber ();
    } else loiv->vec[i] = gettokennumber ();
    if (loiv->handedness) loiv->handedness[i] = 0;
    tok = gettoken (file);
    if (tok == KEY_GT || tok == KEY_LT)
    {
      if (loiv->handedness == 0)
      {
        loiv->handedness = (int *) malloc (loiv->dim*sizeof(int));
        for (j = 0; j <= i; j++) loiv->handedness[j] = 0;
      }
      /* handedness is given by the user */
      if (tok == KEY_GT) loiv->handedness[i] = 1;
      else loiv->handedness[i] = -1;
    } else ungettoken (tok);
    i++;
    if (i >= MAXDTCODELEN)
    {
      printf ("Error: dtcode exceeds maximum allowed length: %d\n", MAXDTCODELEN);
      freeloiv (loiv);
      return (0);
    }
    if ((tok = gettoken (file)) != TOK_COMMA)
    {
      ungettoken (tok);  /* comma separation is optional */
    }
  }
  loiv->len = i;
  assert (tok == TOK_RBRACKET || tok == TOK_RBRACE);
  if (startwithlbracket) assert (tok == TOK_RBRACKET);
  tok = gettoken (file);
  if (tok != TOK_COMMA) ungettoken (tok);
  tok = gettoken (file);
  if (tok == TOK_LBRACE || tok == TOK_LBRACKET)
  {
    ungettoken (tok);
    loiv->next = readnakedvecofintlist (file, type);
    if (loiv->next == 0)
    {
      freeloiv (loiv);
      printf ("Syntax error in second component of code\n");
      return (0);
    }
  } else ungettoken (tok);
  return (loiv);
}

void
freeloiv (struct vecofintlist *loiv)
{
  if (loiv->next) freeloiv (loiv->next);
  if (loiv->handedness) free (loiv->handedness);
  free (loiv);
}

/*
 * convert from dtcode to gausscode
 */

struct vecofintlist *
dtcode2gausscode (struct vecofintlist *loiv)
{
  int i, j;
  int numnodes;
  int node, evenlabel, oddlabel;
  struct vecofintlist *newloiv, *lv, *nlv, *prevnlv;
  int *dt_involution, *dt_negative, *dt_realization;

  assert (loiv->type == LOIV_ISDTCODE || loiv->type == LOIV_ISRDTCODE);

  numnodes = 0;
  prevnlv = 0;
  for (lv = loiv; lv; lv = lv->next)
  {
    numnodes += lv->len;
    nlv = (struct vecofintlist *) malloc (SIZEOFLOIV (2*lv->len));
    nlv->next = 0;
    nlv->len = nlv->dim = 2*lv->len;
    nlv->handedness = 0;
    for (i = 0; i < nlv->len; i++) nlv->vec[i] = 0;
    if (prevnlv)
    {
      prevnlv->next = nlv;
    } else {
      newloiv = nlv;
      newloiv->type = LOIV_ISGAUSSCODE;
    }
    prevnlv = nlv;
  }
  dt_involution = (int *) malloc (2*numnodes*sizeof(int));
  dt_negative = (int *) malloc (2*numnodes*sizeof(int));
  dt_realization = (int *) malloc (2*numnodes*sizeof(int));
  for (i = 0; i < 2*numnodes; i++) dt_realization[i] = 0;
  oddlabel = 0;  // beware: in the printed dtcode this is odd, here we have a difference of one
  for (lv = loiv, nlv = newloiv; lv; lv = lv->next, nlv = nlv->next)
  {
    for (i = 0; i < lv->len; i++)
    {
      evenlabel = abs(lv->vec[i]) - 1;
      dt_involution[oddlabel] = evenlabel;
      dt_involution[evenlabel] = oddlabel;
      dt_negative[oddlabel] = 0;
      if (lv->vec[i] < 0) dt_negative[oddlabel] = 1;
      dt_negative[evenlabel] = 1-dt_negative[oddlabel];
      if (lv->handedness)
      {
        dt_realization[oddlabel] = -lv->handedness[i];   //TODO controllare il segno
        dt_realization[evenlabel] = -dt_realization[oddlabel];
      }
      oddlabel += 2;
    }
  }
  evenlabel = 1;
  node = 1;
  for (lv = loiv, nlv = newloiv; lv; lv = lv->next, nlv = nlv->next)
  {
    for (i = 0; i < lv->len; i++)
    {
      nlv->vec[2*i] = node++;
      if (dt_negative[evenlabel-1]) nlv->vec[2*i] = -nlv->vec[2*i];
      nlv->vec[2*i+1] = dt_involution[evenlabel]/2 + 1;
      if (dt_negative[evenlabel]) nlv->vec[2*i+1] = -nlv->vec[2*i+1];
      if (dt_realization[evenlabel-1] || dt_realization[evenlabel])
      {
        if (nlv->handedness == 0)
        {
          nlv->handedness = (int *) malloc (nlv->len*sizeof(int));
          for (j = 0; j < nlv->len; j++) nlv->handedness[j] = 0;
        }
        nlv->handedness[2*i] = dt_realization[evenlabel-1];
        nlv->handedness[2*i+1] = dt_realization[evenlabel];
      }
      evenlabel += 2;
    }
  }

  free (dt_involution);
  free (dt_negative);
  free (dt_realization);
  return (newloiv);
}

/*
 * read gauss code and convert it into a dtcode
 */

struct vecofintlist *
gausscode2dtcode (struct vecofintlist *loiv)
{
  int numnodes, i, j, k, found;
  struct vecofintlist *newloiv, *lv, *nlv, *prevnlv;
  int *node_label, *dt_involution, *dt_sign, *dt_realization;

  assert (loiv->type == LOIV_ISGAUSSCODE || loiv->type == LOIV_ISRGAUSSCODE);
  if (loiv->next == 0)
  {
    /* special case of knot */
    newloiv = (struct vecofintlist *) malloc (SIZEOFLOIV (loiv->len/2));
    newloiv->type = LOIV_ISDTCODE;
    newloiv->len = newloiv->dim = loiv->len/2;
    newloiv->next = 0;
    newloiv->handedness = 0;
    if (loiv->handedness) newloiv->handedness = (int *) malloc (loiv->len/2 * sizeof(int));
    gauss2dt_knot (loiv, newloiv->vec, newloiv->handedness);

    return (newloiv);
  }

  if (check_gauss_compat(loiv) == 0)
  {
    fprintf (stderr, "Gauss code has starting points incompatible with corresponing DTcode\n");
    return (0);
  }

  numnodes = 0;
  prevnlv = 0;
  for (lv = loiv; lv; lv = lv->next)
  {
    assert ((lv->len % 2) == 0);
    numnodes += lv->len/2;
    nlv = (struct vecofintlist *) malloc (SIZEOFLOIV (lv->len/2));
    nlv->next = 0;
    nlv->len = nlv->dim = lv->len/2;
    nlv->handedness = 0;
    for (i = 0; i < nlv->len; i++) nlv->vec[i] = 0;
    if (prevnlv)
    {
      prevnlv->next = nlv;
    } else {
      newloiv = nlv;
      newloiv->type = LOIV_ISDTCODE;
    }
    prevnlv = nlv;
  }

  node_label = (int *) malloc (numnodes*sizeof(int));
  dt_involution = (int *) malloc (2*numnodes*sizeof(int));
  dt_sign = (int *) malloc (2*numnodes*sizeof(int));
  dt_realization = (int *) malloc (2*numnodes*sizeof(int));
  for (i = 0; i < 2*numnodes; i++) dt_realization[i] = 0;

  /* table of node names by looking at odd positions */
  for (lv = loiv, j = 0; lv; lv = lv->next)
  {
    for (i = 0; i < lv->len; i += 2)
    {
      dt_sign[2*j] = 1;
      if (lv->vec[i] < 0) dt_sign[2*j] = -1;
      if (lv->handedness)
      {
        dt_realization[2*j] = -lv->handedness[i];
        dt_realization[2*j+1] = -lv->handedness[i+1];
      }
      node_label[j++] = abs(lv->vec[i]);
    }
  }

  for (lv = loiv, j = 0; lv; lv = lv->next)
  {
    for (i = 1; i < lv->len; i += 2)
    {
      /* search corresponding node */
      found = 0;
      for (k = 0; k < numnodes; k++)
      {
        if (abs(lv->vec[i]) == node_label[k])
        {
          found = 1;
          dt_involution[2*j + 1] = 2*k;
          dt_involution[2*k] = 2*j + 1;
          break;
        }
      }
      j++;
      assert (found == 1);
    }
  }

  for (i = 0; i < numnodes; i++)
  {
    dt_sign[dt_involution[2*i]] = - dt_sign[2*i];
    if (dt_realization[dt_involution[2*i]]) dt_realization[2*i] = - dt_realization[dt_involution[2*i]];
    //if (dt_realization[2*i]) dt_realization[dt_involution[2*i]] = - dt_realization[2*i];
  }

  for (nlv = newloiv, j = 0; nlv; nlv = nlv->next)
  {
    for (i = 0; i < nlv->len; i++)
    {
      nlv->vec[i] = dt_sign[j]*(dt_involution[j] + 1);
      if (dt_realization[j])
      {
        if (nlv->handedness == 0)
        {
          nlv->handedness = (int *) malloc (nlv->len*sizeof(int));
          for (k = 0; k < nlv->len; k++) nlv->handedness[k] = 0;
        }
        nlv->handedness[i] = dt_realization[j];
      }
      j += 2;
    }
  }

  free (node_label);
  free (dt_involution);
  free (dt_sign);
  free (dt_realization);
  return (newloiv);
}

void
gauss2dt_knot (struct vecofintlist *loiv, int *vecofint, int *handedness)
{
  int *nodeinfo, *nodeinfo2;
  int i, j, isodd, dtcode, nodename;
  int nlabels = loiv->len;
  int nnodes = nlabels/2;

  assert (nlabels == 2*nnodes);
  assert (loiv->next == 0);
  if (loiv->handedness) assert (handedness);

  nodeinfo = (int *) malloc (nnodes*sizeof(int));
  nodeinfo2 = (int *) malloc (nnodes*sizeof(int));

  for (i = 0; i < nnodes; i++) nodeinfo[i] = nodeinfo2[i] = 0;

  for (i = 0; i < nlabels; i++)
  {
    nodename = abs(loiv->vec[i]);
    assert (nodename >= 1 && nodename <= nnodes);
    if (nodeinfo[nodename - 1] == 0)
    {
      nodeinfo[nodename - 1] = i+1;
      if (loiv->vec[i] < 0) nodeinfo[nodename - 1] = -i-1;
      continue;
    }
    assert (nodeinfo2[nodename - 1] == 0);
    if (loiv->vec[i] > 0)
    {
      nodeinfo2[nodename - 1] = i+1;
      assert (nodeinfo[nodename - 1] < 0);
    } else {
      nodeinfo2[nodename - 1] = -i-1;
      assert (nodeinfo[nodename - 1] > 0);
    }
  }
  for (i = 0, j = 0, isodd = 0; i < nlabels; i++)
  {
    isodd = 1 - isodd;
    if (! isodd) continue;
    nodename = abs(loiv->vec[i]);
    dtcode = -nodeinfo[nodename - 1];
    if (abs(nodeinfo[nodename - 1]) == i + 1) dtcode = -nodeinfo2[nodename - 1];
    if (loiv->handedness) handedness[j] = -loiv->handedness[i];
    vecofint[j++] = dtcode;
  }

  free (nodeinfo);
  free (nodeinfo2);

  return;
}

/*
 * check if gauss code (of link) has starting points compatible with the
 * coversion into dtcode (i.e., lead to an odd-even pairint)
 */

int
check_gauss_compat (struct vecofintlist *loiv)
{
  int i, j;
  struct vecofintlist *lv1, *lv2;

  for (lv1 = loiv; lv1; lv1 = lv1->next)
  {
    for (lv2 = lv1; lv2; lv2 = lv2->next)
    {
      if (lv1 == lv2) continue;
      for (i = 0; i < lv1->len; i++)
      {
        for (j = 0; j < lv2->len; j++)
        {
          if (abs(lv1->vec[i]) == abs(lv2->vec[j]))
          {
            if (((i + j) % 2) != 1) return (0);
            if (lv1->vec[i] != -lv2->vec[j]) return (0);
          }
        }
      }
    }
  }
  return (1);
}

/*
 * just "realize" an existing dtcode with partial orientation
 * loiv contains a dtcode
 */

void
realize_loiv (struct vecofintlist *loiv)
{
  int i, *gregionsign;

  if (loiv->type == LOIV_ISRDTCODE) return;
  assert (loiv->type == LOIV_ISDTCODE);
  gregionsign = loiv->handedness;
  if (loiv->handedness == 0)
  {
    gregionsign = loiv->handedness = (int *) malloc (loiv->len*sizeof(int));
    for (i = 0; i < loiv->len; i++) gregionsign[i] = 0;
  }

  realize_loiv_split (loiv->len, loiv->vec, gregionsign);
}

/*
 * main function that reconstructs the correct crossings handedness
 */

static int numnodes;
static int numlabels;
static int numregions;
static int *dt_involution;
static int *dt_realization;
static int *tagged;
static int *marknodes;

void
realize_loiv_split (int lnumnodes, int *vecofint, int *gregionsign)
{
  extern struct global_data globals;
  int i, numconsistent;
  int dt_incomplete;
  int agree, agreeneg;
  char ch;
  int sign, label;

  numnodes = lnumnodes;
  numlabels = 2*numnodes;
  numregions = 2 + numlabels - numnodes;

  assert (numnodes >= 2);

  dt_involution = (int *) malloc (numlabels*(sizeof (int)));
  dt_realization = (int *) malloc (numlabels*(sizeof (int)));
  tagged = (int *) malloc (numlabels * sizeof(int) );
  marknodes = (int *) malloc (numlabels * sizeof(int) );

  for (i = 0; i < numnodes; i++)
  {
    label = abs(vecofint[i]) - 1;
    sign = 1;
    if (vecofint[i] < 0) sign = -1;
    if ((label % 2) != 1)
    {
      printf ("Only even values allowed!\n");
      exit (4);
    }
    dt_involution[2*i] = label;
    dt_involution[label] = 2*i;
  }

  for (i = 0; i < numlabels; i++)
  {
    if (dt_involution[i] < 0 || dt_involution[i] >= numlabels)
    {
      printf ("Must use all even number from 2 to %d\n", numlabels);
      exit (5);
    }
    if (((dt_involution[i] + i) % 2) != 1)
    {
      printf ("Involution is not an even-odd coupling\n");
      exit (6);
    }
    if (i != dt_involution[dt_involution[i]])
    {
      printf ("This is not an involution\n");
      exit (7);
    }
    if (abs(i - dt_involution[i]) <= 1)
    {
      printf ("No tight loop allowed: %d %d\n", i + 1, dt_involution[i] + 1);
      exit (8);
    }
  }

  /* this function wants a zero-based vector */
  dt_realize (dt_involution, dt_realization, numnodes);
  dt_incomplete = 0;
  for (i = 0; i < numlabels; i++) if (dt_realization[i] == 0) dt_incomplete = 1;

  /*
   * this (partial?) reconstruction must be consistent with possible informations given
   * by the user
   */

  agree = 1;
  agreeneg = 1;
  if (gregionsign[0] == 0 && globals.dtcode_realize) gregionsign[0] = globals.dtcode_realize;
  for (i = 0; i < numnodes; i++)
  {
    if (gregionsign[i]*dt_realization[2*i] > 0) agreeneg = 0;
    if (gregionsign[i]*dt_realization[2*i] < 0) agree = 0;
  }
  if ((agree == 0) && (agreeneg == 0))
  {
    printf ("Fatal: the reconstruction is inconsistent with user-supplied values.\n");
    if (verbose)
    {
      for (i = 0; i < numnodes; i++)
      {
        sign = 1;
        if (vecofint[i] < 0) sign = -1;
        printf ("%d%c ", sign*(dt_involution[2*i]+1), (dt_realization[2*i] > 0)?'>':'<');
      }
      printf ("\n");
    }
    exit (4);
  }
  if (agree && agreeneg && (!quiet)
       //&& 0    // TODO: disable it for now, just to pass all tests!
       )
  { /* orientation is not forced by the user, taking the one given by dt_realize */
    start_comment ();
    printf ("Warning: knot orientation is undefined by its dtcode!\n");
    start_comment (); printf ("using the realization given by dt_realize: ");
    printf ("at first crossing (labelled 1)\n");
    start_comment (); printf ("the intersecting path is oriented right-to-left\n");
    start_comment (); printf ("you can force this, and avoid this warning, with the option --right\n");
    start_comment (); printf ("or reverse the orientation with the option --left\n");
  }
  if (agree == 0)
  {
    for (i = 0; i < numlabels; i++) dt_realization[i] = -dt_realization[i];
  }
  for (i = 0; i < numnodes; i++)
    if (dt_realization[2*i]) gregionsign[i] = dt_realization[2*i];
  if (verbose)
  {
    for (i = 0; i < numnodes; i++)
    {
      ch = '?';
      if (dt_realization[2*i] > 0) ch = '>';
      if (dt_realization[2*i] < 0) ch = '<';
      sign = 1;
      if (vecofint[i] < 0) sign = -1;
      printf ("%d%c ", sign*(dt_involution[2*i]+1), ch);
    }
    printf ("\n");
  }

  if (dt_incomplete)
  {
    /* fallback to my old code */
    numconsistent = reconstruct_sign (0, gregionsign);
    if (numconsistent <= 0)
    {
      printf ("No consistent completions found. Perhaps this is a composite knot\n");
      exit (3);
    }
    if (numconsistent > 1)
    {
      printf ("More than one (%d) consistent completion found!\n", numconsistent);
      printf ("maybe this is a compound knot\n");
      exit (2);
    }
    reconstruct_sign (1, gregionsign);
  }
  for (i = 0; i < numnodes; i++)
    if (dt_realization[2*i]) gregionsign[i] = dt_realization[2*i];

  free (dt_realization);
  free (dt_involution);
  free (tagged);
  free (marknodes);
}

static int *choices;
static int choicept;

int maximal_expansion (void);
void make_choice (void);

int
reconstruct_sign (int which, int *gregionsign)
{
  int i, numtagged, countreconstructions;

  choices = (int *) malloc ( numnodes*sizeof(int) );
  choices[0] = 0;
  choicept = 0;
  countreconstructions = 0;

  while (1) /* cycle on all possible reconstructions */
  {
    /* make all regionsigns unknown */
    numtagged = 0;
    for (i = 0; i < numlabels; i++) dt_realization[i] = 0;
    if (gregionsign)
    {
      for (i = 0; i < numnodes; i++)
      {
        if (gregionsign[i] && (dt_realization[2*i] == 0))
        {
          dt_realization[2*i] = gregionsign[i];
          dt_realization[dt_involution[2*i]] = -gregionsign[i];
          numtagged++;
          numtagged += maximal_expansion ();
        }
      }
    }
    while (numtagged < numnodes) /* subsequent waves of reconstruction */
    {
      make_choice ();
      numtagged++;
      numtagged += maximal_expansion ();
      if (debug) printf ("reconstruction: %d, covered %d of %d\n", countreconstructions, numtagged, numnodes);
    }
    if (isconsistent ())
    {
      countreconstructions++;
      if (countreconstructions == which)
      {
        /* found desired reconstruction */
        free (choices);
        return (countreconstructions);
      }
    }
    /* next possible choice... */
    assert (choices[choicept] == 0);
    if (choicept == 0)
    { /* exhausted all possible choices */
      free (choices);
      return (countreconstructions);
    }
    assert (choicept > 0);
    choicept--;
    while (choices[choicept] == -1)
    {
      choices[choicept] = 0;
      if (choicept == 0)
      { /* exhausted all possible choices */
        free (choices);
        return (countreconstructions);
      }
      choicept--;
    }
    assert (choices[choicept] == 1);
    choices[choicept] = -1;
    choices[choicept+1] = 0;
    choicept = 0;
  }
  free (choices);
  assert (0);
  return (countreconstructions);
}

void
make_choice (void)
{
  int i;

  for (i = 0; i < numlabels; i++)
  {
    if (dt_realization[i] == 0)
    {
      /* found an unoriented node */
      if (choices[choicept] == 0)
      {
        choices[choicept] = 1;
        choices[choicept + 1] = 0;
      }
      assert (abs(choices[choicept]) == 1);
      dt_realization[i] = choices[choicept];
      dt_realization[dt_involution[i]] = -choices[choicept];
      choicept++;
      return;
    }
  }
  assert (0);
  return;
}

int
maximal_expansion ()
{
  int insist = 1;
  int totexpansions = 0;
  int i, j, node, expansions;

  while (insist)
  {
    insist = 0;
    for (i = 0; i < numlabels; i++)
    {
      for (j = 0; j < numlabels; j++) tagged[j] = 0;
      /* follow cicle */
      node = i;
      while (1)
      {
        if (tagged[node])
        {
          // printf ("found cycle starting at node: %d\n", dt_involution[node] + 1);
          expansions = inherit (dt_involution[node]);
          if (expansions) insist = 1;
          totexpansions += expansions;
          // printf ("Expanded by %d nodes\n", expansions);
          break;
        }
        assert (tagged[dt_involution[node]] == 0);
        tagged[dt_involution[node]] = 1;
        node = nextlabel (node);
      }
    }
  }
  return (totexpansions);
}

/*
 * build the apparent contour using the realized dtcode
 */

static int *dtsign;
static int *componentoflabel;

struct sketch *
realize_dtcode (int lnumnodes, int *vecofint, int *gregionsign)
{
  struct sketch *sketch;
  int i;
  int sign, label;
  int freegregionsign = 0;

  if (gregionsign == 0)
  {
    freegregionsign = 1;
    gregionsign = (int *) malloc (lnumnodes*sizeof(int));
    for (i = 0; i < lnumnodes; i++) gregionsign[i] = 0;
  }

  realize_loiv_split (lnumnodes, vecofint, gregionsign);

  numnodes = lnumnodes;
  numlabels = 2*numnodes;
  assert (numnodes >= 2);

  dt_involution = (int *) malloc (numlabels*(sizeof (int)));
  dtsign = (int *) malloc (numlabels*(sizeof (int)));
  dt_realization = (int *) malloc (numlabels*(sizeof (int)));
  tagged = (int *) malloc (numlabels * sizeof(int) );
  marknodes = (int *) malloc (numlabels * sizeof(int) );

  for (i = 0; i < numnodes; i++)
  {
    label = abs(vecofint[i]) - 1;
    if ((label % 2) != 1)
    {
      printf ("Only even values allowed!\n");
      exit (4);
    }
    dt_involution[2*i] = label;
    dt_involution[label] = 2*i;
    sign = 1;
    if (vecofint[i] < 0) sign = -1;
    dtsign[2*i] = sign;
    dtsign[label] = -sign;
  }

  numregions = 2 + numlabels - numnodes;

  for (i = 0; i < numnodes; i++)
  {
    dt_realization[2*i] = gregionsign[i];
    dt_realization[dt_involution[2*i]] = -gregionsign[i];
  }

  sketch = newsketch ();
  display_arcs_from_arcs (sketch);
  display_arcs_from_nodes (sketch);
  display_regions (sketch);
  display_regions_from_arcs (sketch);
  display_regions_from_nodes (sketch);
  if (freegregionsign) free (gregionsign);
  free (tagged);
  free (dt_realization);
  free (dtsign);
  free (dt_involution);
  free (marknodes);
  sketch->huffman_labelling = 1;
  if (sketch->regions->next == 0) {
    fprintf (stderr, "Warning: empty sketch!\n");
    sketch->isempty = sketch->huffman_labelling = 1;
  }
  if (debug) printsketch (sketch);
  postprocesssketch (sketch);
  return (sketch);
}

/*
 *
 */

static int *dteven = 0;
static int *dtodd = 0;
static int *dtnode = 0;

int dt_iseven (int);

struct sketch *
orientedgauss2sketch (struct vecofintlist *loiv)
{
  struct vecofintlist *lv;
  struct sketch *sketch;
  int i, nodenum, sum, goon, nf, labelid, numarcs = 0;
  int component;
  int *nodesfound;

  /* setup data */

  link_components = 0;
  for (lv = loiv; lv; lv = lv->next)
  {
    link_components++;
    numarcs += lv->len;
  }

  assert ( (numarcs % 2) == 0);

  numlabels = numarcs;
  numnodes = numlabels/2;
  numregions = 2 + numlabels - numnodes;

  dt_involution = (int *) malloc (numlabels*(sizeof (int)));
  dtsign = (int *) malloc (numlabels*(sizeof (int)));
  dt_realization = (int *) malloc (numlabels*(sizeof (int)));
  tagged = (int *) malloc (numlabels * sizeof(int) );
  componentoflabel = (int *) malloc (numlabels*(sizeof (int)));
  dtnode = (int *) malloc (numlabels*(sizeof (int)));
  dteven = (int *) malloc (numlabels*(sizeof (int)));
  dtodd = (int *) malloc (numlabels*(sizeof (int)));

  start_component = (int *) malloc (link_components * sizeof(int));
  ori_component = (int *) malloc (link_components *sizeof(int));
  for (i = 0; i < link_components; i++) ori_component[i] = 0;
  ori_component[0] = 1;

  nodesfound = (int *) malloc (numarcs/2 * sizeof (int));
  for (i = 0; i < numarcs/2; i++) nodesfound[i] = -1;

  labelid = 0;
  for (lv = loiv, component = 0; lv; lv = lv->next)
  {
    start_component[component] = labelid;
    for (i = 0; i < lv->len; i++)
    {
      componentoflabel[labelid] = component;
      dtsign[labelid] = 1;
      if (lv->vec[i] < 0) dtsign[labelid] = -1;
      dt_realization[labelid] = lv->handedness[i];
      nodenum = abs(lv->vec[i]) - 1;
      dtnode[labelid] = nodenum;
      assert (nodenum < numnodes);
      nf = nodesfound[nodenum];
      if (nf >= 0)
      {
        dt_involution[labelid] = nf;
        dt_involution[nf] = labelid;
      } else {
        nodesfound[nodenum] = labelid;
      }
      labelid++;
    }
    component++;
  }
  free (nodesfound);
  /* give parity to each link component */
  goon = 1;
  while (goon)
  {
    goon = 0;
    for (i = 0; i < numlabels; i++)
    {
      component = componentoflabel[i];
      if (ori_component[component] != 0) continue;
      /* try to orient this component */
      if (ori_component[componentoflabel[dt_involution[i]]])
      {
        sum = i + dt_involution[i];
        ori_component[component] = ori_component[componentoflabel[dt_involution[i]]];
        if ((sum % 2) == 0)
        {
          ori_component[component] = - ori_component[component];
        }
        goon = 1;
      }
    }
  }
  for (component = 0; component < link_components; component++)
  {
    assert (ori_component[component]);
  }
  for (i = 0; i < numnodes; i++) dteven[i] = dtodd[i] = -1;
  for (i = 0; i < numlabels; i++)
  {
    component = componentoflabel[i];
    nodenum = dtnode[i];
    if (dt_iseven(i)) dteven[nodenum] = i;
      else dtodd[nodenum] = i;
  }
  for (i = 0; i < numnodes; i++) assert (dteven[i] >= 0 && dtodd[i] >= 0);

  sketch = newsketch ();
  display_arcs_from_arcs (sketch);
  display_arcs_from_nodes (sketch);
  display_regions (sketch);
  display_regions_from_arcs (sketch);
  display_regions_from_nodes (sketch);
  free (tagged);
  free (dt_realization);
  free (dtnode);
  free (dteven);
  free (dtodd);
  free (dtsign);
  free (dt_involution);
  free (start_component);
  free (ori_component);
  free (componentoflabel);
  sketch->huffman_labelling = 1;
  if (sketch->regions->next == 0) {
    fprintf (stderr, "Warning: empty sketch!\n");
    sketch->isempty = sketch->huffman_labelling = 1;
  }
  if (debug) printsketch (sketch);
  postprocesssketch (sketch);
  return (sketch);
}

/*
 * Read into the knotscape pak files.
 * Please note: the various *.pak files have been created by Morwen Thistlethwaite and Jim Hoste
 * by use with the well-known software knotscape.
 * They are part of the 'knotscape' program.
 *
 * I would like to thank them for their work!
 *
 * The 
 */

#ifndef PKGDATA_DIR
  #define PKGDATA_DIR ""
#endif

#define MAXFILELENGTH 2000

#define TOKENWORDLENGTH 200

static char tokenword[TOKENWORDLENGTH];
static char pathname[MAXFILELENGTH+1];

FILE *
open_pak_file (char *pakname)
{
  FILE *pakfile;
  char *knotscape_homes[]={PKGDATA_DIR, ".", "/home", "/usr/local", "/usr/local/share", 0};
  char *pakpaths[]={"data", "knotTable", "knotscape/knotTable", "knotscape/knotscape_1.01/knotTable", 0};
  int ip1, ip2;

  for (ip1 = 0; knotscape_homes[ip1]; ip1++)
  {
    for (ip2 = 0; pakpaths[ip2]; ip2++)
    {
      strncpy (pathname, knotscape_homes[ip1], MAXFILELENGTH);
      strncat (pathname, "/", MAXFILELENGTH);
      strncat (pathname, pakpaths[ip2], MAXFILELENGTH);
      strncat (pathname, "/", MAXFILELENGTH);
      strncat (pathname, pakname, MAXFILELENGTH);
      strncat (pathname, ".pak", MAXFILELENGTH);
      if (access (pathname, R_OK)) continue;
      pakfile = fopen (pathname, "r");
      return (pakfile);
    }
  }
  printf ("Cannot find knotscape installation.  You should install the knotscape package by\n");
  printf ("Morwen Thistlethwaite and Jim Hoste in /home/knotscape or /usr/local/knotscake\n");
  printf ("The pak files should reside in a folder with a name like /home/knotscape/knotscape_1.01/knotTable/\n");
  printf ("   or /usr/local/share/appcontour/data/\n");
  exit (2);
}

/*
 * drop initial 'K' if present
 */

void
canonify_knotname (void)
{
  int i;

  if (*tokenword == 'K')
  {
    for (i = 0; tokenword[i+1]; i++)
    {
      tokenword[i] = tokenword[i+1];
    }
    tokenword[i] = 0;
  }
  return;
}

struct vecofintlist *
readknotscape (FILE *file, struct sketch **sketchpt)
{
  extern struct global_data globals;
  extern int quiet, verbose;
  char basename[20];
  int nodeid, i, tok, crossings, codelen, alternate, knotnum;
  int sign, ch, tailstart, high, low;
  int sum1 = 0;
  int sum2 = 0;
  int found = 0;
  char *namept;
  int *dtcode, *dtpositive;
  FILE *pakfile, *infile;
  struct vecofintlist *loiv;

  *sketchpt = 0;
  tok = gettoken (file);
  assert (tok == TOK_LBRACE);

  if (getword (file, tokenword, TOKENWORDLENGTH) == TOK_EOF) return (0);

  if (*tokenword == '_')
  {
    assert (tokenword[1] == '1' || tokenword[1] == '2');
    strncpy (tokenword, (tokenword[1] == '1')?globals.knotname1:globals.knotname2, TOKENWORDLENGTH);
  }

  if (*tokenword == 'L')
  {
    /* case of link */
    tok = gettoken (file);
    assert (tok == TOK_RBRACE);
    return (readlinkfromtable (tokenword));
  }

  if (*tokenword == 'H')
  {
    if (tokenword[1] != 'K')
    {
      fprintf (stderr, "Invalid handlebody-knot name: %s\n", tokenword);
      exit (12);
    }
    strncpy (examplesfilename, EXAMPLES_DIR, MAXFILELENGTH);
    strncat (examplesfilename, "/handlebody_knots/hk", MAXFILELENGTH);
    strncat (examplesfilename, &tokenword[2], MAXFILELENGTH);
    strncat (examplesfilename, ".knot", MAXFILELENGTH);
    infile = fopen (examplesfilename, "r");
    if (infile && quiet == 0) fprintf (stderr, "Reading from file %s\n", examplesfilename);
    *sketchpt = readcontour (infile);
    if (*sketchpt == 0) exit (15);
    return (0);
  }

  canonify_knotname ();
  for (i = 0; rolfsen_to_dt[i]; i += 2)
  {
    if (strcmp (tokenword, rolfsen_to_dt[i]) == 0)
    {
      strcpy (tokenword, rolfsen_to_dt[i+1]);
      if (*tokenword == '#')
      {
        printf ("# Warning: there is a difference in the numbering of the last four knots in Rolfsen table\n");
        printf ("# with respect to other references, e.g. the knot atlas.  To avoid this warning you should\n");
        printf ("# postfix the knot name with 'R' (Rolfsen book), or 'KA' (Knot Atlas): e.g. 10_165KA refers\n");
        printf ("# to the last knot in the Knot Atlas and is the same as 10_166R\n");
        strcpy (tokenword, rolfsen_to_dt[i+1]+1);
      }
      if (!quiet)
      {
        printf ("# Converting Rolfsen notation into dt numbering: %s -> %s\n", rolfsen_to_dt[i], tokenword);
      }
    }
  }

  if (!quiet) printf ("# ['pak' files courtesy of Morwen Thistlethwaite and Jim Hoste, in knotscape package]\n");
  //printf ("Knot name: %s\n", tokenword);

  crossings = strtol (tokenword, &namept, 10);
  assert (crossings >= 3 && crossings <= 16);
  switch (*namept++)
  {
    case 'a':
    alternate = 1;
    break;

    case 'n':
    alternate = 0;
    break;

    default:
    printf ("Invalid character after number of crossings: %c (valid values: 'a' or 'n')\n", *namept);
    exit (2);
    break;
  }
  if (*namept == '_') namept++;
  assert (isdigit(*namept));
  knotnum = strtol (namept, &namept, 10);

  tok = gettoken (file);
  assert (tok == TOK_RBRACE);

  sprintf (basename, "%d%c", crossings, (alternate)?'a':'n');

  loiv = (struct vecofintlist *) malloc (SIZEOFLOIV(crossings));
  loiv->len = crossings;
  loiv->type = LOIV_ISDTCODE;
  loiv->next = 0;
  loiv->handedness = 0;
  dtcode = loiv->vec;

  //dtcode = (int *) malloc (crossings *sizeof (int));
  dtpositive = (int *) malloc (crossings *sizeof (int));

  pakfile = open_pak_file (basename);
  codelen = crossings - 1;
  for (i = 0; i < crossings; i++) {dtcode[i] = 0; dtpositive[i] = 1;}
  nodeid = 0;
  while ((ch = fgetc (pakfile)) != EOF)
  {
    nodeid++;
    high = ch/16;
    low = ch - 16*high;
    assert (high == low);
    tailstart = codelen - high;
    for (i = tailstart; i < codelen;)
    {
      ch = fgetc (pakfile);
      high = ch/16;
      low = ch - 16*high;
      dtcode[i++] = high;
      if (i >= codelen) break;
      dtcode[i++] = low;
    }
    if (alternate == 0)
    {
      ch = fgetc (pakfile);
      for (i = 7; i >= 0; i--)
      {
        if (i < crossings) dtpositive[i] = ch & 1;
        ch = ch >> 1;
      }
      ch = fgetc (pakfile);
      for (i = 15; i >= 8; i--)
      {
        if (i < crossings) dtpositive[i] = ch & 1;
        ch = ch >> 1;
      }
    }
    if (nodeid == knotnum) {found = 1; break;}
  }
  if (!found)
  {
    printf ("Cannot find knot %d in pak file %s\n", knotnum, pathname);
    exit (2);
  }

  for (i = 0; i < codelen; i++)
  {
    sign = 2*dtpositive[i]-1;
    sum1 += i;
    sum2 += dtcode[i];
    dtcode[i] = sign*(2*dtcode[i] + 2);
  }
  assert (i == codelen);
  sum1 += codelen;
  sign = 2*dtpositive[codelen]-1;
  dtcode[i] = sign*(2*(sum1 - sum2) + 2);
  codelen++;
  free (dtpositive);
  if (verbose) printloiv (loiv);
  //sketch = realize_dtcode (codelen, dtcode, 0);
  return (loiv);
}

/*
 * search for the data file for link definition
 */

struct vecofintlist *
readlinkfromtable (char *linkname)
{
  extern struct global_data globals;
  struct vecofintlist *loiv;
  FILE *linkfile;
  char *ch, *gc, line[2000];

  char *data_homes[]={".", PKGDATA_DIR, ".", 0};
  char *data_paths[]={"data", 0};
  int ip1, ip2;

  linkfile = 0;
  for (ip1 = 0; data_homes[ip1]; ip1++)
  {
    for (ip2 = 0; data_paths[ip2]; ip2++)
    {
      strncpy (pathname, data_homes[ip1], MAXFILELENGTH);
      strncat (pathname, "/", MAXFILELENGTH);
      strncat (pathname, data_paths[ip2], MAXFILELENGTH);
      strncat (pathname, "/", MAXFILELENGTH);
      strncat (pathname, "links_gausscodes.txt", MAXFILELENGTH);
      if (access (pathname, R_OK)) continue;
      linkfile = fopen (pathname, "r");
      break;
    }
    if (linkfile) break;
  }

  if (linkfile == 0)
  {
    printf ("Fatal: Cannot open file %s, check installation\n", pathname);
    exit (1);
  }

  while (fgets (line, 2000, linkfile))
  {
    for (ch = line; ch; ch++)
    {
      if (*ch == ':')
      {
        gc = ch+1;
        *ch = 0;
        if (strcmp (line, linkname) == 0)
        {
          loiv = read_gausscode_from_string (gc);
          return (loiv);
        }
        break;
      }
    }
  }

  printf ("Cannot find link name: %s in file %s\n", linkname, pathname);
  fclose (linkfile);
  return (0);
}

struct vecofintlist *
read_gausscode_from_string (char *gc)
{
  struct vecofintlist *loiv;
  char *chpt = gc;
  int sign, i;

  while (isspace (*chpt)) chpt++;
  if (*chpt != '{') return (0);
  chpt++;
  loiv = (struct vecofintlist *) malloc (SIZEOFLOIV(MAXDTCODELEN));
  loiv->type = LOIV_ISGAUSSCODE;
  loiv->dim = MAXDTCODELEN;
  loiv->next = 0;
  loiv->handedness = 0;
  i = 0;
  while (isspace (*chpt)) chpt++;
  while (*chpt != '}')
  {
    sign = 1;
    if (*chpt == '-') {sign = -1; chpt++;}
    if (*chpt == '+') {sign = 1; chpt++;}
    assert (i < MAXDTCODELEN);
    loiv->vec[i] = strtol (chpt, &chpt, 10);
    loiv->vec[i] = sign*loiv->vec[i];
    i++;
    while (isspace (*chpt)) chpt++;
  }
  loiv->len = i;
  chpt++;
  loiv->next = read_gausscode_from_string (chpt);

  return (loiv);
}

/*
 * node arcs are numbered from 1 to numlabels, each produces two
 * apparent contour arcs, one on the right, same orientation, i -> 2*i - 1
 * one on the left, with opposite orientation, i -> 2*i
 *
 * convention: arcs are numbered from 0 to 2n-1 (n is the number of nodes)
 * after a checkerboard coloring of the regions, arcs are divided into:
 * even: white on right
 * odd: white on left
 *
 * nodes are numbered from 0 to n-1.  The function node = node(arc) indicates the
 * arrival node of the given arc.  After a checkerboard coloring, exactly two arcs, one
 * even and one odd are the preimages of the "node" function, we shall call them
 * arc = even(node) and arc = odd(node).
 * If originated by a dtcode, the three functions satisfy:
 *
 * even(i) = 2*i
 * odd(i) = 2*i+1
 * node(i) = i/2 if i is even, = (i-1)/2 if i is odd
 *
 * since this could be used in a different context/numbering we provide explicit functions
 * (to be changed in the future) for these mappings
 */

int dt_even (int i)
{
  if (dteven == 0) return (2*i);
  return (dteven[i]);
}

int dt_odd (int i)
{
  if (dtodd == 0) return (2*i + 1);
  return (dtodd[i]);
}

int dt_node (int i)
{
  assert (dt_iseven (i));
  if (dtnode == 0) return (i/2);
  return (dtnode[i]);
}

int dt_iseven (int i)
{
  int isevennumber;

  isevennumber = ((i % 2) == 0);
  if (ori_component == 0) return (isevennumber);
  if (ori_component[componentoflabel[i]] > 0) return (isevennumber);
  return (1 - isevennumber);
}

/*
 * nextlabel and prevlabel should cicle through the arcs of each component
 * of a link
 */

int
nextlabel (int label)
{
  int component;
  int startnextcomponent = numlabels;
  int startcomponent = 0;

  if (link_components > 1)
  {
    component = componentoflabel[label];
    startcomponent = start_component[component];
    if (component < link_components - 1) startnextcomponent = start_component[component+1];
  }
  label++;
  if (label >= startnextcomponent) label -= startnextcomponent - startcomponent;
  return (label);
}

int
prevlabel (int label)
{
  int component;
  int startnextcomponent = numlabels;
  int startcomponent = 0;

  if (link_components > 1)
  {
    component = componentoflabel[label];
    startcomponent = start_component[component];
    if (component < link_components - 1) startnextcomponent = start_component[component+1];
  }
  label--;
  if (label < startcomponent) label += startnextcomponent - startcomponent;
  return (label);
}

/*
 * display_arcs_from_arcs
 * pair of arcs parallel to arcs of the diagram
 * the first numbered (+1) is the one on the right
 * the second numbered (+2) is the one on the left
 */

void
display_arcs_from_arcs (struct sketch *s)
{
  struct arc *arc1, *arc2;
  int i;

  for (i = 0; i < numlabels; i++)
  {
    arc1 = newarc (s);
    arc1->tag = 2*i+1;
    if (arc1->next && arc1->tag > arc1->next->tag)
    {
      s->arcs = arc1->next;
      arc1->next = 0;
      insert_arc_in_list (arc1, s->arcs);
    }
    arc2 = newarc (s);
    arc2->tag = 2*i+2;
    if (arc2->next && arc2->tag > arc2->next->tag)
    {
      s->arcs = arc2->next;
      arc2->next = 0;
      insert_arc_in_list (arc2, s->arcs);
    }
    arc1->depths = (int *) malloc (sizeof (int));
    arc2->depths = (int *) malloc (sizeof (int));
    arc1->depths[0] = arc2->depths[0] = 0;
    arc1->depthsdim = arc2->depthsdim = 1;
    arc1->cusps = arc2->cusps = 0;
    arc1->endpoints = arc2->endpoints = 2;
    //printf ("Arc %d: [0];\n", 2*i - 1);
    //printf ("Arc %d: [0];\n", 2*i);
  }
}

/*
 * display_arcs_from_nodes
 */

void
display_arcs_from_nodes (struct sketch *s)
{
  struct arc *arc1, *arc2, *arc3, *arc4;
  int i, overpass;
  int offset = 2*numlabels;

  //printf ("# start of node arcs\n");
  for (i = 0; i < numnodes; i++)
  {
    arc1 = newarc (s);
    arc1->tag = offset + 4*i + 1;
    if (arc1->next && arc1->tag > arc1->next->tag)
    {
      s->arcs = arc1->next;
      arc1->next = 0;
      insert_arc_in_list (arc1, s->arcs);
    }
    arc2 = newarc (s);
    arc2->tag = offset + 4*i + 2;
    if (arc2->next && arc2->tag > arc2->next->tag)
    {
      s->arcs = arc2->next;
      arc2->next = 0;
      insert_arc_in_list (arc2, s->arcs);
    }
    arc3 = newarc (s);
    arc3->tag = offset + 4*i + 3;
    if (arc3->next && arc3->tag > arc3->next->tag)
    {
      s->arcs = arc3->next;
      arc3->next = 0;
      insert_arc_in_list (arc3, s->arcs);
    }
    arc4 = newarc (s);
    arc4->tag = offset + 4*i + 4;
    if (arc4->next && arc4->tag > arc4->next->tag)
    {
      s->arcs = arc4->next;
      arc4->next = 0;
      insert_arc_in_list (arc4, s->arcs);
    }
    arc1->depths = (int *) malloc (sizeof (int));
    arc2->depths = (int *) malloc (sizeof (int));
    arc3->depths = (int *) malloc (sizeof (int));
    arc4->depths = (int *) malloc (sizeof (int));
    overpass = 0;
    if (dtsign[dt_even(i)] < 0) overpass = 2;
    arc1->depths[0] = arc3->depths[0] = 2 - overpass;
    arc2->depths[0] = arc4->depths[0] = overpass;
    arc1->depthsdim = arc2->depthsdim = arc3->depthsdim = arc4->depthsdim = 1;
    arc1->cusps = arc2->cusps = arc3->cusps = arc4->cusps = 0;
    arc1->endpoints = arc2->endpoints = arc3->endpoints = arc4->endpoints = 2;
    //printf ("Arc %d: [%d];\n", offset + 4*i - 3, 2 - overpass);
    //printf ("Arc %d: [%d];\n", offset + 4*i - 2, overpass);
    //printf ("Arc %d: [%d];\n", offset + 4*i - 1, 2 - overpass);
    //printf ("Arc %d: [%d];\n", offset + 4*i - 0, overpass);
  }
}

/*
 * display_regions
 * the first red region starts left of arc labelled 0
 */

void
display_regions (struct sketch *s)
{
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *a;
  int *redtagged, *bluetagged;
  int redregionnum = 0;
  int blueregionnum = 0;
  int i, arc, atag;

  redtagged = (int *) malloc ( numlabels * sizeof(int) );
  bluetagged = (int *) malloc ( numlabels * sizeof(int) );
  for (i = 0; i < numlabels; i++) redtagged[i] = bluetagged[i] = 0;

  //printf ("Displaying red regions...\n");

  for (i = 0; i < numlabels; i++)
  {
    if (redtagged[i]) continue;
    redregionnum++;  /* found new region */
    blast = 0;
    //printf ("New red region: %d:\n", redregionnum);
    region = newregion (s);
    region->tag = redregionnum - 1;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }
    bl = newborderlist_tail (region);
    if (redregionnum == 1)
    {
      bl->isexternal = 1;
      bl = newborderlist_tail (region);
    }
    //printf ("Region %d: ", redregionnum - 1);
    //if (redregionnum == 1) printf ("() ");
    //printf ("(");
    for (arc = i;;)
    {
      redtagged[arc] = redregionnum;
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

      if ( dt_iseven (arc) )
      { /* "odd" arc */
        if (dt_realization[arc] > 0)
        {
          arc = nextlabel(dt_involution[arc]);
          atag = 2*arc + 2;
          //printf ("-a%d ", 2*arc + 2);
        } else {
          arc = dt_involution[arc];
          atag = 2*arc + 1;
          //printf ("-a%d ", 2*arc + 1);
        }
      } else {  /* "even" arc */
        if (dt_realization[prevlabel(arc)] > 0)
        {
          arc = dt_involution[prevlabel(arc)];
          atag = 2*arc + 1;
          //printf ("-a%d ", 2*arc + 1);
        } else {
          arc = nextlabel(dt_involution[prevlabel(arc)]);
          atag = 2*arc + 2;
          //printf ("-a%d ", 2*arc + 2);
        }
      }
      for (a = s->arcs; a; a = a->next)
      {
        if (a->tag == atag)
        {
          b->info = a;
          break;
        }
      }
      if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

      if (arc == i) break;
      //{
      //  printf (");\n");
      //  break;
      //}
    }
  }

  //printf ("Displaying blue regions (counterclockwise)...\n");

  for (i = 0; i < numlabels; i++)
  {
    if (bluetagged[i]) continue;
    blueregionnum++;  /* found new region */
    blast = 0;
    //printf ("New blue region: %d:\n", blueregionnum);
    region = newregion (s);
    region->tag = redregionnum + blueregionnum - 1;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }
    bl = newborderlist_tail (region);
    //printf ("Region %d: ", redregionnum + blueregionnum - 1);
    //printf ("(");
    for (arc = i;;)
    {
      bluetagged[arc] = blueregionnum;
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
      if ( dt_iseven (arc) )
      { /* "odd" arc */
        if (dt_realization[prevlabel(arc)] > 0)
        {
          arc = dt_involution[prevlabel(arc)];
          atag = 2*arc + 1;
          //printf ("-a%d ", 2*arc + 1);
        } else {
          arc = nextlabel(dt_involution[prevlabel(arc)]);
          atag = 2*arc + 2;
          //printf ("-a%d ", 2*arc + 2);
        }
      } else {  /* "even" arc */
        if (dt_realization[arc] > 0)
        {
          arc = nextlabel(dt_involution[arc]);
          atag = 2*arc + 2;
          //printf ("-a%d ", 2*arc + 2);
        } else {
          arc = dt_involution[arc];
          atag = 2*arc + 1;
          //printf ("-a%d ", 2*arc + 1);
        }
      }
      for (a = s->arcs; a; a = a->next)
      {
        if (a->tag == atag)
        {
          b->info = a;
          break;
        }
      }
      if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

      if (arc == i) break;
      //{
      //  printf (");\n");
      //  break;
      //}
    }
  }

  assert (redregionnum + blueregionnum == numregions);
}

/*
 * display_regions_from_arcs
 * nodes are indexed by using the "odd" label i: i/2
 * then we have 4 arcs for each node, so there is a multiplying factor 4
 * for a final 2*(i-1) + h
 * h = 0,1,2,3
 */

void
display_regions_from_arcs (struct sketch *s)
{
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *a;
  int i, atag;
  int offset = numregions;
  int arcsoffset = 2*numlabels + 1;
  int arcahead, arcbehind;

  //printf ("# elongated regions corresponding to arcs\n");
  for (i = 0; i < numlabels; i++)
  {
    blast = 0;
    region = newregion (s);
    region->tag = offset + i;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }

    if ( dt_iseven (i) )
    {
      arcahead = 4*dt_node(i);
      arcbehind = 4*dt_node((dt_involution[prevlabel(i)])) + 3;
    } else {
      arcahead = 4*dt_node(dt_involution[i]) + 1;
      arcbehind = 4*dt_node(prevlabel(i)) + 2;
    }
    //printf ("Region %d: ", offset + i + 1);
    bl = newborderlist_tail (region);
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
    atag = 2*i + 1;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
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
    atag = arcsoffset + arcahead;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
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
    atag = 2*i + 2;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
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
    atag = arcsoffset + arcbehind;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
    //printf ("(+a%d -a%d +a%d -a%d);\n", 2*i - 1, arcsoffset + arcahead, 2*i, arcsoffset + arcbehind);
  }
}

/*
 * display_regions_from_nodes
 */

void
display_regions_from_nodes (struct sketch *s)
{
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *a;
  int i, atag, oddarc, ori;
  int offset = numregions + numlabels;
  int arcsoffset = 2*numlabels + 1;

  //printf ("# small nodal regions\n");
  for (i = 0; i < numnodes; i++)
  {
    blast = 0;
    region = newregion (s);
    region->tag = offset + i;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }
    bl = newborderlist_tail (region);
    oddarc = dt_even(i);
    ori = dt_realization[oddarc];
    oddarc = arcsoffset + 4*i;

    atag = oddarc;
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
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    atag = oddarc + 2 - ori;
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
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    atag = oddarc + 2;
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
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    atag = oddarc + 2 + ori;
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
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    //printf ("Region %d: (+a%d +a%d +a%d +a%d);\n", offset + i, oddarc,
    //                                                           oddarc + 2 - ori,
    //                                                           oddarc + 2,
    //                                                           oddarc + 2 + ori);
  }
}

/*
 * reconstruct dt_realization of nodes
 * if nonzero, gregionsign points to a vector of given handedness:
 * 0: value unassigned
 * +1: if arriving at the odd-numbered label, the crossed strand
 *     is oriented from right to left
 * -1: if arriving at the odd-numbered label, the crossed strand
 *     is oriented from left to right
 *
 * This doesn't work well (see Haken satellite example), perhaps we should try
 * to incorporate the routine "realize" by Jim Hoste in decode_new_DT.c (knotscape)
 */

int
inherit (int startlabel)
{
  int isinside;
  int endlabel, label;
  int expansions = 0;

  // printf ("searching inheritance for cycle starting at %d\n", startlabel + 1);

  /* there is a cicle starting ad startlabel + 1 and ending
   * at dt_involution[startlabel] + 1
   */

  endlabel = dt_involution[startlabel];
  assert (nextlabel(startlabel) != dt_involution[startlabel]);  /* no tight loops allowed */
  assert (nextlabel(endlabel) != dt_involution[endlabel]);  /* no tight loops allowed */

  /* now find where the path continuing from dt_involution[startlabel] + 1
   * first intersects the loop
   */

  isinside = 0;
  for (label = nextlabel(endlabel); label != startlabel; label = nextlabel(label))
  {
    if (isinpath (dt_involution[label], startlabel, endlabel))
    {
      if (isinside)
      {
        // printf ("Sign Inheritance between %d and %d (EXITING)\n", startlabel + 1, label + 1);
        expansions += propagate (-dt_realization[startlabel], label);
        expansions += propagate (-dt_realization[label], startlabel);
      } else {
        // printf ("Sign Inheritance between %d and %d (ENTERING)\n", startlabel + 1, label + 1);
        expansions += propagate (dt_realization[startlabel], label);
        expansions += propagate (dt_realization[label], startlabel);
      }
      isinside = 1-isinside;
    }
  }
  return (expansions);
}

/*
 *
 */

int
isconsistent (void)
{
  int i;
  int countregions = 0;

  for (i = 0; i < numlabels; i++) tagged[i] = 0;

  /* each arc should be walked once positively and once negatively
   * tag 0
   * tag 1 +
   * tag 2 -
   * tag 3 +-
   */

  while (1)
  {
    for (i = 0; i < numlabels; i++)
    {
      if (tagged[i] < 3) break;
    }
    if (i >= numlabels) break; /* complete tour! everything fine */

    countregions++;
    /* i is an arc that has not been walked twice */
    if (tagged[i] == 0 || tagged[i] == 2)
    {
      if (!tour_of_region (i, 1)) return (0);  /* something went wrong during tour */
    } else {
      if (!tour_of_region (i, -1)) return (0);  /* something went wrong during tour */
    }
  }
  if (countregions != numregions) return (0);
  //printf ("correct number of regions!\n");
  /* we NEED to find a test that works for sure! This is not enough */
  return (1);
}

/*
 * return 0 if something goes wrong
 */

int
tour_of_region (int startlabel, int startvelocity)
{
  int label = startlabel;
  int velocity = startvelocity;
  int i, node_to_tag;

  for (i = 0; i < numlabels; i++) marknodes[i] = 0;  /* will mark traversed labels */

  while (1)
  {
    node_to_tag = label;
    if (velocity < 0) node_to_tag = prevlabel(label);
    if (marknodes[node_to_tag] != 0) return (0);
    marknodes[node_to_tag] = 1;
    if (velocity > 0)
    {
      if ((tagged[label] & 1)) return (0); /* should not be already visited */
      tagged[label] += 1;
    } else {
      if ((tagged[label] & 2)) return (0);
      tagged[label] += 2;
    }
    walk_left (&label, &velocity);
    if (label == startlabel && velocity == startvelocity) return (1);
  }
  return (0);
}

/*
 * walk around a region conterclockwise
 */

void
walk_left (int *labelpt, int *velocitypt)
{
  if (*velocitypt > 0)
  {
    if (dt_realization[*labelpt] > 0)
    {
      /* velocity remains positive */
      *labelpt = nextlabel(dt_involution[*labelpt]);
    } else {
      *velocitypt = - *velocitypt;
      *labelpt = dt_involution[*labelpt];
    }
  } else {
    if (dt_realization[prevlabel(*labelpt)] > 0)
    {
      /* velocity remains negative */
      *labelpt = dt_involution[prevlabel(*labelpt)];
    } else {
      *velocitypt = - *velocitypt;
      *labelpt = nextlabel(dt_involution[prevlabel(*labelpt)]);
    }
  }
}

int
isinpath (int i, int start, int end)
{
  if (end >= start) return (i >= start && i <= end);

  if (i >= start) return (1);
  if (i <= end) return (1);
  return (0);
}

int
propagate (int sign, int label)
{
  if (sign == 0) return (0);
  if (dt_realization[label] != 0) return (0);

  dt_realization[label] = sign;
  dt_realization[dt_involution[label]] = -sign;
  return (1);
}

/*
 * this is taken verbatim from decode_new_DT.c of knotscape
 * based on Jim Hoste previous code
 */

static void dt_realize(
    int     *anInvolution,
    int     *aRealization,
    int     aNumCrossings)
{
    /*
     *  I'm hoping to read through -- and understand -- Dowker
     *  and Thistlethwaite's paper.  But for the moment I'll
     *  try to splice in some of Jim's code and hope I'm interpreting
     *  it correctly.
     */

    int N;
    int i,j;
    int *modTWO_N;
    int *modN;
    int *seq;
    int *emb;
    int *A,*D;
    int Aempty,Dempty;
    int OkSoFar;
    int *phi;
    int x;

    N = aNumCrossings;

    /*
     *  Allocate local arrays.
     */
    modTWO_N    = (int *) malloc(4 * N * sizeof(int));
    modN        = (int *) malloc(2 * N * sizeof(int));
    seq         = (int *) malloc(4 * N * sizeof(int));
    emb         = (int *) malloc(2 * N * sizeof(int));
    A           = (int *) malloc(2 * N * sizeof(int));
    D           = (int *) malloc(2 * N * sizeof(int));
    phi         = (int *) malloc(2 * N * sizeof(int));

    /*create the modTWO_N array*/
    for(i=0;i<2*N;i++){
        modTWO_N[i]=i;
        modTWO_N[i+2*N]=i;
    }
    /*create the modN array*/
    for(i=0;i<N;i++){
        modN[i]=i;
        modN[i+N]=i;
    }
    /* get seq and height from DT code*/
    /* seq is two copies of full DT involution on crossings numbered 0 to 2N-1 */

    /*
     *  Concatenate two copies of theInvolution[] to obtain seq[].
     */
    for(i = 0; i < 2*N; i++)
    {
        seq[i      ] = anInvolution[i];
        seq[i + 2*N] = anInvolution[i];
    }

    /* begin realizability routine to recover embedding of projection */
    /*zero emb, A, and D. A and D will only contain zeroes and ones*/
    for(i=0;i<2*N;i++){
        emb[i]=A[i]=D[i]=0;
    }
    /*set initial conditions*/
    OkSoFar=A[0]=A[seq[0]]=1;
    emb[0]=1;
    emb[seq[0]]=-1;

    /* see if A is empty, ie is all zeroes*/
    for(j=0;j<2*N-1 && !A[j];j++)
        /*nothing*/;
    Aempty=!A[j];

    while(!Aempty && OkSoFar){
        /* let i be least member of A*/
        for(i=0; !A[i]; i++)
            /*nothing*/;
        /*determine phi for this value of i*/
        phi[i]=1;
        for(j=i+1;j<i+2*N;j++){
            phi[modTWO_N[j]]=(seq[j]>=i && seq[j]<=seq[i]) ? -phi[modTWO_N[j-1]]:phi[modTWO_N[j-1]];
        }
        /*establish D*/
        for(j=0; j<i; j++)
            D[j]=1;
        for(j=seq[i]+1; j<2*N; j++)
            D[j]=1;
        /* see if D is empty, ie is all zeroes*/
        for(j=0;j<2*N-1 && !D[j];j++)
            ;
        Dempty=!D[j];
        while(!Dempty && OkSoFar){
            /*let x be least member of D*/
            for(x=0; !D[x]; x++)
                /*nothing*/;
            if(x<i){
                if(seq[x]<i || seq[x]>seq[i]){
                    if(phi[x]*phi[seq[x]]==1){
                        D[x]=D[seq[x]]=0;
                    }
                    else{
                        OkSoFar=0;
                    }
                }
                else{
                    if(emb[x] != 0){/* emb[x] is already defined*/
                        if(phi[x]*phi[seq[x]]*emb[i]==emb[x])
                            D[x]=0;
                        else
                            OkSoFar=0;
                    }
                    else{/* emb[x] is not yet defined. Hence x<>0 and x<>seq[0]*/
                        emb[x]=phi[x]*phi[seq[x]]*emb[i];
                        emb[seq[x]]=-emb[x];
                        D[x]=0;
                        if( modTWO_N[abs(seq[x]-seq[x-1])]==1){
                            /*nothing*/
                        }
                        else{
                            A[x]=A[seq[x]]=1;
                        }
                    }
                }
            }
            else{/*x>seq[i]*/
                if(seq[x]<i || seq[x]>seq[i]){
                    D[x]=D[seq[x]]=0;
                }
                else{
                    if(emb[x]!=0){
                        D[x]=0;
                    }
                    else{
                        emb[x]=phi[x]*phi[seq[x]]*emb[i];
                        emb[seq[x]]=-emb[x];
                        D[x]=0;
                        if( modTWO_N[abs(seq[x]-seq[x-1])]==1){
                            /*nothing*/
                        }
                        else{
                            A[x]=A[seq[x]]=1;
                        }
                    }
                }
            }
            /* see if D is empty, ie is all zeroes*/
            for(j=0;j<2*N-1 && !D[j];j++)
                ;
            Dempty=!D[j];
        }/*end of while*/
        A[i]=0;
        A[seq[i]]=0;
        /* see if A is empty, ie is all zeroes*/
        for(j=0;j<2*N-1 && !A[j];j++)
            /*nothing*/;
        Aempty=!A[j];
    }/*end of while*/

    /*end of realizability routine*/

    if(!OkSoFar){/* sequence is not realizable*/
        uFatalError("realize", "decode_DT");
    }

    /*
     *  Convert emb[] to aRealization[].
     */
    for(i = 0; i < 2*N; i++)
        //aRealization[i] = (emb[i] == +1);
        aRealization[i] = emb[i];

    /*
     *  Free local arrays.
     */
    free(modTWO_N);
    free(modN);
    free(seq);
    free(emb);
    free(A);
    free(D);
    free(phi);
}

/*
 * print loiv data using dtcode or gausscode syntax
 */

void
printloiv (struct vecofintlist *loiv)
{
  int i;
  struct vecofintlist *lv;

  switch (loiv->type)
  {
    case LOIV_ISDTCODE:
    case LOIV_ISRDTCODE:
      printf ("dtcode {");
      for (lv = loiv; lv; lv = lv->next)
      {
        printf ("[");
        for (i = 0; i < lv->len; i++)
        {
          if (i > 0) printf (" ");
          printf ("%d", lv->vec[i]);
          if (lv->handedness && lv->handedness[i]) printf ("%c", (lv->handedness[i] > 0)?'>':'<');
        }
        printf ("]");
      }
      printf ("}\n");
      break;

    case LOIV_ISGAUSSCODE:
      printf ("gausscode {");
      for (lv = loiv; lv; lv = lv->next)
      {
        printf ("{");
        for (i = 0; i < lv->len; i++)
        {
          if (i > 0) printf (" ");
          printf ("%d", lv->vec[i]);
          if (lv->handedness && lv->handedness[i]) printf ("%c", (lv->handedness[i] > 0)?'>':'<');
        }
        printf ("}");
      }
      printf ("}\n");
      break;

    default:
      printf ("unknown type for loiv: %d\n", loiv->type);
      break;
  }
}
