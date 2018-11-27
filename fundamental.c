/*
 * computation of the fundamental group of the interior of
 * the surface
 */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "contour.h"
#include "laurent.h"
#include "fundamental.h"
#include "fox.h"
#include "representations.h"
#include "alexander.h"
#include "parser.h"

#define MAXGENERATORS 26

extern struct global_data globals;
extern int debug;
extern int quiet;
extern int verbose;
extern int interactive;

void
compute_fundamental (struct ccomplex *cc, int action)
{
  int res, ccnum, errors = 0;
  struct ccomplexcc *cccc;

  complex_melt (cc);
  if (debug) cellcomplex_checkconsistency (cc);

  if (debug) printf ("Constructing spanning tree\n");
  ccnum = find_spanning_tree (cc);
  if (debug) printf ("Found %d connected components\n", ccnum);

  if (ccnum == 0 && !quiet) printf ("Empty set, there is no fundamental group.\n");
  assert (ccnum >= 0);

  for (cccc = cc->cc; cccc; cccc = cccc->next)
  {
    if (ccnum > 1 && !quiet) printf ("\nConnected component %d:\n", cccc->tag);
    cccc->p = compute_fundamental_single (cc, cccc);
    switch (action)
    {
      case ACTION_FUNDAMENTAL:
        fundamental_group (cccc->p);
        break;
      case ACTION_AFUNDAMENTAL:
        abelianized_fundamental_group (cccc->p);
        break;
      case ACTION_FOXJACOBIAN:
        if (globals.simplifypresentation) simplify_presentation (cccc->p);
        foxjacobian (cccc->p);
        break;
      case ACTION_ALEXANDER:
        if (globals.simplifypresentation) simplify_presentation (cccc->p);
        res = alexander (cccc->p);
        if (res == 0) errors++;
        break;
      case ACTION_LINKINGNUMBER:
        if (globals.simplifypresentation) simplify_presentation (cccc->p);
        linkingnumber (cccc->p);
        break;
      case ACTION_CCCOUNTSL2ZP:
        if (globals.simplifypresentation) simplify_presentation (cccc->p);
        cccountsl2zp (cccc->p);
        break;
      case ACTION_CCCOUNTPSL2Q:
        if (globals.simplifypresentation) simplify_presentation (cccc->p);
        cccountpsl2q (cccc->p);
        break;
      case ACTION_CCCOUNTSN:
        if (globals.simplifypresentation) simplify_presentation (cccc->p);
        cccountsn (cccc->p);
        break;
      default:
        printf ("Invalid action: %d\n", action);
        break;
    }
  }
  if (errors) exit (1);
}

void
fundamental_group (struct presentation *p)
{
  if (verbose >= 2) print_presentation (p);
  if (interactive >= 2) fg_interactive (p);
  if (globals.simplifypresentation) simplify_presentation (p);
  if (globals.preabelian) topreabelian (p);
  if (interactive) fg_interactive (p);
  print_presentation (p);
  if (verbose) print_exponent_matrix (p);
}

void
abelianized_fundamental_group (struct presentation *p)
{
  if (debug) print_presentation (p);
  if (interactive >= 2) fg_interactive (p);
  simplify_presentation (p);
  if (globals.preabelian) topreabelian (p);
  if (interactive) fg_interactive (p);
  print_invariant_factors (p);
  if (verbose) print_exponent_matrix (p);
}

/*
 * compute presentation of fundamental group of one component
 */

struct presentation *
compute_fundamental_single (struct ccomplex *cc, struct ccomplexcc *cccc)
{
  struct presentation *p;
  struct presentationrule *rule;
  struct ccomplexface *faces = cc->faces;
  struct ccomplexface *face;
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexarc *arc;
  struct ccomplexnode *nodes = cc->nodes;
  struct ccomplexnode *node;
  int gennum, length, narc, j, n, m, i, s;
  int rulesnum = 0;
  int *ivec, *variable;

  assert (cccc->betti_2 >= 0);
  p = (struct presentation *) malloc (sizeof (struct presentation));
  p->gennum = 0;
  p->rules = p->elements = 0;
  p->characteristic = 1;  /* the spanning tree has characteristic 1 */

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type == CC_REMOVED) continue;
    if (arc->isinspanningtree) continue;
    node = cc->nodes + arc->enda;
    if (node->cc != cccc) continue;
    p->characteristic--;  /* each further arc counts -1 */
    p->gennum++;
  }

  if (p->gennum > 0)
  {
    variable = (int *) malloc (p->gennum*sizeof (int));
    for (gennum = 0, n = 0; n < cc->arcdim; n++)
    {
      arc = cc->arcs + n;
      if (arc->type == CC_REMOVED) continue;
      if (arc->isinspanningtree) continue;
      node = cc->nodes + arc->enda;
      if (node->cc != cccc) continue;
      variable[gennum++] = n;
    }

    /* loop on the list of rules */

    for (m = 0; m < cc->facedim; m++)
    {
      face = faces + m;
      if (face->type == CC_REMOVED) continue;
      ivec = face->faceborder;
      arc = arcs + onarc2narc (ivec[0]);
      node = nodes + arc->enda;
      if (node->cc != cccc) continue;
      /* we have a new rule */
      p->characteristic++;  /* each face counts +1 */
      rulesnum++;
      for (length = 0, i = 0; i < face->facebordernum; i++)
      {
        narc = onarc2narc (ivec[i]);
        arc = arcs + narc;
        if (!arc->isinspanningtree) length++;
      }
      rule = (struct presentationrule *) malloc (length*sizeof (int) + sizeof (struct presentationrule));
      rule->length = length;
      rule->next = 0;
      if (p->rules != 0) rule->next = p->rules;
      p->rules = rule;

      j = 0;
      for (i = 0; i < face->facebordernum; i++)
      {
        narc = onarc2narc (ivec[i]);
        arc = arcs + narc;
        if (arc->isinspanningtree) continue;
        for (s = 0; s < p->gennum; s++)
        {
          if (variable[s] == narc) break;
        }
        assert (s < p->gennum);
        s++;
        if (ivec[i] < 0) s *= -1;
        rule->var[j++] = s;
      }
    }
    free (variable);
  }

  p->espected_deficiency = 1 - p->characteristic + cccc->spherical_voids;
  assert (cccc->betti_2 >= 0);
  if (verbose)
  {
    /* b1 = b0 + b2 - C */
    printf ("Euler characteristic: %d\n", p->characteristic);
    printf ("Betti numbers: b0 = 1, b1 = %d, b2 = %d\n", 1 + cccc->betti_2 - p->characteristic, cccc->betti_2);
    printf ("Number of spherical voids in S3: %d\n", cccc->spherical_voids);
    printf ("Initial deficiency: %d\n", gennum - rulesnum);
    printf ("Espected deficiency: %d\n", p->espected_deficiency);
  }
  return (p);
}

/*
 * simplify presentation
 */

/*
 * local profiles
 */

int simplify_presentation2 (struct presentation *p);
int sp_eliminatevar (struct presentation *p);
void sp_do_eliminatevar (struct presentation *p, int);
int sp_removeemptyrules (struct presentation *p);
int sp_simplifyword (struct presentation *p);
int sp_findcommonsubstring (struct presentation *p, int bidirectional);
int sp_rotate_to_canonical (struct presentation *p);
int sp_rotate_to_canonical_single (struct presentationrule *r, int tryrevert);
int sp_findduplicatewords (struct presentation *p);
int sp_tryeliminatevariable (struct presentation *p);
struct presentationrule *sp_do_substitute (struct presentationrule *r, int subvar, int *replaceword, int optsublen);
struct presentationrule *sp_remove_marked_elements (struct presentationrule *r);
void tietze_mulright (struct presentation *p,
                       struct presentationrule *r1,
                       struct presentationrule *r2, int expon);
int preabelian_step (struct presentation *p, int row, struct presentationrule *rcol);

int
simplify_presentation (struct presentation *p)
{
  int goon = 1;
  int count = 0;

  while (goon)
  {
    goon = 0;
    goon += simplify_presentation2 (p);
    goon += sp_tryeliminatevariable (p);
    count += goon;
  }
  return (count);
}

int
simplify_presentation2 (struct presentation *p)
{
  int goon = 1;
  int count = 0;

  while (goon)
  {
    if (debug) printf ("Simplifying presentation: %d generators\n", p->gennum);
    goon = 0;
    goon += sp_eliminatevar (p);
    goon += sp_simplifyword (p);
    goon += sp_findcommonsubstring (p, 1);
    //goon += sp_rotate_to_canonical (p);
    //goon += sp_findduplicatewords (p);
    goon += sp_removeemptyrules (p);
    count += goon;
  }
  count += sp_rotate_to_canonical (p);
  return (count);
}

int
sp_tryeliminatevariable (struct presentation *p)
{
  struct presentationrule *r;
  struct presentationrule *optsubrule;
  int goon = 1;
  int count = 0;
  int optsublen = INT_MAX;
  int i, j, v, lpos, lcount;
  int optsubpos, optsubvar, subvar;
  int *replaceword;

  while (goon)
  {
    if (debug) printf ("in sp_tryeliminatevariable: %d generators\n", p->gennum);
    goon = 0;
    optsublen = INT_MAX;
    optsubvar = 0;
    for (v = 1; v <= p->gennum; v++)
    {
      for (r = p->rules; r; r = r->next)
      {
        if (r->length < 2) continue;
        if (r->length >= optsublen) continue;
        for (lcount = 0, i = 0; i < r->length; i++)
        {
          if (abs(r->var[i]) == v) {lcount++; lpos = i;}
        }
        if (lcount != 1) continue; // variable 'v' cannot be eliminated
        optsublen = r->length;
        optsubvar = v;
        optsubrule = r;
        optsubpos = lpos;
      }
    }
    if (optsubvar > 0)
    {
      goon++;
      count++;
      optsublen--;
      replaceword = (int *) malloc (optsublen * sizeof (int));
      subvar = - optsubrule->var[optsubpos];
      assert (optsubvar == abs(subvar));
      for (i = 0; i < optsublen; i++)
      {
        j = i + optsubpos + 1;
        if (j >= optsubrule->length) j -= optsubrule->length;
        replaceword[i] = optsubrule->var[j];
        assert (abs (replaceword[i]) != optsubvar);
      }
      for (r = p->rules; r; r = r->next) r = sp_do_substitute (r, subvar, replaceword, optsublen);
      p->rules = sp_remove_marked_elements (p->rules);
      for (r = p->elements; r; r = r->next) r = sp_do_substitute (r, subvar, replaceword, optsublen);
      p->elements = sp_remove_marked_elements (p->elements);
      sp_do_eliminatevar (p, subvar);
      free (replaceword);
      count += simplify_presentation2 (p);
    }
  }

  return (count);
}

/*
 * This actually *adds* the modified rule as a new rule pointed by the old one.
 * Rhe old rule gets marked by setting the length to an (invalid) negative value.
 */

struct presentationrule *
sp_do_substitute (struct presentationrule *r, int subvar, int *replaceword, int wordlen)
{
  struct presentationrule *newr;
  int i, j, k, newlength, lcount;

  assert (subvar != 0);
  assert (wordlen >= 1);
  if (wordlen == 1) // in-place substitution
  {
    for (i = 0; i < r->length; i++)
    {
      if (r->var[i] == subvar) r->var[i] = replaceword[0];
      if (r->var[i] == - subvar) r->var[i] = - replaceword[0];
    }
    return (r);
  } else {
    // count number of occurrences of subvar
    for (lcount = 0, i = 0; i < r->length; i++)
      if (r->var[i] == subvar || r->var[i] == - subvar) lcount++;
    newlength = r->length - lcount + lcount*wordlen;
    newr = (struct presentationrule *) malloc (newlength*sizeof (int) + sizeof (struct presentationrule));
    newr->length = newlength;
    newr->next = r->next;
    r->next = newr;
    for (i = 0, j = 0; i < r->length; i++)
    {
      if (r->var[i] != subvar && r->var[i] != - subvar) newr->var[j++] = r->var[i];
      if (r->var[i] == subvar)
      {
        for (k = 0; k < wordlen; k++) newr->var[j++] = replaceword[k];
      }
      if (r->var[i] == - subvar)
      {
        for (k = wordlen - 1; k >= 0; k--) newr->var[j++] = -replaceword[k];
      }
    }
    assert (j == newlength);
    // last operation: mark old rule
    r->length = -1;
    return (newr);
  }
}

struct presentationrule *
sp_remove_marked_elements (struct presentationrule *r)
{
  struct presentationrule *rnext;

  if (r == 0) return (0);
  if (r->length < 0)
  {
    rnext = r->next;
    free (r);
    return (sp_remove_marked_elements (rnext));
  }
  r->next = sp_remove_marked_elements (r->next);
  return (r);
}

int
sp_simplifyword (struct presentation *p)
{
  struct presentationrule *r;
  int goon = 1;
  int count = 0;
  int i, j, inext, ifound1, ifound2;

  while (goon)
  {
    goon = 0;
    for (r = p->rules; r; r = r->next)
    {
      if (r->length < 2) continue;
      ifound1 = -1;
      for (i = 0; i < r->length; i++)
      {
        inext = i+1;
        if (inext >= r->length) inext = 0;
        if (r->var[i] + r->var[inext] == 0)
        {
          /* trovato due caratteri che si semplificano */
          ifound1 = i;
          ifound2 = inext;
          break;
        }
      }
      if (ifound1 >= 0)
      {
        goon++;
        count++;
        for (i = 0, j = 0; i < r->length; i++)
        {
          if (i != ifound1 && i != ifound2) r->var[j++] = r->var[i];
        }
        r->length -= 2;
        assert (j == r->length);
      }
    }
    for (r = p->elements; r; r = r->next)
    {
      if (r->length < 2) continue;
      ifound1 = -1;
      for (i = 0; i < r->length - 1; i++)
      {
        if (r->var[i] + r->var[i+1] == 0)
        {
          /* trovato due caratteri che si semplificano */
          for (j = i; j + 2 < r->length; j++) r->var[j] = r->var[j+2];
          goon++;
          count++;
          r->length -= 2;
          break;
        }
      }
    }
  }
  return (count);
}

int
sp_rotate_to_canonical (struct presentation *p)
{
  struct presentationrule *r;
  int count = 0;

  for (r = p->rules; r; r = r->next)
  {
    count += sp_rotate_to_canonical_single (r, 1);
  }
  return (count);
}

int
sp_rotate_to_canonical_single (struct presentationrule *r, int tryrevert)
{
  int count = 0;
  int i, j, iopt, p1, p2;
  int firstvar, jsaved;
  int betterinverse;

  if (r->length <= 0) return (0);  /* already canonical */
  if (r->length == 1)
  {
    if (r->var[0] > 0) return (0); /* already canonical */
    r->var[0] = -r->var[0];        /* invert single variable */
    return (1);
  }

  iopt = 0;
  for (i = 1; i < r->length; i++)
  {
    /* lexicographical comparison */
    for (j = 0; j < r->length; j++)
    {
      p1 = (iopt + j)%r->length;
      p2 = (i + j)%r->length;
      if (r->var[p2]*r->var[p1] < 0)  /* different sign */
      {
        if (r->var[p2] > 0) iopt = i; /* positive wins */
        break;
      }
      if (r->var[p2] != r->var[p1])
      {
        if (abs(r->var[p2]) < abs(r->var[p1])) iopt = i;
        break;
      }
    }
  }

  if (iopt > 0)
  {
    count++;
    for (i = 0; i < iopt; i++)
    {
      firstvar = r->var[0];
      for (j = 1; j < r->length; j++)
        r->var[j-1] = r->var[j];
      r->var[r->length - 1] = firstvar;
    }
  }

  if (tryrevert == 0) return (count);

  /* is there a better inverse? */
  betterinverse = 0;
  for (i = 0; i < r->length; i++)
  {
    /* compare word with the inverse starting from i */
    /* remember to mentally change sign of var of reverted (inverse) */
    for (j = 0; j < r->length; j++)
    {
      // p1 = j;
      p2 = (i - j + r->length)%r->length;
      if (r->var[p2]*r->var[j] > 0)
      {
        if (r->var[p2] < 0) betterinverse = 1;
        break;
      }
      if (r->var[p2] + r->var[j] != 0)
      {
        if (abs(r->var[p2]) < abs(r->var[j])) betterinverse = 1;
        break;
      }
    }
  }

  if (betterinverse == 0) return (count);

  count++;
  /* invert! */
 
  for (i = 0, j = r->length - 1; i <= j; i++, j--)
  {
    /* exchange and invert i and j */
    if (i == j)
    {
      r->var[i] = -r->var[i];
      continue;
    }
    jsaved = r->var[j];
    r->var[j] = -r->var[i];
    r->var[i] = -jsaved;
  }
 
  return (count + sp_rotate_to_canonical_single (r, 0));
}

/*
 * A better simplification would be to find the largest
 * common sequence between two words, if it is long enough
 * we can decrease the length of the longest word.
 *
 * Simply removing duplicate words works for the simple
 * examples that we tested.
 */

int sp_fcs_r1r2 (struct presentationrule *r1,
                 struct presentationrule *r2, 
                 int bidirectional, int kminforrotation);

int
sp_findcommonsubstring (struct presentation *p, int bidirectional)
{
  int goon;
  int count = 0;
  struct presentationrule *r1, *r2, *lr, *sr;

  /* TODO: perhaps we should try a similar action between rules and selected elements
   */

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (r1 = p->rules; r1; r1 = r1->next)
    {
      for (r2 = r1->next; r2; r2 = r2->next)
      {
        lr = r1;
        sr = r2;
        if (r1->length < r2->length) {lr = r2; sr = r1;}
        goon = sp_fcs_r1r2 (lr, sr, bidirectional, sr->length/2 + 1);
        if (goon)
        {
          assert (2*abs(goon) > sr->length);
          count++;
          if (goon > 0)
          {
            tietze_mulright (p, lr, sr, -1);
          } else {
            tietze_mulright (p, lr, sr, 1);
          }
        }
        if (goon) break;
      }
      if (goon) break;
    }
    if (goon) sp_simplifyword (p);
    count += goon;
  }
  return (count);
}

/*
 * find the longest common subsequence in r1 and r2 (also in
 * reverse if bidirectional = 1)
 * if the length of the common subsequence is not less than
 * kminforrotation then also rotate the two strings such that
 * the common part is at the end.
 * if matching is reverse, then the matching part in r2
 * is at the beginning.
 * returns length of common subsequence (negative if reversed
 * in r2) if rotation takes place, otherwise return 0
 */

int
sp_fcs_r1r2 (struct presentationrule *r1,
             struct presentationrule *r2, 
             int bidirectional, int kminforrotation)
{
  //struct presentationrule *rsave;
  int i, j, k, ik, jk, iopt, jopt, direction;
  int saved, numofrot, kmax = 0;
  int minlen, retcode;

  minlen = r1->length;
  if (r2->length < minlen) minlen = r2->length;
  //if (r1->length < r2->length)
  //{
  //  rsave = r1;
  //  r1 = r2;
  //  r2 = rsave;
  //}
  for (i = 0; i < r1->length; i++)
  {
    for (j = 0; j < r2->length; j++)
    {
      for (k = 0; k <= minlen; k++)
      {
        ik = (i + k) % r1->length;
        jk = (j + k) % r2->length;
        if (k == minlen || r1->var[ik] != r2->var[jk])
        { /* sottostringa uguale di lunghezza k */
          if (k > kmax) /* migliore della precedente */
          {
            iopt = i;
            jopt = j;
            kmax = k;
            direction = 1;
          }
          break;
        }
      }
      if (bidirectional)
      {
        for (k = 0; k <= minlen; k++)
        {
          ik = (i + k) % r1->length;
          jk = (j - k + r2->length) % r2->length;
          if (k == minlen || r1->var[ik] != - r2->var[jk])
          { /* sottostringa inversa uguale di lunghezza k */
            if (k > kmax) /* migliore della precedente */
            {
              iopt = i;
              jopt = j;
              kmax = k;
              direction = -1;
            }
            break;
          }
        }
      }
    }
  }

  //printf ("********* kmax = %d; i = %d, j = %d, direction = %d\n", kmax, iopt, jopt, direction);
  retcode = 0;
  if (kmax >= kminforrotation)
  {
    /* perform rotation */
    /* r1: actual start iopt, should be r1->length - kmax */
    /* rotate left iopt - r1->length + kmax times */
    numofrot = (iopt + kmax)%r1->length;
    for (i = 0; i < numofrot; i++)
    {
      saved = r1->var[0];
      for (k = 1; k < r1->length; k++)
        r1->var[k-1] = r1->var[k];
      r1->var[r1->length - 1] = saved;
    }
    retcode = -kmax;
    numofrot = (jopt - kmax + r2->length + 1)%r2->length;
    if (direction > 0)
    {
      retcode = kmax;
      numofrot = (jopt + kmax)%r2->length;
    }
    for (j = 0; j < numofrot; j++)
    {
      saved = r2->var[0];
      for (k = 1; k < r2->length; k++)
        r2->var[k-1] = r2->var[k];
      r2->var[r2->length - 1] = saved;
    }
  }

  return (retcode);
}

int
sp_findduplicatewords (struct presentation *p)
{
  struct presentationrule *r, *s;
  int count = 0;
  int i, equal;

  for (r = p->rules; r; r = r->next)
  {
    if (r->next == 0) break;
    for (s = r->next; s; s = s->next)
    {
      if (r->length == 0) continue;
      if (r->length != s->length) continue;
      equal = 1;
      for (i = 0; i < r->length; i++)
      {
        if (r->var[i] != s->var[i]) {equal = 0; break;}
      }
      if (equal)
      {
        count++;
        s->length = 0;
      }
    }
  }
  return (count);
}

int
sp_removeemptyrules (struct presentation *p)
{
  struct presentationrule *r, *pr;
  int goon = 1;
  int count = 0;

  while (goon)
  {
    goon = 0;
    pr = 0;
    for (r = p->rules; r; r = r->next)
    {
      if (r->length == 0)
      {
        if (pr) pr->next = r->next;
         else p->rules = r->next;
	free (r);
        goon++;
        count++;
        break;
      }
      pr = r;
    }
  }
  return (count);
}

int
sp_eliminatevar (struct presentation *p)
{
  /* find a rule of length one */
  struct presentationrule *r;

  for (r = p->rules; r; r = r->next)
  {
    if (r->length != 1) continue;
    sp_do_eliminatevar (p, r->var[0]);
    return (1);
  }
  return (0);
}

void
sp_do_eliminatevar (struct presentation *p, int generator)
{
  struct presentationrule *r;
  int letter, i, j;

  if (generator < 0) generator *= -1;

  p->gennum--;

  for (r = p->rules; r; r = r->next)
  {
    for (i = 0, j = 0; i < r->length; i++)
    {

      letter = r->var[i];
      if (letter < 0) letter *= -1;
      if (letter == generator) continue;
      if (letter > generator) letter--;
      if (r->var[i] < 0) letter *= -1;
      r->var[j++] = letter;
    }
    r->length = j;
  }
  for (r = p->elements; r; r = r->next)
  {
    for (i = 0, j = 0; i < r->length; i++)
    {
      letter = r->var[i];
      if (letter < 0) letter *= -1;
      if (letter == generator) continue;
      if (letter > generator) letter--;
      if (r->var[i] < 0) letter *= -1;
      r->var[j++] = letter;
    }
    r->length = j;
  }
}

/*
 * print presentation of fundamental group
 */

void
print_presentation (struct presentation *p)
{
  struct presentationrule *r;
  int rulenum, elementsnum, g;
  extern int outformat;

  if (outformat == OUTFORMAT_APPCONTOUR) printf ("fpgroup {\n");
  if (p->gennum == 0)
  {
    if (!quiet)
    {
      printf ("Trivial group");
      if (p->elements)
      {
        printf (" with selected element");
        if (p->elements->next) printf ("s");
      }
      printf ("\n");
    }
    printf ("<");
    if (p->elements)
    {
      printf (";;");
      for (r = p->elements; r; r = r->next)
      {
        printf ("1");
        if (r->next) printf (", ");
      }
    }
    if (quiet) printf (">\n");
     else printf (">\n");
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
    return;
  }
  for (rulenum = 0, r = p->rules; r; r = r->next) rulenum++;
  for (elementsnum = 0, r = p->elements; r; r = r->next) elementsnum++;
  if (!quiet)
  {
    if (rulenum == 0)
    {
      printf ("Free group of rank %d", p->gennum);
      if (elementsnum > 0) printf (" and %d selected element%s", elementsnum, (elementsnum>1)?"s":"");
      printf ("\n");
    } else {
      printf ("Finitely presented group with %d generator%s", p->gennum,
                (p->gennum == 1)?"":"s");
      if (elementsnum > 0) printf (" and %d selected element%s", elementsnum, (elementsnum>1)?"s":"");
      printf ("\n");
    }
  }

  if (p->gennum > MAXGENERATORS)
  {
    printf ("Too many generators to print the rules\n");
    return;
  }
  printf ("<");
  for (g = 0; g < p->gennum; g++)
  {
    if (g > 0) printf (",");
    printf ("%c", 'a' + g);
  }
  printf ("; ");
  print_rule_list (p->rules, p->gennum);
  if (p->elements)
  {
    /* there are selected elements present */
    printf ("; ");
    print_rule_list (p->elements, p->gennum);
    printf (">\n");
  } else printf (">\n");
  if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
}

void
print_rule_list (struct presentationrule *r, int gennum)
{
  if (r == 0) return;

  print_single_rule (r, gennum);
  if (r->next) printf (", ");
  print_rule_list (r->next, gennum);
}

void
print_single_rule (struct presentationrule *r, int gennum)
{
  int i;
  int generator;
  char var;

  if (r->length == 0) printf ("1");
  for (i = 0; i < r->length; i++)
  {
    generator = r->var[i];
    var = 'a';
    if (generator < 0)
    {
      var = 'A';
      generator *= -1;
    }
    generator--;
    var += generator;
    if (generator >= gennum) var = '?';
    printf ("%c", var);
  }
}

/*
 * print invariant factors
 */

void
print_invariant_factors (struct presentation *p)
{
  struct presentationrule *r;
  int sum, ones, rank, matrixrank;
  int i;

  topreabelian (p);

  ones = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    if (abs(sum) == 1) ones++;
    if (sum) matrixrank = i;
  }
  rank = p->gennum - matrixrank;

  if (matrixrank == ones)
  {
    /* torsion free */
    if (quiet) printf ("rank: %d\n", rank);
      else {
      if (rank) printf ("Torsion-free abelian group of rank %d\n", rank);
        else printf ("Trivial group\n");
    }
    return;
  }

  printf ("Abelian group of rank %d with torsion; invariant factors:", rank);

  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    if (abs(sum) == 1) continue;
    if (sum == 0) break;
    printf (" %d", abs(sum));
  }
  printf ("\n");

  return;
}

/*
 * print exponent matrix
 */

void
print_exponent_matrix (struct presentation *p)
{
  struct presentationrule *r;
  int i, sum;

  printf ("Exponent sum matrix:\n");
  for (i = 1; i <= p->gennum; i++)
  {
    printf ("[");
    for (r = p->rules; r; r = r->next)
    {
      sum = get_exp_sum (r, i);
      printf ("%3d", sum);
    }
    printf ("]\n");
  }
}

/*
 * interact with user about finitely presented group
 */ 

/* local prototipes */

struct presentationrule *find_nth_relator (struct presentation *p, int n);
void tietze_invrel (struct presentationrule *r);
void rotate_relator (struct presentationrule *r, int rots);
struct presentationrule *tietze_exchange_relators (struct presentationrule *list,
                         struct presentationrule *r1,
                         struct presentationrule *r2);
void tietze_invgen (struct presentation *p, int n);
void tietze_exchange_generators (struct presentation *p, int m, int n);
void tietze_exchgen_in_rule (struct presentationrule *r, int m, int n);
void tietze_xisabtok (struct presentation *p, int m, int n, int k);
struct presentationrule *tietze_xisabtokl (struct presentationrule *r, int m, int n, int posk, int signk);
void addcommutator (struct presentation *p, int m, int n);

/* user interaction */

void
fg_interactive (struct presentation *p)
{
  struct presentationrule *r, *r1, *r2;
  char commandline[100];
  char *chpt, *cmd;
  int m, n, nsaved, rots;
  int print = 1;
  int autorot = 1;
  int count, kmax, direction;

  assert (p->elements == 0);  /* TODO: take into account the presence of selected elements */
  while (1)
  {
    if (print)
    {
      print_presentation (p);
      print_exponent_matrix (p);
    }
    print = 0;
    printf ("--> ");
    if (fgets (commandline, 99, stdin) == 0) break;
    chpt = commandline;
    while (*chpt && isspace (*chpt)) chpt++;
    cmd = strsep (&chpt, " \t\n");
    if (strcasecmp (cmd, "quit") == 0) break;
    if (strcasecmp (cmd, "help") == 0)
    {
      printf ("Valid commands:\n");
      printf (" quit\n");
      printf (" help\n");
      printf ("Column transformations:\n");
      printf (" negcol <n>\n");
      printf (" xchcol <m> <n>\n");
      printf (" addcol <m> <n> : add col n to col m (right mult)\n");
      printf (" subcol <m> <n> : subtract col n to col m (right mult with inverse)\n");
      printf ("Row transformations:\n");
      printf (" negrow <n>\n");
      printf (" xchrow <m> <n>\n");
      printf (" addrow <m> <n> : add row n to row m\n");
      printf (" subrow <m> <n> : subtract row n to row m\n");
      printf ("Alexander polynomial:\n");
      printf (" alexander : compute Alexander polynomial\n");
      printf (" corank1 : links... (experimental)\n");
      printf ("Relators and simplification:\n");
      printf (" rotrel <n> <rots> : rotate relator <n> left <rots> times\n");
      printf (" simplify : perform a complete simplification (might change signs\n");
      printf (" test r1 r2: find common substring\n");
      printf (" preabelian: transform to preabelian presentation\n");
      printf (" preabelianstep: perform a preabelian step\n");
      printf (" commute <m> <n>: add a commutator for variables m and n\n");
      printf (" auto 0/1: disable/enable automatic rotation before addcol/subcol\n");
      continue;
    }
    if (strcasecmp (cmd, "auto") == 0)
    {
      n = atoi (chpt);
      if (n == 0 && autorot == 0) {printf ("Auto is already off\n"); continue;}
      if (n == 1 && autorot == 1) {printf ("Auto is already on\n"); continue;}
      if (n < 0 || n > 1) {printf ("Invalid value %d\n", n); continue;}
      autorot = n;
    }
    if (strcasecmp (cmd, "negcol") == 0)
    {
      n = atoi (chpt);
      r = find_nth_relator(p, n);
      if (r == 0) {printf ("Invalid column number %d\n", n); continue;}
      tietze_invrel (r);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "rotrel") == 0)
    {
      n = (int) strtol (chpt, &chpt, 10);
      rots = (int) strtol (chpt, &chpt, 10);
      r = find_nth_relator(p, n);
      if (r == 0) {printf ("Invalid column number %d\n", n); continue;}
      rotate_relator (r, rots);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "addcol") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      r1 = find_nth_relator(p, m);
      r2 = find_nth_relator(p, n);
      if (r1 == 0 || r2 == 0) {printf ("Invalid column number %d or %d\n", m, n); continue;}
      if (r1 == r2) {printf ("Same column %d = %d\n", m, n); continue;}
      if (autorot)
      {
        tietze_invrel (r2);
        kmax = sp_fcs_r1r2 (r1, r2, 0, 1);
        tietze_invrel (r2);
        printf ("Found common substring of length %d\n", kmax);
        //print_presentation (p);
      }
      tietze_mulright (p, r1, r2, 1);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "subcol") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      r1 = find_nth_relator(p, m);
      r2 = find_nth_relator(p, n);
      if (r1 == 0 || r2 == 0) {printf ("Invalid column number %d or %d\n", m, n); continue;}
      if (r1 == r2) {printf ("Same column %d = %d\n", m, n); continue;}
      if (autorot)
      {
        kmax = sp_fcs_r1r2 (r1, r2, 0, 1);
        printf ("Found common substring of length %d\n", kmax);
      }
      tietze_mulright (p, r1, r2, -1);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "xchcol") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      if (n < m) {nsaved = n; n = m; m = nsaved;}
      r1 = find_nth_relator(p, m);
      r2 = find_nth_relator(p, n);
      if (r1 == 0 || r2 == 0) {printf ("Invalid column number %d or %d\n", m, n); continue;}
      if (r1 == r2) {printf ("Same column %d = %d\n", m, n); continue;}
      p->rules = tietze_exchange_relators (p->rules, r1, r2);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "negrow") == 0)
    {
      n = atoi (chpt);
      if (n <= 0 || n > p->gennum) {printf ("Invalid row number %d\n", n); continue;}
      tietze_invgen (p, n);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "xchrow") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      if (m <= 0 || m > p->gennum) {printf ("Invalid row number %d\n", m); continue;}
      if (n <= 0 || n > p->gennum) {printf ("Invalid row number %d\n", n); continue;}
      if (m == n) {printf ("Same row %d = %d\n", m, n); continue;}
      tietze_exchange_generators (p, m, n);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "addrow") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      if (m <= 0 || m > p->gennum) {printf ("Invalid row number %d\n", m); continue;}
      if (n <= 0 || n > p->gennum) {printf ("Invalid row number %d\n", n); continue;}
      if (m == n) {printf ("Same row %d = %d\n", m, n); continue;}
      tietze_xisabtok (p, n, m, -1);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "subrow") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      if (m <= 0 || m > p->gennum) {printf ("Invalid row number %d\n", m); continue;}
      if (n <= 0 || n > p->gennum) {printf ("Invalid row number %d\n", n); continue;}
      if (m == n) {printf ("Same row %d = %d\n", m, n); continue;}
      tietze_xisabtok (p, n, m, 1);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "preabelianstep") == 0)
    {
      count = preabelian_step (p, 1, p->rules);
      if (count)
      {
        printf ("%d modifications\n", count);
        print = 1;
      } else printf ("Nothing done\n");
      continue;
    }
    if (strcasecmp (cmd, "preabelian") == 0)
    {
      topreabelian (p);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "alexander") == 0)
    {
      alexander (p);
      continue;
    }
    if (strcasecmp (cmd, "corank1") == 0)
    {
      corank_one_alexander (p);
      continue;
    }
    if (strcasecmp (cmd, "commute") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      if (m <= 0 || m > p->gennum) {printf ("Invalid row number %d\n", m); continue;}
      if (n <= 0 || n > p->gennum) {printf ("Invalid row number %d\n", n); continue;}
      if (m == n) {printf ("Same row %d = %d\n", m, n); continue;}
      addcommutator (p, m, n);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "read") == 0)
    {
      read_group_presentation (stdin, p);
      print = 1;
      continue;
    }
    print = 1;
    if (strcasecmp (cmd, "simplify") == 0)
    {
      simplify_presentation (p);
      continue;
    }
    if (strcasecmp (cmd, "common") == 0)
    {
      sp_findcommonsubstring (p, 1);
      continue;
    }
    if (strcasecmp (cmd, "test") == 0)
    {
      m = (int) strtol (chpt, &chpt, 10);
      n = (int) strtol (chpt, &chpt, 10);
      r1 = find_nth_relator(p, m);
      r2 = find_nth_relator(p, n);
      if (r1 == 0 || r2 == 0) {printf ("Invalid column number %d or %d\n", m, n); continue;}
      if (r1 == r2) {printf ("Same column %d = %d\n", m, n); continue;}
      direction = 1;
      kmax = sp_fcs_r1r2 (r1, r2, 1, 1);
      if (kmax < 0)
      {
        direction = -1;
        kmax = -kmax;
      }
      if (kmax) printf ("found common subsequence of length %d, direction = %d\n", kmax, direction);
      continue;
    }
    print = 0;
    if (*cmd == 0) continue;
    if (*cmd == EOF) break;
    printf ("Invalid command: %s\n", cmd);
  }
}

/*
 * find relator by column number
 */

struct presentationrule *
find_nth_relator (struct presentation *p, int n)
{
  int i;
  struct presentationrule *r;

  if (n <= 0) return (0);
  for (i = 1, r = p->rules; i < n && r; i++, r = r->next);
  return (r);
}

void
tietze_invrel (struct presentationrule *r)
{
  int i, j, saved;

  if (r->length == 0) return;
  for (i = 0, j = r->length - 1; i <= j; i++, j--)
  {
    if (i == j) {r->var[i] = -r->var[i]; continue;}
    saved = r->var[i];
    r->var[i] = -r->var[j];
    r->var[j] = -saved;
  }
}

void
rotate_relator (struct presentationrule *r, int rots)
{
  int i, j;
  int saved;

  if (r->length <= 1) return;
  while (rots < 0) rots += r->length;
  if (rots == 0) return;
  for (i = 0; i < rots; i++)
  {
    saved = r->var[0];
    for (j = 1; j < r->length; j++)
    {
      r->var[j-1] = r->var[j];
    }
    r->var[r->length - 1] = saved;
  }
}

struct presentationrule *
tietze_exchange_relators (struct presentationrule *list,
                           struct presentationrule *r1,
                           struct presentationrule *r2)
{
  struct presentationrule *r, *r1next, *r2next;
  assert (list != r2);
  assert (list && r1 && r2);

  if (list != r1)
  {
    list->next = tietze_exchange_relators (list->next, r1, r2);
    return (list);
  }

  /* r1 is the first element */
  r1next = r1->next;
  r2next = r2->next;

  if (r1next == r2)
  {
    r1->next = r2next;
    r2->next = r1;
    return (r2);
  }
  for (r = r1; r; r = r->next)
  {
    if (r->next == r2)
    {
      r->next = r1;
      break;
    }
  }
  r1->next = r2next;
  r2->next = r1next;

  return (r2);
}

void
tietze_mulright (struct presentation *p,
                  struct presentationrule *r1,
                  struct presentationrule *r2,
                  int expon)
{
  int i, j, k, totlen;
  struct presentationrule *r, *rr;

  assert (p && r1 && r2);
  if (expon == 0) return;

  totlen = r1->length + abs(expon)*r2->length;

  r = (struct presentationrule *) malloc (sizeof (struct presentationrule) + totlen*sizeof(int));
  r->next = r1->next;
  r->length = totlen;

  for (i = 0; i < r1->length; i++)
  {
    r->var[i] = r1->var[i];
  }
  assert (i == r1->length);

  for (k = 0; k < abs(expon); k++)
  {
    for (j = 0; j < r2->length; j++)
    {
      if (expon > 0)
      {
        r->var[i++] = r2->var[j];
      } else {
        r->var[i++] = - r2->var[r2->length - j - 1];
      }
    }
  }
  assert (i == r->length);

  r->next = r1->next;
  if (p->rules == r1)
  {
    p->rules = r;
  } else {
    for (rr = p->rules; rr; rr = rr->next)
    {
      if (rr->next == r1)
      {
        rr->next = r;
        break;
      }
    }
  }

  free (r1);
  sp_simplifyword (p);
  return;
}

void
tietze_invgen (struct presentation *p, int g)
{
  struct presentationrule *r;
  int i;

  assert (g > 0 && g <= p->gennum);

  for (r = p->rules; r; r = r->next)
  {
    for (i = 0; i < r->length; i++)
    {
      if (abs(r->var[i]) == g) r->var[i] = - r->var[i];
    }
  }
  for (r = p->elements; r; r = r->next)
  {
    for (i = 0; i < r->length; i++)
    {
      if (abs(r->var[i]) == g) r->var[i] = - r->var[i];
    }
  }
}

void
tietze_exchange_generators (struct presentation *p, int m, int n)
{
  struct presentationrule *r;

  assert (m > 0 && m <= p->gennum);
  assert (n > 0 && n <= p->gennum);

  if (m == n) return;

  for (r = p->rules; r; r = r->next) tietze_exchgen_in_rule (r, m, n);
  for (r = p->elements; r; r = r->next) tietze_exchgen_in_rule (r, m, n);
}

void
tietze_exchgen_in_rule (struct presentationrule *r, int m, int n)
{
  int i;

  for (i = 0; i < r->length; i++)
  {
    if (abs(r->var[i]) == m)
    {
      if (r->var[i] > 0) r->var[i] = n;
        else r->var[i] = -n;
      continue;
    }
    if (abs(r->var[i]) == n)
    {
      if (r->var[i] > 0) r->var[i] = m;
        else r->var[i] = -m;
    }
  }
}

void
tietze_xisabtok (struct presentation *p, int m, int n, int k)
{
  int posk, signk;

  assert (m > 0 && m <= p->gennum);
  assert (n > 0 && n <= p->gennum);
  assert (m != n);

  if (k == 0) return;

  signk = 1;
  posk = k;
  if (k < 0) {signk = -1; posk = -k;}

  p->rules = tietze_xisabtokl (p->rules, m, n, posk, signk);
  p->elements = tietze_xisabtokl (p->elements, m, n, posk, signk);

  sp_simplifyword (p);

  return;
}

struct presentationrule *
tietze_xisabtokl (struct presentationrule *r, int m, int n, int posk, int signk)
{
  int i, j, kk, numvarm;
  struct presentationrule *newr;

  assert (posk > 0);
  assert (signk == 1 || signk == -1);
  if (r == 0) return (0);

  /* devo calcolare la nuova lunghezza
   * che aumenta di k*s dove s e' il numero di
   * occorrenze della variabile m
   */

  numvarm = 0;
  for (i = 0; i < r->length; i++)
  {
    if (abs(r->var[i]) == m) numvarm++;
  }
  if (numvarm == 0)
  {
    r->next = tietze_xisabtokl (r->next, m, n, posk, signk);
    return (r);
  }
  newr = (struct presentationrule *) malloc (sizeof (struct presentationrule)
                                       + sizeof(int)*(r->length + posk*numvarm));

  for (i = 0, j = 0; i < r->length; i++)
  {
    if (abs(r->var[i]) != m)
    {
      newr->var[j++] = r->var[i];
      continue;
    }
    if (r->var[i] > 0)
    {
      newr->var[j++] = r->var[i];
      for (kk = 0; kk < posk; kk++)
      {
        newr->var[j++] = - signk*n;
      }
    } else {
      for (kk = 0; kk < posk; kk++)
      {
        newr->var[j++] = signk*n;
      }
      newr->var[j++] = r->var[i];
    }
  }
  newr->length = j;
  assert (j == r->length + posk*numvarm);

  newr->next = r->next;
  r->next = 0;
  free (r);
  newr->next = tietze_xisabtokl (newr->next, m, n, posk, signk);
  return (newr);
}

/*
 * transformation of a presentation to a preabelian form
 */

void
topreabelian (struct presentation *p)
{
  if (verbose) printf ("Computing preabelian presentation...\n");
  assert (p->elements == 0);   /* TODO: take into account selected elements */
  while (preabelian_step (p, 1, p->rules));
}

int
preabelian_step (struct presentation *p, int row, struct presentationrule *rcol)
{
  int sum, pivot, optval, optrow, nmulrow;
  int icol, optcol, nmulcol;
  int i, isdivisor, divisor;
  struct presentationrule *r, *optrcol, *nmulrcol;

  if (row > p->gennum || rcol == 0) return (0);
  /* scan the exponent submatrix */

  pivot = get_exp_sum (rcol, row);
  optval = abs(pivot);
  isdivisor = 1;
  nmulrow = optrow = 0;
  nmulrcol = optrcol = 0;
  optcol = nmulcol = 0;

  for (i = row; i <= p->gennum; i++)
  {
    for (r = rcol, icol = row; r; r = r->next, icol++)
    {
      sum = get_exp_sum (r, i);
      if (sum)
      {
        if (optval == 0 || abs(sum) < optval)
        {
          optval = abs(sum);
          optrow = i;
          optrcol = r;
          optcol = icol;
        }
        if ( (abs(pivot) > 0 && (abs(sum) % abs(pivot)) != 0) ||
             (abs(pivot) == 0 && sum != 0) )
        { /* element is not a multiple of pivot */
          isdivisor = 0;
          nmulrow = i;
          nmulrcol = r;
          nmulcol = icol;
        }
      }
    }
  }
  if (optrow)
  {
    assert (optrcol);
    assert (optrow != row || optrcol != rcol);
    if (interactive || verbose)
    {
      printf ("Found optimal value at row %d, col %d\n", optrow, optcol);
      if (optrow != row) printf ("  Exchange row %d with row %d\n", row, optrow);
      if (optrcol != rcol) printf ("  Exchange col %d with col %d\n", row, optcol);
    }
    /* exchange rows and columns */
    if (optrow != row) tietze_exchange_generators (p, row, optrow);
    if (optrcol != rcol) p->rules = tietze_exchange_relators (p->rules, rcol, optrcol);
    return (1);
  }
  /* at this point we have the minimum value at the pivot position */
  /* working on the first column... */
  for (i = row+1; i <= p->gennum; i++)
  {
    sum = get_exp_sum (rcol, i);
    if (sum)
    {
      assert (pivot);
      /* Euclid's division between sum and pivot */
      divisor = abs(sum)/abs(pivot);
      assert (divisor);
      /* if same sign, then subtract, else add */
      if ((sum > 0 && pivot > 0) || (sum < 0 && pivot < 0)) divisor = -divisor;
      tietze_xisabtok (p, row, i, -divisor);
      if (interactive || verbose) printf ("  Row %d times %d added to row %d\n", row, divisor, i);
      return (1);
    }
  }
  /* if we are here, the first column is zero */
  /* working on the first row... */
  for (r = rcol->next, icol = row + 1; r; r = r->next, icol++)
  {
    sum = get_exp_sum (r, row);
    if (sum)
    {
      assert (pivot);
      /* Euclid's division between sum and pivot */
      divisor = abs(sum)/abs(pivot);
      assert (divisor);
      /* if same sign, then subtract, else add */
      if ((sum > 0 && pivot > 0) || (sum < 0 && pivot < 0)) divisor = -divisor;
      if (divisor > 0) tietze_invrel (rcol);
      sp_fcs_r1r2 (r, rcol, 0, 1);
      if (divisor > 0) tietze_invrel (rcol);
      tietze_mulright (p, r, rcol, divisor);
      if (interactive || verbose) printf ("  Column %d times %d added to column %d\n", row, divisor, icol);
      return (1);
    }
  }
  /* if we are here, the first column and first row are zero */
  if (isdivisor == 0)
  {
    assert (nmulrcol);
    /* add column nmulrcol to first column */
    tietze_invrel (nmulrcol);
    sp_fcs_r1r2 (rcol, nmulrcol, 0, 1);
    tietze_invrel (nmulrcol);
    tietze_mulright (p, rcol, nmulrcol, 1);
    if (interactive || verbose) printf ("Found nonmultiple value at row %d, column %p\n", nmulrow, nmulrcol);
    return (1);
  }
  /* if we are here, the first column and first row are zero and all other elements are multiple */
  if (pivot < 0)
  {
    tietze_invgen (p, row);
    if (interactive || verbose) printf ("  Inverted generator %d\n", row);
    return (1);
  }
  return (preabelian_step (p, row + 1, rcol->next));
}

int
get_exp_sum (struct presentationrule *r, int n)
{
  int sum = 0;
  int j;

  for (j = 0; j < r->length; j++)
  {
    if (abs(r->var[j]) != n) continue;
    if (r->var[j] > 0) sum++; else sum--;
  }

  return (sum);
}

/*
 * compute the rank of the fundamental group
 * (it is assumed that the presentation is in
 * preabelian form!!!)
 */

int
compute_fg_rank (struct presentation *p)
{
  struct presentationrule *r;
  int i, ones = 0, zeros = 0, dim = 0, divisor;

  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    divisor = get_exp_sum (r, i);
    assert (divisor >= 0);
    dim++;
    if (divisor == 1) ones++;
    if (divisor == 0) zeros++;
    if (divisor == 1) assert (zeros == 0);
  }

  if (ones + zeros != dim) return (-1);

  return (p->gennum - ones);
}

void
addcommutator (struct presentation *p, int m, int n)
{
  struct presentationrule *r, *rr;

  assert (m > 0 && m <= p->gennum);
  assert (n > 0 && n <= p->gennum);
  if (m == n) return;

  r = (struct presentationrule *) malloc (sizeof (struct presentationrule) +
        4*sizeof(int));

  r->length = 4;
  r->var[0] = m;
  r->var[1] = n;
  r->var[2] = -m;
  r->var[3] = -n;
  r->next = 0;

  if (p->rules == 0)
  {
    p->rules = r;
    return;
  }
  for (rr = p->rules; rr; rr = rr->next)
  {
    if (rr->next == 0)
    {
      rr->next = r;
      return;
    }
  }
}

/*
 * reading a group presentation from stdin
 */

/* local prototypes */
struct presentationrule *read_relators_list (FILE *file, char *generator_names, int gennum);

struct presentationlist *
read_group_presentation_list (FILE *file)
{
  struct presentationlist *pstlist;
  struct presentationlist *pstl;
  struct presentationlist *pstnew;
  struct presentation *pst;
  int tok;

  tok = gettoken (file);
  assert (tok == TOK_FPGROUPLIST);
  tok = gettoken (file);
  assert (tok == TOK_LBRACE);

  pstlist = pstl = 0;
  while (1)
  {
    tok = gettoken (file);
    ungettoken (tok);
    if (tok != KEY_LT) break;

    pst = (struct presentation *) malloc (sizeof (struct presentation));
    pst->rules = pst->elements = 0;
    read_group_presentation (file, pst);
    pstnew = (struct presentationlist *) malloc (sizeof (struct presentationlist));
    pstnew->p = pst;
    pstnew->next = 0;
    if (pstl)
    {
      pstl->next = pstnew;
      pstl = pstnew;
    } else {
      pstlist = pstl = pstnew;
    }
  }

  tok = gettoken (file);
  assert (tok == TOK_RBRACE);

  return (pstlist);
}

void
read_group_presentation (FILE *file, struct presentation *p)
{
  int tok;
  int wantrightbrace = 0;
  char generator_names[MAXGENERATORS];

  remove_all_relators (p);
  p->gennum = 0;

  tok = gettoken (file);

  if (tok == TOK_FPGROUP)
  {
    tok = gettoken (file);
    if (tok == TOK_LBRACE)
    {
      wantrightbrace = 1;
      tok = gettoken (file);
    }
  }
  if (tok != KEY_LT)
  {
    printf ("'<' expected, got token %d instead\n", tok);
    return;
  }

  p->gennum = read_generators_list (file, generator_names, MAXGENERATORS);

  tok = gettoken (file);
  if (tok != TOK_SEMICOLON && tok != KEY_GT)
  {
    printf ("Invalid syntax, tok = %d; semicolon expected\n", tok);
    return;
  }
  if (tok == KEY_GT) ungettoken (tok);
  p->rules = read_relators_list (file, generator_names, p->gennum);

  tok = gettoken (file);
  if (tok == TOK_SEMICOLON)
  {
    if (verbose > 1) printf ("Reading selected elements in group\n");
    p->elements = read_relators_list (file, generator_names, p->gennum);
    tok = gettoken (file);
  }

  if (tok != KEY_GT)
  {
    printf ("Invalid syntax, tok = %d; rangle expected\n", tok);
    return;
  }

  if (wantrightbrace)
  {
    tok = gettoken (file);
    if (tok != TOK_RBRACE)
    {
      printf ("Expected right brace at end\n");
      ungettoken (tok);
    }
  }
  return;
}

int
read_generators_list (FILE *file, char *gennames, int maxgennum)
{
  char ch;

  ch = mygetchar (file);

  if (isalnum(ch))
  {
    if (! islower(ch))
    {
      printf ("Generators must be a single lower-case letter\n");
      return (0);
    }
    *gennames++ = ch;
    ch = mygetchar (file);
    if (ch == ',') return (read_generators_list (file, gennames, maxgennum-1) + 1);
    ungetc (ch, file);
    return (1);
  }
  ungetc (ch, file);
  return (0);
}

#define BUFSIZE 2000

struct presentationrule *
read_relators_list (FILE *file, char *generator_names, int gennum)
{
  int i, j;
  int sign;
  char ch, buf[BUFSIZE + 1];
  struct presentationrule *r;

  ch = mygetchar (file);

  if (isalpha(ch))
  {
    i = 0;
    buf[i++] = ch;
    while (isalpha(ch = mygetchar (file)))
      if (i < BUFSIZE) buf[i++] = ch;
    buf[i] = 0;
    ungetc (ch, file);
    if (i >= BUFSIZE) fprintf (stderr, "Relator is too long, truncated!\n");
    r = (struct presentationrule *) malloc (sizeof (struct presentationrule) +
           strlen(buf)*sizeof(int));
    r->length = strlen(buf);
    r->next = 0;
    for (i = 0; i < r->length; i++)
    {
      sign = 1;
      if (buf[i] >= 'A' && buf[i] <= 'Z')
      {
        sign = -1;
        buf[i] += 'a' - 'A';
      }
      r->var[i] = 0;  //XXX per ora!
      for (j = 0; j < gennum; j++)
      {
        if (buf[i] == generator_names[j]) r->var[i] = sign*(j + 1);
      }
      if (r->var[i] == 0)
      {
        r->length = 0;
        printf ("Cannot find generator: %c\n", buf[i]);
      }
    }
    ch = mygetchar (file);
    if (ch == ',')
    {
      r->next = read_relators_list (file, generator_names, gennum);
      return (r);
    }
    ungetc (ch, file);
    return (r);
  }
  ungetc (ch, file);

  return (0);
}

void
remove_all_relators (struct presentation *p)
{
  struct presentationrule *r;

  r = p->rules;
  if (r == 0) return;

  p->rules = r->next;
  free (r);
  remove_all_relators (p);
}

/*
 * procedures for the simplification of the cell complex
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
 * perform varius melt operations iteratively
 */

int
complex_melt (struct ccomplex *cc)
{
  int goon = 1;
  int count = 0;

  while (goon)
  {
    goon = 0;
    goon += complex_facemelt (cc);
    goon += complex_faceremovekink (cc);
    count += goon;
  }
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
      face1 = faces + nface1;
      face2 = faces + nface2;
      ivec1 = face1->faceborder;
      ivec2 = face2->faceborder;
      if (face1->facebordernum + face2->facebordernum <= 2) continue;
      /*
       * We can melt the two faces together
       */
      if (ivec1[i1] < 0) cc_revert_face (cc, nface1);
      if (ivec2[i2] > 0) cc_revert_face (cc, nface2);
      goon = 1;
      count++;
      complex_do_melt_faces (cc, nface1, nface2, n);
      if (debug) cellcomplex_checkconsistency (cc);
    }
  }
  if (debug) printf ("Melted %d faces\n", count);
  return (count);
}

/*
 * remove consecutive equal arcs in a face
 * (with opposite orientation), the face must
 * have at least three arcs.
 * the case with two arcs creates a "bubble"
 * that can be removed without changing the
 * fundamental group.
 */

int
complex_faceremovekink (struct ccomplex *cc)
{
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexface *faces = cc->faces;
  struct ccomplexarc *arc;
  struct ccomplexface *face;
  struct ccomplexnode *nodes = cc->nodes;
  int count = 0;
  int goon = 1;
  int n, i;
  int *ivec;
  int nnode, inext;

  while (goon)
  {
    goon = 0;
    for (n = 0; n < cc->facedim; n++)
    {
      face = faces + n;
      if (face->type == CC_REMOVED) continue;
      if (face->facebordernum <= 2) continue;
      ivec = face->faceborder;
      for (i = 0; i < face->facebordernum; i++)
      {
        arc = arcs + onarc2narc (ivec[i]);
        if (arc->refcount != 2) continue;
        nnode = arc->endb;
        if (ivec[i] < 0) nnode = arc->enda;
        if (nodes[nnode].refcount != 1) continue;
        inext = i + 1;
        if (inext >= face->facebordernum) inext = 0;
        if (ivec[i] + ivec[inext] == 0)
        {
          goon = 1;
          count++;
          complex_do_removekink (cc, n, i, inext, nnode);
          if (debug) cellcomplex_checkconsistency (cc);
          break;
        }
      }
    }
  }

  if (debug) printf ("Removed %d face kinks\n", count);
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
  struct ccomplexarc *arc;
  int *ivec1, *ivec2, *newivec;
  int found1 = 0;
  int found2 = 0;
  int i, j, k, j2;
  int newsize;

  arc = cc->arcs + narc;
  face1 = cc->faces + nface1;
  ivec1 = face1->faceborder;
  face2 = cc->faces + nface2;
  ivec2 = face2->faceborder;
  newsize = face1->facebordernum + face2->facebordernum - 2;
  newivec = (int *) malloc (newsize*sizeof (int));

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
  assert (k == newsize);
  arc->refcount -= 2;
  complex_remove_arc (cc, narc);
  face1->facebordernum += face2->facebordernum;
  face1->facebordernum -= 2;
  free (ivec1);
  face1->faceborder = newivec;
  complex_remove_face_nd (cc, nface2);
}

void
complex_do_removekink (struct ccomplex *cc, int n, int iprev, int inext, int nnode)
{
  struct ccomplexface *faces = cc->faces;
  struct ccomplexface *face;
  struct ccomplexarc *arcs = cc->arcs;
  struct ccomplexarc *arc;
  int i, k, *ivec, narc;

  face = faces + n;
  ivec = face->faceborder;
  narc = onarc2narc (ivec[iprev]);
  arc = arcs + narc;
  arc->refcount -= 2;

  for (i = 0, k = 0; i < face->facebordernum; i++)
  {
    if (i == iprev || i == inext) continue;
    ivec[k++] = ivec[i];
  }
  face->facebordernum -= 2;
  complex_remove_arc (cc, narc);
  complex_remove_node (cc, nnode);
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
complex_remove_face_nd (struct ccomplex *cc, int nface)
{
  struct ccomplexface *faces = cc->faces;

  free (faces[nface].faceborder);
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
 *
 * It is defined up to homotopy type (in particular, 3D cells are all
 * retracted suitably to end up with a 2D complex.
 *
 * If global variable "focus_on_fundamental" is defined, then we are allowed
 * to first apply surgeries that do not change the fundamental group
 */

struct ccomplex *
compute_cellcomplex (struct sketch *s, int fg_type)
{
  struct ccomplex *cc;
  struct region *region;
  int euler, surfeuler, realeuler;
  int i, res, status, stratum;

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

  surfeuler = euler_characteristic (s);
  assert ((surfeuler % 2) == 0);

  /*
   * the cell complex of the external body is simply computed by
   * placing everything into a big sphere
   * NOTE: this modifies the "sketch" structure
   */

  if (fg_type == FG_EXTERNAL)
  {
    status = put_in_s1 (s);
    assert (status);
    postprocesssketch (s);
    surfeuler += 2;
  }
  if (globals.finfinity != 0) fprintf (stderr, "Value of f at infinity (%d) must be zero\n",
                              globals.finfinity);
  assert (globals.finfinity == 0);
  computefvalue (s, s->regions, 0 /* should be finfinity */);
  if (fg_type != FG_SURFACE && globals.focus_on_fundamental && globals.autosurgery)
  {
    while (suggest_p_surgery (s, &region, &stratum))
    {
      res = add_s1 (s, region, stratum, -1);
      surfeuler -= 2;
      if (verbose) printf ("Applying punchhole surgery on region %d, strata %d - %d\n",
                            region->tag, stratum, stratum+1);
      free_connected_components (s);
      assert (res);
    }
  }

  switch (fg_type)
  {
    case FG_SURFACE:
      realeuler = surfeuler;
      break;
    case FG_INTERNAL:
    case FG_EXTERNAL:
      realeuler = surfeuler/2;
      break;
    default:
      realeuler = -9999;
      break;
  }
  cc = (struct ccomplex *) malloc (sizeof (struct ccomplex));
  cc->type = fg_type;
  cc->sketch = s;
  cc->cc = 0;
  cc->cc_characteristics = 0;
  if (s->isempty)
  {
    cc->nodenum = cc->arcnum = cc->facenum = 0;
    cc->nodes = 0;
    cc->arcs = 0;
    cc->faces = 0;
    return (0);
  }
  cc_euler_characteristic (s);
  cc->cc_characteristics = (int *) malloc (s->ccnum * sizeof(int));
  for (i = 0; i < s->ccnum; i++) cc->cc_characteristics[i] = s->cc_characteristics[i];
  cc->surfccnum = s->ccnum;
  assert (s->cc_tagged || s->isempty);

  cc->nodenum = cc->nodedim = fundamental_countnodes (s);
  cc->arcnum = cc->arcdim = fundamental_countarcs (cc->sketch, cc->type);
  cc->facenum = cc->facedim = fundamental_countfaces (cc->sketch, cc->type);
  euler = cc->nodenum - cc->arcnum + cc->facenum;
  if (debug) printf ("Euler characteristic: %d = %d nodes - %d arcs + %d faces.\n",
             euler, cc->nodenum, cc->arcnum, cc->facenum);
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
  /*
   * it is necessary to call find_spanning_tree before any simplification
   * in order to gather information about the cavities
   */
  find_spanning_tree (cc);

  if (debug)
  {
    cellcomplex_print (cc, 2);
    status = cellcomplex_checkconsistency (cc);
    assert (status);
    printf ("Connected components: %d\n", cc->ccnum);
  }
  return (cc);
}

/*
 * find a spanning tree for the 1D cell complex.  As a byproduct
 * we also have a separation of the various connected components
 * of the 1D skeleton, for each connected component a base node is
 * selected.
 * Moreover, for each component the number of adjacent surfaces of the contour
 * is computed, whence the betti number n. 2
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
  int computebetti = 0;
  int *surfcomptag;

  /*
   * this is a little messy, however
   * finding the connected components produces a spanning tree as a 
   * side effect, so we do that here anyway.  However we do not want to
   * lose the betti numbers informations possibly already computed
   */
  /* initializations */
  if (cc->cc)
  {
    /* update che ccomplexcc list to point to still existing nodes
     * in the same component
     */
    for (cccc = cc->cc; cccc; cccc = cccc->next)
    {
      for (i = 0; i < cc->nodedim; i++)
      {
        if (nodes[i].type == CC_REMOVED) continue;
        if (nodes[i].cc == cccc)
        {
          cccc->basenode = i;
          break;
        }
      }
    }
  } else {
    computebetti = 1;
    cc->cc = 0;
    cc->ccnum = 0;
    assert (cc->surfccnum > 0);
    surfcomptag = (int *) malloc (cc->surfccnum * sizeof (int));
  }
  for (i = 0; i < cc->nodedim; i++) nodes[i].cc = 0;
  for (i = 0; i < cc->arcdim; i++) arcs[i].isinspanningtree = 0;

  cccc = cc->cc;
  while ( (bn = (computebetti ?
                 cc_find_new_base_node (cc) :
                 (cccc?cccc->basenode:(-1))
                                                 )) >= 0)
  {
    if (computebetti)
    {
      for (i = 0; i < cc->surfccnum; i++) surfcomptag[i] = 0;
      cccc = (struct ccomplexcc *) malloc (sizeof (struct ccomplexcc));
      cccc->tag = cc->ccnum;
      cccc->basenode = bn;
      assert (nodes[bn].surfcc >= 0);
      surfcomptag[nodes[bn].surfcc]++;
      cccc->p = 0;
      cccc->next = 0;
      cc->ccnum++;
      if (lastcc == 0)
        cc->cc = cccc;
       else
        lastcc->next = cccc;
      lastcc = cccc;
    }
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
          if (computebetti) surfcomptag[nodes[n2].surfcc]++;
        } else {
          assert (nodes[n1].cc == 0);
          nodes[n1].cc = nodes[n2].cc;
          if (computebetti) surfcomptag[nodes[n1].surfcc]++;
        }
      }
    }
    if (computebetti)
    {
      cccc->betti_2 = cccc->spherical_voids = 0;
      for (i = 0; i < cc->surfccnum; i++)
      {
        if (surfcomptag[i] > 0)
        {
          cccc->betti_2++;
          if (cc->cc_characteristics[i] == 2) cccc->spherical_voids++;
        }
      }
      /*
       * the betti number 2 is the number of cavities, i.e. the number of
       * adjacent surfaces minus 1
       */
      cccc->betti_2--;
    } else {
      cccc = cccc->next;
    }
  }
  if (computebetti) free (surfcomptag);
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
  struct region *rleft, *rright;
  int fleft, strata, stratum;
  int vecdim = cc->nodedim;
  int surfccid, dne, dse, dmin, dmax, delta, i;

  vecpt = vec;
  assert (s->cc_tagged);

  for (a = s->arcs; a; a = a->next)
  {
    if (a->cusps > 0)
    {
      rright = a->regionright->border->region;
      surfccid = rright->strati[a->depths[1]]; /* tutte le cuspidi stanno forzatamente nella stessa superficie */
      for (i = 1; i <= a->cusps; i++)
      {
        assert (vecpt < vec + vecdim);
        vecpt->type = CC_NODETYPE_CUSP;
        vecpt->stratum = a->depths[i];
        vecpt->ne = vecpt->se = a;
        vecpt->cusp = i;
        vecpt->surfcc = surfccid;
        vecpt++;
      }
    }
    if (a->endpoints == 0)
    {
      rleft = a->regionleft->border->region;
      fleft = rleft->f;
      strata = fleft -1;
      for (stratum = 0; stratum < strata; stratum++)
      {
        assert (vecpt < vec + vecdim);
        vecpt->type = CC_NODETYPE_VIRTUALCUT;
        if (stratum == a->depths[0]) vecpt->type = CC_NODETYPE_VIRTUALFOLD;
        vecpt->stratum = stratum;
        vecpt->ne = vecpt->se = a;
        if (stratum <= a->depths[0])
          vecpt->surfcc = rleft->strati[stratum];
         else
          vecpt->surfcc = rleft->strati[stratum + 1];
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
      rleft = ane->regionleft->border->region;
      strata = ane->regionright->border->region->f;
      dne = ane->depths[0];
      dse = ase->depths[0];
      if (dne >= dse + 2)
      {
        dne--;
        dmin = dse;
        dmax = dne;
      } else {
        dse++;
        dmin = dne;
        dmax = dse;
      }
      for (stratum = 0; stratum < strata; stratum++)
      {
        delta = 0;
        if (stratum >= dmin) delta++;
        if (stratum >= dmax) delta++;
        assert (vecpt < vec + vecdim);
        vecpt->type = CC_NODETYPE_CUT;
        if (stratum == dse || stratum == dne) vecpt->type = CC_NODETYPE_FOLD;
        vecpt->stratum = stratum;
        vecpt->ne = ane;
        vecpt->se = ase;
        vecpt->surfcc = rleft->strati[stratum + delta];
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
    } else {
      for (stratum = 0; stratum < strata; stratum++)
      {
        if ((stratum % 2) == 1 && cc->type != FG_SURFACE) continue;
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
  /*
   * the target point for a virtual arc is by convention the
   * starting point for the arc corresponding to bord (using
   * the arc orientation, not the "bord" list orientation)
   * this choice has to be compensated when generating the
   * horizontal face in case bord->orientation is negative
   */

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
    columns -= acnodes;
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
    } else {
      if (fg_type != FG_SURFACE) strata /= 2;
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

      default:
        orizfaces += strata/2;
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
        fright += 1;
      }
      for (i = 1; i <= a->cusps; i++)
      {
        dmod2 = a->depths[i] % 2;
        if (dmod2 == 1) cuspwalls++;
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
  int stratum, strata, arcnum, i, dcusp;
  //int d;
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
    for (stratum = 0; stratum < strata; stratum++)
    {
      if ((stratum % 2) == 1 && cc->type != FG_SURFACE) continue;
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
                assert (r->f > 0);     // be sure we cannot be here
                na = fund_findarc (cc, a, astratum, i+1, startcusp);
                assert (na >= 0);
                ivec[arcnum++] = -na-1;
                astratum = stratum;
                if (stratum >= a->depths[i]) astratum++;
                startcusp = i+1;
              }
            }
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
   * only for FG_INTERNAL / EXTERNAL
   */
  if (cc->type != FG_SURFACE)
  {
    for (a = s->arcs; a; a = a->next)
    {
      strata = a->regionleft->border->region->f;  // strata with fold lines counting twice
      //d = a->depths[0];
      parity = 0;
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
cellcomplex_printnodes (struct ccomplex *cc, int lverbose, int *noderemap)
{
  int n;
  struct ccomplexnode *node;

  for (n = 0; n < cc->nodedim; n++)
  {
    node = cc->nodes + n;
    if (node->type == CC_REMOVED) continue;
    printf ("node %d", noderemap[n]);
    if (lverbose) printf (" ref %d", node->refcount);
    if (lverbose >= 2)
    {
      printf (" surfcc %d", node->surfcc);
      if (node->type == CC_NODETYPE_CUSP) 
        printf (" of CUSP type, ne-arc %d, cusp %d", node->ne->tag, node->cusp);
       else
        printf (" of type %d, ne-arc %d, stratum %d", node->type, node->ne->tag, node->stratum);
    }
    printf ("\n");
  }
}

void
cellcomplex_printarcs (struct ccomplex *cc, int lverbose, int *noderemap, int *arcremap)
{
  int n;
  struct ccomplexarc *arc;

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type == CC_REMOVED) continue;
    printf ("arc %d [%d %d]", arcremap[n], noderemap[arc->enda], noderemap[arc->endb]);
    if (lverbose) printf (" ref %d", arc->refcount);
    if (lverbose >= 2)
    {
      printf (", stratum %d, type %d", arc->stratum, arc->type);
      if (arc->type == CC_ARCTYPE_CUT || arc->type == CC_ARCTYPE_FOLD)
        printf (" (%d-%d), tag %d", arc->cusp1, arc->cusp2, arc->arc->tag);
    }
    printf ("\n");
  }
}

void
cellcomplex_printfaces (struct ccomplex *cc, int lverbose, int *arcremap)
{
  int n, i, na, nr;
  int *ivec;
  struct ccomplexface *face;

  for (n = 0, nr = 0; n < cc->facedim; n++)
  {
    face = cc->faces + n;
    if (face->type == CC_REMOVED) continue;
    ivec = face->faceborder;
    printf ("face %d [", nr++);
    for (i = 0; i < face->facebordernum; i++)
    {
      na = onarc2narc(ivec[i]);
      printf ("%c", (ivec[i]>0)?'+':'-');
      printf ("%d ", arcremap[na]);
    }
    printf ("]");
    if (lverbose >= 2) printf (", stratum %d", face->stratum);
    printf ("\n");
  }
}

void
cellcomplex_print (struct ccomplex *cc, int lverbose)
{
  int n, nr, *noderemap, *arcremap;
  struct ccomplexnode *node;
  struct ccomplexarc *arc;

  noderemap = (int *) malloc (cc->nodedim*sizeof(int));
  arcremap = (int *) malloc (cc->arcdim*sizeof(int));

  for (n = 0, nr = 0; n < cc->nodedim; n++)
  {
    noderemap[n] = -1;
    node = cc->nodes + n;
    if (node->type != CC_REMOVED) noderemap[n] = nr++;
  }
  for (n = 0, nr = 0; n < cc->arcdim; n++)
  {
    arcremap[n] = -1;
    arc = cc->arcs + n;
    if (arc->type != CC_REMOVED) arcremap[n] = nr++;
  }
  cellcomplex_printnodes (cc, lverbose, noderemap);
  cellcomplex_printarcs (cc, lverbose, noderemap, arcremap);
  cellcomplex_printfaces (cc, lverbose, arcremap);
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

/*
 * compute the Euler characteristic of a cell complex
 */

int
complex_characteristic (struct ccomplex *cc)
{
  if (debug) cellcomplex_checkconsistency (cc);
  if (debug) printf ("Nodes: %d, arcs: %d, faces: %d\n", cc->nodenum, cc->arcnum, cc->facenum);

  return (cc->nodenum - cc->arcnum + cc->facenum);
}

