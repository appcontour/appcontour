/*
 * computation of the fundamental group of the interior of
 * the surface
 */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "contour.h"
#include "fundamental.h"

#define MAXGENERATORS 26

extern int debug;
extern int quiet;
extern int verbose;
extern int interactive;

void
compute_fundamental (struct ccomplex *cc)
{
  int ccnum;
  //int count;
  struct ccomplexcc *cccc;

  //count = complex_melt (cc);
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
    if (debug) print_presentation (cccc->p);
    if (interactive >= 2) fg_interactive (cccc->p);
    simplify_presentation (cccc->p);
    if (interactive) fg_interactive (cccc->p);
    print_presentation (cccc->p);
    if (verbose) print_exponent_matrix (cccc->p);
  }
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
  int *ivec, *variable;

  p = (struct presentation *) malloc (sizeof (struct presentation));
  p->gennum = 0;
  p->rules = 0;

  for (n = 0; n < cc->arcdim; n++)
  {
    arc = cc->arcs + n;
    if (arc->type == CC_REMOVED) continue;
    if (arc->isinspanningtree) continue;
    node = cc->nodes + arc->enda;
    if (node->cc != cccc) continue;
    p->gennum++;
  }
  if (p->gennum == 0) return (p);

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
void nielsen_mulright (struct presentation *p,
                       struct presentationrule *r1,
                       struct presentationrule *r2, int expon);

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
      for (r = p->rules; r; r = r->next)
        r = sp_do_substitute (r, subvar, replaceword, optsublen);
      sp_do_eliminatevar (p, subvar);
      free (replaceword);
      count += simplify_presentation2 (p);
    }
  }

  return (count);
}

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
    // last operation: empty old rule
    r->length = 0;
    return (newr);
  }
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
            nielsen_mulright (p, lr, sr, -1);
          } else {
            nielsen_mulright (p, lr, sr, 1);
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
}

/*
 * print presentation of fundamental group
 */

void
print_presentation (struct presentation *p)
{
  struct presentationrule *r;
  int generator, rulenum, i, g;
  char var;

  if (p->gennum == 0)
  {
    if (quiet) printf ("<>\n");
     else printf ("Trivial group\n<>\n");
    return;
  }
  for (rulenum = 0, r = p->rules; r; r = r->next) rulenum++;
  if (!quiet)
  {
    if (rulenum == 0)
      printf ("Free group of rank %d\n", p->gennum);
     else
      printf ("Finitely presented group with %d generator%s\n", p->gennum,
                (p->gennum == 1)?"":"s");
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
  for (r = p->rules; r; r = r->next)
  {
    if (r != p->rules) printf (", ");
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
      if (generator >= p->gennum) var = '?';
      printf ("%c", var);
    }
  }
  printf (">\n");
}

/*
 * print exponent matrix
 */

void
print_exponent_matrix (struct presentation *p)
{
  struct presentationrule *r;
  int i, j, sum;

  printf ("Exponent sum matrix:\n");
  for (i = 1; i <= p->gennum; i++)
  {
    printf ("[");
    for (r = p->rules; r; r = r->next)
    {
      sum = 0;
      for (j = 0; j < r->length; j++)
      {
        if (abs(r->var[j]) != i) continue;
        if (r->var[j] > 0) sum++; else sum--;
      }
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
void nielsen_invrel (struct presentationrule *r);
void rotate_relator (struct presentationrule *r, int rots);
struct presentationrule *nielsen_exchange_relators (struct presentationrule *list,
                         struct presentationrule *r1,
                         struct presentationrule *r2);
void nielsen_invgen (struct presentation *p, int n);
void nielsen_exchange_generators (struct presentation *p, int m, int n);
void nielsen_xisabtok (struct presentation *p, int m, int n, int k);
struct presentationrule *nielsen_xisabtokl (struct presentationrule *r, int m, int n, int posk, int signk);

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
  int kmax, direction;

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
      printf ("Relators and simplification:\n");
      printf (" rotrel <n> <rots> : rotate relator <n> left <rots> times\n");
      printf (" simplify : perform a complete simplification (might change signs\n");
      printf (" test r1 r2: find common substring\n");
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
      nielsen_invrel (r);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "rotcol") == 0)
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
        nielsen_invrel (r2);
        kmax = sp_fcs_r1r2 (r1, r2, 0, 1);
        nielsen_invrel (r2);
        printf ("Found common substring of length %d\n", kmax);
        //print_presentation (p);
      }
      nielsen_mulright (p, r1, r2, 1);
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
      nielsen_mulright (p, r1, r2, -1);
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
      p->rules = nielsen_exchange_relators (p->rules, r1, r2);
      print = 1;
      continue;
    }
    if (strcasecmp (cmd, "negrow") == 0)
    {
      n = atoi (chpt);
      if (n <= 0 || n > p->gennum) {printf ("Invalid row number %d\n", n); continue;}
      nielsen_invgen (p, n);
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
      nielsen_exchange_generators (p, m, n);
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
      nielsen_xisabtok (p, n, m, -1);
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
      nielsen_xisabtok (p, n, m, 1);
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
nielsen_invrel (struct presentationrule *r)
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
nielsen_exchange_relators (struct presentationrule *list,
                           struct presentationrule *r1,
                           struct presentationrule *r2)
{
  struct presentationrule *r, *r1next, *r2next;
  assert (list != r2);
  assert (list && r1 && r2);

  if (list != r1)
  {
    list->next = nielsen_exchange_relators (list->next, r1, r2);
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
nielsen_mulright (struct presentation *p,
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
nielsen_invgen (struct presentation *p, int g)
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
}

void
nielsen_exchange_generators (struct presentation *p, int m, int n)
{
  struct presentationrule *r;
  int i;

  assert (m > 0 && m <= p->gennum);
  assert (n > 0 && n <= p->gennum);

  if (m == n) return;

  for (r = p->rules; r; r = r->next)
  {
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
}

void
nielsen_xisabtok (struct presentation *p, int m, int n, int k)
{
  int posk, signk;

  assert (m > 0 && m <= p->gennum);
  assert (n > 0 && n <= p->gennum);
  assert (m != n);

  if (k == 0) return;

  signk = 1;
  posk = k;
  if (k < 0) {signk = -1; posk = -k;}

  p->rules = nielsen_xisabtokl (p->rules, m, n, posk, signk);

  sp_simplifyword (p);

  return;
}

struct presentationrule *
nielsen_xisabtokl (struct presentationrule *r, int m, int n, int posk, int signk)
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
    r->next = nielsen_xisabtokl (r->next, m, n, posk, signk);
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
  assert (j == r->length + numvarm);

  newr->next = r->next;
  r->next = 0;
  free (r);
  newr->next = nielsen_xisabtokl (newr->next, m, n, posk, signk);
  return (newr);
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
 */

struct ccomplex *
compute_cellcomplex (struct sketch *s, int fg_type)
{
  extern int finfinity;
  struct ccomplex *cc;
  int euler, surfeuler, realeuler;
  int status;

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
  }
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
  if (debug) printf ("Euler characteristic: %d = %d nodes - %d arcs + %d faces.\n",
             euler, cc->nodenum, cc->arcnum, cc->facenum);
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
    cccc->p = 0;
    cccc->next = 0;
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

