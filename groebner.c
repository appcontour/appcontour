/*
 * implementing algorithm for Groebner basis computation
 * based on ideas in KanKap:88 adapted for the Laurent case
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "laurent.h"
#include "alexander.h"
#include "groebner.h"

struct alexanderideal *
groebner1 (struct alexanderideal *ai)
{
  int i, reductions, newgenerators;
  extern int experimental;
  struct stemideal *si;

  if (ai->indets < 1) return (ai);
  assert (ai->indets == 1);

  if (experimental)
  {
    printf ("l1num = %d, max_generators_num = %d\n", ai->l1num, ai->max_generators_num);
    for (i = 0; i < ai->l1num; i++)
    {
      print_laurentpoly (ai->l[i], "t");
      printf ("\n");
    }
  }

  si = ai2si (ai);
  while (1)
  {
    reductions = groebner1_tryreduce (si);
if (experimental)
{
  printf ("Performed %d reductions\n", reductions);
  printout_si (si);
}
    assert (reductions >= 0);
    if (experimental || reductions > 0) groebner1_dropzeros (si);
if (experimental)
{
  printf ("after zero polynomials drop\n");
  printout_si (si);
}
    newgenerators = groebner1_add_spolynomials (si);
    assert (newgenerators >= 0);
    if (reductions + newgenerators <= 0) break;
  }

  ai = si2ai (si, ai);
  free_stemideal (si);
  return (ai);
}

/*
 * remove polynomials that got reduced to zero
 */

int
groebner1_dropzeros (struct stemideal *si)
{
  int i, j;
  extern int verbose, experimental;

  for (i = 0, j = 0; i < si->num; i++)
  {
    if (si->stem[i] != 0)
    {
      si->stem[j++] = si->stem[i];
    } else {
      if (verbose) printf ("Ideal simplification: zero polynomial removed\n");
    }
  }
  if (experimental) printf ("i - j = %d\n", i - j);
  si->num = j;
  return (i - j);
}

/*
 * reduction process
 */

int
groebner1_tryreduce (struct stemideal *si)
{
  int reduce, totreduce;

  totreduce = 0;
  do
  {
    reduce = groebner1_tryreduce_once (si);
    totreduce += reduce;
  } while (reduce > 0);

  return (totreduce);
}

/*
 * reduction process cycle
 */

int
groebner1_tryreduce_once (struct stemideal *si)
{
  int i, j, retcode;
  int reductions = 0;

  for (i = 0; i < si->num; i++)
  {
    if (si->stem[i] == 0) continue;
    assert (si->stem[i]->degree >= 0);
    /* use i-th polynomial to try and reduce all others */
    for (j = 0; j < si->num; j++)
    {
      if (i == j) continue;
      si->stem[j] = groebner1_reduce_using_rule (si->stem[j], si->stem[i], &retcode);
      if (retcode > 0) reductions += retcode;
    }
  }
  return (reductions);
}

/*
 * reduce a laurent polynomial using another laurent polynomial
 */

struct stem *
groebner1_reduce_using_rule (struct stem *p, struct stem *rule, int *statuspt)
{
  int bshift, tshift, j, j1, j2, ruledegree, pdegree;
  Stemint tcoeff, bcoeff, quotient;
  int bzeros, tzeros, dummy;
  Stemint maxcsize;

  *statuspt = 0;
  assert (rule && rule->degree >= 0);
  if (p == 0) return (0);
printf ("reducing ");
printout_stem (p);
printf (" using rule: ");
printout_stem (rule);
printf ("\n");
  ruledegree = rule->degree;
  pdegree = p->degree;

  tcoeff = rule->coef[rule->degree];
  bcoeff = rule->coef[0];
  maxcsize = stem_linf (rule);
  for (bshift = 0, tshift = pdegree - ruledegree; tshift >= bshift; bshift++, tshift--)
  {
    quotient = gb_int_div (p->coef[pdegree - bshift], tcoeff);
    if (quotient)
    {
printf ("TOP, quotient = %lld: ", quotient);
      /* TODO: here we need a check on integer overflow */
      if (ll_safety_check (quotient, maxcsize) == 0) return (p);
      (*statuspt)++;
      for (j1 = 0, j2 = pdegree - bshift - ruledegree; j1 <= ruledegree; j1++, j2++)
      {
        p->coef[j2] -= quotient*rule->coef[j1];
      }
printout_stem (p);
printf ("\n");
      if (stem_linf (p) >= LLONG_MAX/2)
      {
        /* backtrack if coefficients grow too much */
        for (j1 = 0, j2 = pdegree - bshift - ruledegree; j1 <= ruledegree; j1++, j2++)
        {
          p->coef[j2] += quotient*rule->coef[j1];
        }
        return (p);
      }
    }
    if (tshift > bshift)
    {
      quotient = gb_int_div (p->coef[bshift], bcoeff);
      if (quotient)
      {
printf ("BOTTOM, quotient = %lld: ", quotient);
        if (ll_safety_check (quotient, maxcsize) == 0) return (p);
        (*statuspt)++;
        for (j1 = 0, j2 = bshift; j1 <= ruledegree; j1++, j2++)
        {
          p->coef[j2] -= quotient*rule->coef[j1];
        }
printout_stem (p);
printf ("\n");
        if (stem_linf (p) >= LLONG_MAX/2)
        {
          /* backtrack if coefficients grow too much */
          for (j1 = 0, j2 = bshift; j1 <= ruledegree; j1++, j2++)
          {
            p->coef[j2] += quotient*rule->coef[j1];
          }
          return (p);
        }
      }
    }
  }

  bzeros = 0;
  for (j = 0; j <= pdegree; j++)
  {
    if (p->coef[j] == 0) bzeros++;
      else break;
  }
  if (bzeros >= pdegree + 1)
  {
    free (p);
    return (0);
  }
  tzeros = 0;
  for (j = pdegree; j >= 0; j--)
  {
    if (p->coef[j] == 0) tzeros++;
      else break;
  }

  if (bzeros > 0)
  {
    for (j1 = 0, j2 = bzeros; j2 <= pdegree; j1++, j2++) p->coef[j1] = p->coef[j2];
  }
  pdegree -= bzeros + tzeros;
  p->degree = pdegree;

  assert (p->coef[pdegree]);
  if (p->coef[pdegree] < 0)
  {
printf ("CHANGE SIGN of "); printout_stem (p); printf ("\n");
    /* normalize sign of leading term */
    for (j = 0; j <= pdegree; j++)
    {
      p->coef[j] = -p->coef[j];
    }
    /* recursively reduce... the leading term will not change and
     * we don't risk infinite recursion
     */
    p = groebner1_reduce_using_rule (p, rule, &dummy);
printf ("  ---> "); printout_stem (p); printf ("\n");
  }
  return p;
}

/*
 * compute and add S-polynomials
 */

int
groebner1_add_spolynomials (struct stemideal *si)
{
  printf ("ADD-S-POLYNOMIALS NOT IMPLEMENTED\n");
  return (0);
}

/*
 * convert from alexander ideal to the simpler stemideal
 */

struct stemideal *
ai2si (struct alexanderideal *ai)
{
  struct stemideal *si;
  int i;
  int numalloc = ai->l1num + GB_EXTRAROOM(ai->l1num);

  si = (struct stemideal *) malloc (STEMIDEALSIZE(numalloc));
  si->dim = numalloc;
  si->num = ai->l1num;
  for (i = 0; i < ai->l1num; i++) si->stem[i] = lp2stem (ai->l[i]);

  return (si);
}

struct stem *
lp2stem (struct laurentpoly *lp)
{
  struct stem *stem;
  int i;

  stem = (struct stem *) malloc (STEMSIZE(lp->stemdegree + 1));
  stem->dim = lp->stemdegree + 1;
  stem->degree = lp->stemdegree;

  for (i = 0; i <= lp->stemdegree; i++) stem->coef[i] = lp->stem[i].l0;

  return (stem);
}

struct alexanderideal *
si2ai (struct stemideal *si, struct alexanderideal *ai)
{
  int i;

  assert (ai->max_generators_num >= si->num);

  for (i = 0; i < ai->l1num; i++) {free_laurentpoly (ai->l[i]); ai->l[i] = 0;}
  for (i = 0; i < si->num; i++)
  {
    ai->l[i] = stem2lp (si->stem[i]);
  }
  ai->l1num = si->num;
  return (ai);
}

struct laurentpoly *
stem2lp (struct stem *stem)
{
  int i;
  struct laurentpoly *lp;

  lp = (struct laurentpoly *) malloc (POLYSIZE (stem->degree + 1));
  lp->stemdegree = stem->degree;
  lp->minexpon = 0;
  lp->indets = 1;
  for (i = 0; i <= stem->degree; i++)
  {
    assert (INT_MIN <= stem->coef[i] && stem->coef[i] <= INT_MAX);
    lp->stem[i].l0 = stem->coef[i];
  }

  return (lp);
}

void
free_stemideal (struct stemideal *si)
{
  int i;

  for (i = 0; i < si->num; i++)
  {
    assert (si->stem[i]);
    free (si->stem[i]);
  }

  free (si);
}

void
printout_si (struct stemideal *si)
{
  int i;

  printf ("Stem ideal generated by %d stems:\n", si->num);
  for (i = 0; i < si->num; i++)
  {
    printout_stem (si->stem[i]);
    printf ("\n");
  }
}

void
printout_stem (struct stem *stem)
{
  int i;

  if (stem == 0) {printf ("0"); return;}

  for (i = 0; i <= stem->degree; i++)
  {
    printf ("%+llds^%d", stem->coef[i], i);
  }
}

Stemint
stem_linf (struct stem *stem)
{
  Stemint maxc = 0;
  int i;

  if (stem == 0) return (0);
  for (i = 0; i <= stem->degree; i++)
  {
    if (llabs(stem->coef[i]) > maxc) maxc = llabs(stem->coef[i]);
  }
  return maxc;
}

Stemint
gb_int_div (Stemint dividend, Stemint divisor)
{
  Stemint rest, quotient;

  assert (divisor);
  if (divisor < 0) {divisor = -divisor; dividend = -dividend;}
  quotient = dividend/divisor;

  rest = dividend - quotient*divisor;
  if (rest < 0)
  {
    rest += divisor;
    quotient--;
  }

  if (2*rest > divisor) quotient++;

  return (quotient);
}
