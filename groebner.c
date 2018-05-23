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
  extern int verbose, factorideal;
  struct stemideal *si;
  struct stem *gcd = 0;

  if (ai->indets < 1) return (ai);
  assert (ai->indets == 1);

  if (verbose > 2)
  {
    printf ("Computing groebner base for Laurent in one indet,\n");
    printf ("ideal generated by %d polynomials (max = %d):\n", ai->l1num, ai->max_generators_num);
    for (i = 0; i < ai->l1num; i++)
    {
      print_laurentpoly (ai->l[i], "t");
      printf ("\n");
    }
  }

  si = ai2si (ai);
  groebner1_dropzeros (si);
  while (1)
  {
    reductions = groebner1_tryreduce (si);
    assert (reductions >= 0);
    if (reductions > 0) groebner1_dropzeros (si);
    if (verbose > 2)
    {
      printf ("after %d reductions:\n", reductions);
      printout_si (si);
      printf ("\n");
    }
    newgenerators = groebner1_add_spolynomials (si);
    assert (newgenerators >= 0);
    if (newgenerators > 0) groebner1_dropzeros (si);
    if (verbose > 2)
    {
      printf ("Added %d new S-polynomials to obtain:\n", newgenerators);
      printout_si (si);
      printf ("\n");
    }
    if (reductions + newgenerators <= 0) break;
  }

  if (factorideal)
  {
    if (verbose) printf ("Computing gcd\n");
    gcd = stemideal_gcd (si);
  }
  ai = si2ai (si, ai);
  if (gcd)
  {
    ai->gcd = stem2lp (gcd);
    free (gcd);
  }
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
  extern int verbose;

  for (i = 0, j = 0; i < si->num; i++)
  {
    if (si->stem[i] != 0)
    {
      si->stem[j++] = si->stem[i];
    } else {
      if (verbose) printf ("Ideal simplification: zero polynomial removed\n");
    }
  }
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
//printf ("Trying to reduce [%d,%d] - ", i, j); printout_stem (si->stem[j]);
//printf (" using "); printout_stem (si->stem[i]); printf ("\n");
      si->stem[j] = groebner1_reduce_using_rule (si->stem[j], si->stem[i], &retcode);
//printf (" -->[retcode = %d] ", retcode); printout_stem (si->stem[j]); printf ("\n");
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
  int bzeros, tzeros, previous_are_zero;
  //int dummy;
  Stemint maxcsize;

  *statuspt = 0;
  assert (rule && rule->degree >= 0);
  if (p == 0) return (0);
  ruledegree = rule->degree;
  pdegree = p->degree;

  tcoeff = rule->coef[rule->degree];
  bcoeff = rule->coef[0];
  maxcsize = stem_linf (rule);
  previous_are_zero = 1;
  for (bshift = 0, tshift = pdegree - ruledegree; tshift >= bshift; bshift++, tshift--)
  {
    quotient = gb_int_div (p->coef[pdegree - bshift], tcoeff);
    if (quotient)
    {
      if (ll_safety_check (quotient, maxcsize) == 0) return (p);
      (*statuspt)++;
      for (j1 = 0, j2 = pdegree - bshift - ruledegree; j1 <= ruledegree; j1++, j2++)
      {
        p->coef[j2] -= quotient*rule->coef[j1];
      }
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
    if (previous_are_zero  && p->coef[pdegree - bshift] < 0)
    {
      /* normalize sign of leading term */
      for (j = 0; j <= pdegree; j++)
      {
        p->coef[j] = -p->coef[j];
      }
    }
    if (p->coef[pdegree - bshift]) previous_are_zero = 0;
    if (tshift > bshift)
    {
      quotient = gb_int_div (p->coef[bshift], bcoeff);
      if (quotient)
      {
        if (ll_safety_check (quotient, maxcsize) == 0) return (p);
        (*statuspt)++;
        for (j1 = 0, j2 = bshift; j1 <= ruledegree; j1++, j2++)
        {
          p->coef[j2] -= quotient*rule->coef[j1];
        }
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
    /* normalize sign of leading term */
    for (j = 0; j <= pdegree; j++)
    {
      p->coef[j] = -p->coef[j];
    }
  }
  return p;
}

/*
 * compute and add S-polynomials
 */

int
groebner1_add_spolynomials (struct stemideal *si)
{
  struct stem *spoltop, *spolbottom;
  int i, j, startnum;
  int added = 0;

  startnum = si->num;
  assert (si->dim > startnum);  /* make sure we have space for at least one new polynomial */
  for (i = 0; i < startnum - 1; i++)
  {
    for (j = i+1; j < startnum; j++)
    {
      build_S_pols (si->stem[i], si->stem[j], &spoltop, &spolbottom);
      if (spoltop) spoltop = reduce_pol_si (spoltop, si);
      if (spoltop)
      {
        added++;
        assert (si->num < si->dim);
        si->stem[si->num++] = spoltop;
        if (si->num >= si->dim)
        {
          if (spolbottom) free (spolbottom);   /* avoid memory leak */
          return (added);
        }
      }
      if (spolbottom) spolbottom = reduce_pol_si (spolbottom, si);
      if (spolbottom)
      {
        added++;
        assert (si->num < si->dim);
        si->stem[si->num++] = spolbottom;
        if (si->num >= si->dim) return (added);
      }
    }
  }
  return (added);
}

/*
 *
 */

int
build_S_pols (struct stem *p1, struct stem *p2, struct stem **spoltoppt, struct stem **spolbottompt)
{
  int i, j, p1deg, p2deg;
  Stemint a, b, maxcsize;
  struct stem *spoltop, *spolbottom;

  if (p2->degree > p1->degree) return (build_S_pols (p2, p1, spoltoppt, spolbottompt));

  p1deg = p1->degree;
  p2deg = p2->degree;
  spoltop = spolbottom = *spoltoppt = *spolbottompt = 0;
  if (p1deg == p2deg) return (0);  /* this should lead to reduction in the reduction stage */
  assert (p1deg > p2deg);

  /* now p1 is the polynomial with larger degree */
  maxcsize = stem_linf (p1);
  /* compute the S-pol from the top */
  a = gb_int_div (p2->coef[p2deg], p1->coef[p1deg]);
  if (a && ll_safety_check (a, maxcsize))
  {
    b = p2->coef[p2deg] - a*p1->coef[p1deg];

    spoltop = (struct stem *) malloc (STEMSIZE(p1deg + 1));
    spoltop->dim = p1deg + 1;
    spoltop->degree = p1deg;
    spoltop->coef[p1deg] = b;
    for (i = 0; i < p1deg; i++) spoltop->coef[i] = - a * p1->coef[i];
    for (i = 0, j = p1deg-p2deg; i < p2deg; i++, j++) spoltop->coef[j] += p2->coef[i];
    spoltop = stem_normalize (spoltop);
    if (stem_linf (spoltop) >= LLONG_MAX/2)
    {
      free (spoltop);
      signal_int_overflow();
      spoltop = 0;
    }
  }
  /* now compute the S-pol from the bottom */
  a = gb_int_div (p2->coef[0], p1->coef[0]);
  if (a && ll_safety_check (a, maxcsize))
  {
    b = p2->coef[0] - a*p1->coef[0];

    spolbottom = (struct stem *) malloc (STEMSIZE(p1deg + 1));
    spolbottom->dim = p1deg + 1;
    spolbottom->degree = p1deg;
    spolbottom->coef[0] = b;
    for (i = 1; i <= p1deg; i++) spolbottom->coef[i] = - a * p1->coef[i];
    for (i = 1; i <= p2deg; i++) spolbottom->coef[i] += p2->coef[i];
    spolbottom = stem_normalize (spolbottom);
    if (stem_linf (spolbottom) >= LLONG_MAX/2)
    {
      free (spolbottom);
      signal_int_overflow();
      spolbottom = 0;
    }
  }

  *spoltoppt = spoltop;
  *spolbottompt = spolbottom;
  if (spoltop && spolbottom) return (2);
  if (spoltop || spolbottom) return (1);
  return (0);
}

/*
 *
 */

struct stem *
reduce_pol_si (struct stem *spol, struct stemideal *si)
{
  int status;

  while (1)
  {
    status = 0;
    spol = reduce_pol_si_cycle (spol, si, &status);
    if (status == 0) break;
  }
  return spol;
}

/*
 *
 */

struct stem *
reduce_pol_si_cycle (struct stem *spol, struct stemideal *si, int *statuspt)
{
  int i;

  for (i = 0; i < si->num; i++)
  {
    spol = groebner1_reduce_using_rule (spol, si->stem[i], statuspt);
  }

  return spol;
}

/*
 * stem_normalize: normalize with respect to a unit
 */

struct stem *
stem_normalize (struct stem *stem)
{
  int j, j1, j2, bzeros, tzeros, degree;

  bzeros = 0;
  degree = stem->degree;
  for (j = 0; j <= degree; j++)
  {
    if (stem->coef[j] == 0) bzeros++;
      else break;
  }
  if (bzeros >= degree + 1)
  {
    free (stem);
    return (0);
  }

  tzeros = 0;
  for (j = degree; j >= 0; j--)
  {
    if (stem->coef[j] == 0) tzeros++;
      else break;
  }

  if (bzeros > 0)
  {
    for (j1 = 0, j2 = bzeros; j2 <= degree; j1++, j2++) stem->coef[j1] = stem->coef[j2];
  }

  degree -= bzeros + tzeros;
  stem->degree = degree;

  assert (stem->coef[degree]);
  if (stem->coef[degree] < 0)
  {
    /* normalize sign of leading term */
    for (j = 0; j <= degree; j++)
    {
      stem->coef[j] = -stem->coef[j];
    }
  }
  return (stem);
}

/*
 * compute gcd of the set in a stemideal
 */

struct stem *
stemideal_gcd (struct stemideal *si)
{
  int i;
  struct stem *newgcd, *oldgcd;

  assert (si->num > 1);

  newgcd = stem_euclid (si->stem[0], si->stem[1]);
  for (i = 2; i < si->num; i++)
  {
    oldgcd = newgcd;
    newgcd = stem_euclid (oldgcd, si->stem[i]);
    free (oldgcd);
  }
  return (newgcd);
}

/*
 * compute the gcd between two stems
 */

struct stem *
stem_euclid (struct stem *p1, struct stem *p2)
{
  struct stem *gcd;

/* XXX copiare qui laurentpoly_euclid come idea di lavoro */

  printf ("stem_gcd not yet implemented, returning 1\n");
  gcd = (struct stem *) malloc (STEMSIZE(1));
  gcd->dim = 1;
  gcd->degree = 0;
  gcd->coef[0] = 1;

  return (gcd);
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
  int i, sign;

  if (lp == 0) return (0);
  stem = (struct stem *) malloc (STEMSIZE(lp->stemdegree + 1));
  stem->dim = lp->stemdegree + 1;
  stem->degree = lp->stemdegree;

  sign = 1;
  if (lp->stem[lp->stemdegree].l0 < 0) sign = -1;
  for (i = 0; i <= lp->stemdegree; i++) stem->coef[i] = sign*lp->stem[i].l0;

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
