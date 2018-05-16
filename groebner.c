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
  int j, ruledegree;

  *statuspt = 0;
  assert (rule && rule->degree >= 0);
  if (p == 0) return (0);
  ruledegree = rule->degree;

  for (j = p->degree; j >= ruledegree; j--)
  {
    printf ("RULE higher: NOT IMPLEMENTED\n");
  }

  assert ( *statuspt == 0 );  // TODO reminder, here we MUST canonify stem

  for (j = 0; p->degree >= ruledegree + j; j++)
  {
    printf ("RULE lower: NOT IMPLEMENTED\n");
  }

  assert ( *statuspt == 0 );  // TODO reminder, here we MUST canonify stem
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
  for (i = 0; i <= stem->degree; i++) lp->stem[i].l0 = stem->coef[i];

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
  printf ("printout_si non implemented\n");
}
