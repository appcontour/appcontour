/*
 * implementing algorithm for Groebner basis computation
 * based on ideas in KanKap:88 adapted for the Laurent case
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "laurent.h"
#include "alexander.h"
#include "groebner.h"

struct alexanderideal *
groebner1 (struct alexanderideal *ai)
{
  int i, reductions, newgenerators;
  extern int experimental;

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
  while (1)
  {
    reductions = groebner1_tryreduce (ai->l, ai->l1num);
if (experimental)
{
  printf ("Performed %d reductions\n", reductions);
  printout_ideal1 (ai, 0);
}
    assert (reductions >= 0);
    if (experimental || reductions > 0) ai->l1num -= groebner1_dropzeros (ai->l, ai->l1num);
if (experimental)
{
  printf ("after zero polynomials drop\n");
  printout_ideal1 (ai, 0);
}
    newgenerators = groebner1_add_spolynomials (ai->l, ai->l1num, ai->max_generators_num);
    assert (newgenerators >= 0);
    if (newgenerators > 0) ai->l1num += newgenerators;
    if (reductions + newgenerators <= 0) break;
  }

  return (ai);
}

/*
 * remove polynomials that got reduced to zero
 */

int
groebner1_dropzeros (struct laurentpoly *l[], int lnum)
{
  int i, j;
  extern int verbose, experimental;

  for (i = 0, j = 0; i < lnum; i++)
  {
    if (l[i] != 0)
    {
      l[j++] = l[i];
    } else {
      if (verbose) printf ("Ideal simplification: zero polynomial removed\n");
    }
  }
  if (experimental) printf ("i - j = %d\n", i - j);
  return (i - j);
}

/*
 * reduction process
 */

int
groebner1_tryreduce (struct laurentpoly *l[], int lnum)
{
  int reduce, totreduce;

  totreduce = 0;
  do
  {
    reduce = groebner1_tryreduce_once (l, lnum);
    totreduce += reduce;
  } while (reduce > 0);

  return (totreduce);
}

/*
 * reduction process cycle
 */

int
groebner1_tryreduce_once (struct laurentpoly *l[], int lnum)
{
  int i, j, retcode;
  int reductions = 0;

  for (i = 0; i < lnum; i++)
  {
    if (l[i] == 0) continue;
    assert (l[i]->stemdegree >= 0);
    /* use i-th polynomial to try and reduce all others */
    for (j = 0; j < lnum; j++)
    {
      if (i == j) continue;
      l[j] = groebner1_reduce_using_rule (l[j], l[i], &retcode);
      if (retcode > 0) reductions += retcode;
    }
  }
  return (reductions);
}

/*
 * reduce a laurent polynomial using another laurent polynomial
 */

struct laurentpoly *
groebner1_reduce_using_rule (struct laurentpoly *p, struct laurentpoly *rule, int *statuspt)
{
  int j, ruledegree;

  *statuspt = 0;
  assert (rule && rule->stemdegree >= 0);
  if (p == 0) return (0);
  ruledegree = rule->stemdegree;

  for (j = p->stemdegree; j >= ruledegree; j--)
  {
    printf ("RULE higher: NOT IMPLEMENTED\n");
  }

  assert ( *statuspt == 0 );  // TODO reminder, here we MUST canonify stem

  for (j = 0; p->stemdegree >= ruledegree + j; j++)
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
groebner1_add_spolynomials (struct laurentpoly *l[], int lnum, int ldim)
{
  printf ("ADD-S-POLYNOMIALS NOT IMPLEMENTED\n");
  return (0);
}
