/*
 * computation of Alexander polinomial of
 * a finitely presented group with infinite ciclic
 * commutator quotient
 */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "contour.h"
#include "fundamental.h"
#include "alexander.h"
#include "parser.h"

int
alexander (struct presentation *p)
{
  struct presentationrule *r;
  int i, sum, rank, matrixrank;

  topreabelian (p);

  rank = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    assert (sum >= 0);
    if (sum) matrixrank = i;
    if (sum && sum != 1)
    {
      printf ("Cannot compute Alexander polynomial for groups with torsion\n");
      return (0);
    }
  }
  rank = p->gennum - matrixrank;
  if (rank != 1)
  {
    printf ("Cannot compute Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }

  printf ("Not implemented\n");
  return (0);
}
