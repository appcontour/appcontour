/*
 * computation of the jacobian of the presentation
 * (transformed in a preabelian form)
 * using Fox calculus.
 * The result is then mapped through the abelianizer map
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <fundamental.h>
#include <fox.h>

void
foxjacobian (struct presentation *p)
{
  extern int verbose;

  topreabelian (p);
  if (verbose) print_presentation (p);
  printf ("Not yet implemented\n");
}
