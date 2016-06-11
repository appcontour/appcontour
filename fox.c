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
#include <laurent.h>
#include <fox.h>

void
foxjacobian (struct presentation *pr)
{
  extern int verbose, abelianize, preabelian;
  struct presentationrule *r;
  struct laurentpoly *p;
  int rank = 1, i, j, e;

  if (abelianize == 1 && preabelian == 0)
  {
    printf ("In order to abelianize we need a preabelian presentation, forcing its computation...\n");
    printf ("You can avoid this message by including the --preabelian option\n");
    preabelian = 1;
  }
  if (preabelian) topreabelian (pr);
  if (verbose) print_presentation (pr);

  if (preabelian) rank = compute_fg_rank (pr);
  if (abelianize == 2) rank = 0;

  assert (rank >= 0);
  for (i = 0, r = pr->rules; r; i++, r = r->next)
  {
    for (j = 0; j < pr->gennum; j++)
    {
      //printf ("(%d,%d): ", i+1, j+1);
      if (j > 0) printf (",\t");
      if (rank == 0)
      {
        e = get_exp_sum (r, j+1);
        printf ("%d", e);
        continue;
      }
      p = foxderivative (r, j, pr->gennum - rank, rank);
      if (abelianize == 0) continue;
      if (rank == 1)
      {
        print_laurentpoly (p, "t");
      } else {
        print_laurentpoly (p, "uvwxyz");
      }
      //printf ("\n");
    }
    printf (";\n");
  }
}

/*
 * compute the fox derivative of word r with respect to generator gen (zero-based)
 * then project the result onto the last "rank" generators (presumably the abelianizer
 * map)
 */

#define MAX_RANK_ALLOWED 20
static int exponvec[MAX_RANK_ALLOWED];

struct laurentpoly *
foxderivative (struct presentationrule *r, int gen, int offset, int rank)
{
  extern int abelianize;
  struct laurentpoly *p = 0;
  int j, k;

  assert (rank <= MAX_RANK_ALLOWED);

  for (j = 0; j < r->length; j++)
  {
    if (abs(r->var[j]) != gen+1) continue;
    if (r->var[j] > 0)
    {
      if (abelianize == 0)
      {
        printf ("+");
        if (j == 0) printf ("1");
        for (k = 0; k < j; k++)
        {
          if (r->var[k]>0) printf ("%c", r->var[k]-1+'a');
          else printf ("%c", -r->var[k]-1+'A');
        }
      } else {
        count_and_map (&r->var[0], j, offset, rank, exponvec);
        p = laurentpoly_addmonom (p, rank, exponvec, 1);
      }
    } else {
      if (abelianize == 0)
      {
        printf ("-");
        for (k = 0; k <= j; k++)
        {
          if (r->var[k]>0) printf ("%c", r->var[k]-1+'a');
          else printf ("%c", -r->var[k]-1+'A');
        }
      } else {
        count_and_map (&r->var[0], j + 1, offset, rank, exponvec);
        p = laurentpoly_addmonom (p, rank, exponvec, -1);
      }
    }
  }
  //if (abelianize == 0) printf ("\n");
  return (p);
}

/*
 * count_and_map scans an integer vector and sums the exponents correctly
 */

void
count_and_map (int *vec, int len, int offset, int indets, int *exponvec)
{
  int i, sign;

  for (i = 0; i < indets; i++) exponvec[i] = 0;
  for (i = 0; i < len; i++)
  {
    if (abs(vec[i]) - 1 < offset) continue;
    if (abs(vec[i]) - 1 >= offset + indets) continue;
    sign = 1;
    if (vec[i] < 0) sign = -1;
    exponvec[abs(vec[i]) - 1 - offset] += sign;
  }
}
