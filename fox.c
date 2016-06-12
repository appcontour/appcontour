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
  extern int verbose, abelianize, preabelian, interactive;
  struct presentationrule *r;
  struct laurentpoly *p;
  int rank = 1, i, j, e, numrelators = 0;
  int maxlen, rid, k, *lincomb1, *lincomb2, gen;

  if (abelianize == 1 && preabelian == 0)
  {
    printf ("In order to abelianize we need a preabelian presentation, forcing its computation...\n");
    printf ("You can avoid this message by including the --preabelian option\n");
    preabelian = 1;
  }
  if (preabelian) {topreabelian (pr); rank = compute_fg_rank (pr);}
  if (abelianize == 2) rank = 0;
  if (verbose)
  {
    print_presentation (pr);
    if (preabelian) printf ("rank = %d\n", rank);
  }

  assert (rank >= 0);

  maxlen = 0;
  for (i = 0, r = pr->rules; r; i++, r = r->next) if (r->length > maxlen) maxlen = r->length;

  lincomb1 = (int *) malloc ((maxlen+1)*sizeof (int));
  lincomb2 = (int *) malloc ((maxlen+1)*sizeof (int));

  for (i = 0, r = pr->rules; r; i++, r = r->next)
  {
    numrelators++;
    for (j = 1; j <= pr->gennum; j++)
    {
      if (j > 1) printf (",\t");
      for (k = 0; k < r->length; k++) lincomb1[k] = 0;
      lincomb1[r->length] = 1;
      foxderivative (r->length, r->var, lincomb1, lincomb2, j);
      if (rank == 0)
      {
        e = map_to_trivial (r->var, lincomb2, r->length);
        printf ("%d", e);
        continue;
      }
      if (abelianize == 0) print_groupring_el (r->var, lincomb2, r->length);
      if (abelianize == 0) continue;
      p = map_to_abelian (r->var, lincomb2, r->length, pr->gennum - rank, rank);
      if (rank == 1)
      {
        print_laurentpoly (p, "t");
      } else {
        print_laurentpoly (p, "uvwxyz");
      }
    }
    printf (";\n");
  }

  if (interactive)
  {
    while (1)
    {
      printf ("Interactive mode: \n");
      print_presentation (pr);
      if (preabelian) printf ("rank = %d\n", rank);
      printf ("enter relator (1 to %d): ", numrelators);
      scanf ("%d", &rid);
      if (rid <= 0) break;
      if (rid > numrelators) {printf ("Invalid relator\n"); continue;}

      for (i = 1, r = pr->rules; r && i < rid; i++, r = r->next);
      for (k = 0; k <= r->length; k++) lincomb1[k] = 0;
      lincomb1[r->length] = 1;
      while (1)
      {
        print_groupring_el (r->var, lincomb1, r->length);
        printf ("\ntrivialized into: %d\n", map_to_trivial (r->var, lincomb1, r->length));
        if (preabelian)
        {
          p = map_to_abelian (r->var, lincomb1, r->length, pr->gennum - rank, rank);
          printf ("abelianized into: ");
          print_laurentpoly (p, "uvwxyz");
          printf ("\n");
        }
        printf ("enter generator num (0 to end): ");
        scanf ("%d", &gen);
        if (gen <= 0) break;
        printf ("Differentiating with respect to generator %d:\n", gen);
	foxderivative (r->length, r->var, lincomb1, lincomb2, gen);
        for (k = 0; k <= r->length; k++) lincomb1[k] = lincomb2[k];
      }
    }
  }
  free (lincomb1);
  free (lincomb2);
}

/*
 * map an element of the group-ring onto the last "rank" generators
 * (presumably the abelianizer map)
 */

#define MAX_RANK_ALLOWED 20
static int exponvec[MAX_RANK_ALLOWED];

struct laurentpoly *
map_to_abelian (int *rvar, int *lincomb, int len, int offset, int rank)
{
  struct laurentpoly *p = 0;
  int j;

  assert (rank <= MAX_RANK_ALLOWED);

  for (j = 0; j <= len; j++)
  {
    if (lincomb[j] == 0) continue;
    count_and_map (rvar, j, offset, rank, exponvec);
    p = laurentpoly_addmonom (p, rank, exponvec, lincomb[j]);
  }
  return (p);
}

/*
 * special case rank = 0: trivializer
 */

int
map_to_trivial (int *rvar, int *lincomb, int len)
{
  int j, val = 0;

  for (j = 0; j <= len; j++) val += lincomb[j];
  return (val);
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

void
print_groupring_el (int *rvar, int *lincomb, int len)
{
  int j, nothingprinted = 1;

  for (j = 0; j <= len; j++)
  {
    if (lincomb[j] == 0) continue;
    nothingprinted = 0;
    if (lincomb[j] > 0) printf ("+"); else printf ("-");
    if (abs(lincomb[j]) != 1 || j == 0) printf ("%d", abs(lincomb[j]));
    print_relator (rvar, j);
  }
  if (nothingprinted) printf ("+0");
}

void
print_relator (int *rvar, int len)
{
  int k;

  for (k = 0; k < len; k++)
  {
    if (rvar[k]>0) printf ("%c", rvar[k]-1+'a');
    else printf ("%c", -rvar[k]-1+'A');
  }
}

/*
 * compute the fox derivative of a (special) group-ring element  r with respect to generator gen
 * then project the result onto the last "rank" generators (presumably the abelianizer
 * map)
 */

void
foxderivative (int len, int *rvar, int *lincomb1, int *lincomb2, int gen)
{
  int k, kk;

  for (k = 0; k <= len; k++) lincomb2[k] = 0;

  for (k = 0; k <= len; k++)
  {
    /* differentiate the initial segment with respect to gen */
    if (lincomb1[k] == 0) continue;
    for (kk = 0; kk < k; kk++)
    {
      if (abs(rvar[kk]) != gen) continue;
      if (rvar[kk] > 0) lincomb2[kk] += lincomb1[k];
       else lincomb2[kk+1] -= lincomb1[k];
    }
  }
}
