/*
 * computation of the jacobian of the presentation
 * (transformed in a preabelian form)
 * using Fox calculus.
 * The result is then mapped through the abelianizer map
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "contour.h"
#include <fundamental.h>
#include <laurent.h>
#include <alexander.h>
#include <fox.h>

void
foxjacobian (struct presentation *pr)
{
  extern struct global_data globals;
  extern int verbose, interactive;
  struct presentationrule *r;
  struct laurentpoly *p;
  int rank = 1, i, j, e, numrelators = 0;
  int maxlen, rid, k, *lincomb1, *lincomb2, gen;

  if (globals.abelianize == 1 && globals.preabelian == 0)
  {
    printf ("In order to abelianize we need a preabelian presentation, forcing its computation...\n");
    printf ("You can avoid this message by including the --preabelian option\n");
    globals.preabelian = 1;
  }
  if (globals.preabelian) {topreabelian (pr); rank = compute_fg_rank (pr);}
  if (globals.abelianize == 2) rank = 0;
  if (verbose)
  {
    print_presentation (pr);
    if (globals.preabelian) printf ("rank = %d\n", rank);
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
      if (globals.abelianize == 0) print_groupring_el (r->var, lincomb2, r->length);
      if (globals.abelianize == 0) continue;
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
      if (globals.preabelian) printf ("rank = %d\n", rank);
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
        if (globals.preabelian)
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

/*
 *
 */

struct alexanderideal *
three_components_link (struct presentation *p)
{
  struct alexanderideal *ai;
  struct laurentmatrix *jacobian;
  struct laurentpoly *ltemp, *rest, *rest2, **column, **av, **aw, **bw; /* av is in one-less indet! */
  struct laurentpoly *bj, *cj, *bmod, *ltemp2, *detu, *detv, *detw, *res;
  struct presentationrule *r;
  int *lincomb1, *lincomb2;
  int numrelators = 0, maxlen = 0;
  int i, j, k, aminexpon, rank = 3;

  /* compute the first (n-3) columns of the jacobian (if any) */

  for (r = p->rules; r; r = r->next) numrelators++;

  /* the presentation is assumed to be preabelian, we compute the fox derivatives
   * of the relators with respect to the first n-3 generators
   */

  for (i = 0, r = p->rules; r; i++, r = r->next) if (r->length > maxlen) maxlen = r->length;
  lincomb1 = (int *) malloc ((maxlen+1)*sizeof (int));
  lincomb2 = (int *) malloc ((maxlen+1)*sizeof (int));
  jacobian = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
  jacobian->numrows = numrelators;
  jacobian->numcols = p->gennum;
  jacobian->columns = (struct laurentpoly ***) malloc (p->gennum * sizeof (struct laurentpoly **));
  assert (p->gennum == numrelators + 1);
  assert (p->gennum >= 3);
  for (i = 0; i < p->gennum - 3; i++)
  {
    jacobian->columns[i] = column = (struct laurentpoly **) malloc (numrelators * sizeof (struct laurentpoly *));
    for (j = 0, r = p->rules; r; r = r->next, j++)
    {
      for (k = 0; k < r->length; k++) lincomb1[k] = 0;
      lincomb1[r->length] = 1;
      foxderivative (r->length, r->var, lincomb1, lincomb2, i+1);
      column[j] = map_to_abelian (r->var, lincomb2, r->length, p->gennum - rank, rank);
    }
  }

  /*
   * computing Av and Aw column.  Av is actually composed of polynomials in 2 indets
   */

  aw = (struct laurentpoly **) malloc (numrelators * sizeof (struct laurentpoly *));
  av = (struct laurentpoly **) malloc (numrelators * sizeof (struct laurentpoly *));
  bw = (struct laurentpoly **) malloc (numrelators * sizeof (struct laurentpoly *));
  for (j = 0, r = p->rules; r; r = r->next, j++)
  {
    for (k = 0; k < r->length; k++) lincomb1[k] = 0;
    lincomb1[r->length] = 1;
    foxderivative (r->length, r->var, lincomb1, lincomb2, p->gennum - 2);
    ltemp = map_to_abelian (r->var, lincomb2, r->length, p->gennum - rank, rank);
    aminexpon = 0;
    if (ltemp) aminexpon = ltemp->minexpon;
    aw[j] = laurent_divide_by_1minusw (ltemp, &rest);
    av[j] = 0;
    if (rest)
    {
      av[j] = laurent_divide_by_1minusw (rest, &rest2);
      free_laurentpoly (rest);
      assert (rest2 == 0);
      /*
       * remember that av[j] should be multiplied by w^e, with e = a[j]->minexpon
       */
    }
    /* bw is computed by dividing B+Av(1-u) by (1-w) */
    for (k = 0; k < r->length; k++) lincomb1[k] = 0;
    lincomb1[r->length] = 1;
    foxderivative (r->length, r->var, lincomb1, lincomb2, p->gennum - 1);
    bj = map_to_abelian (r->var, lincomb2, r->length, p->gennum - rank, rank);
    ltemp = laurent_dup (av[j]);
    laurent_mulu (ltemp);
    ltemp2 = laurent_add_scal (av[j],ltemp, -1); /* one less indet */
    free_laurentpoly (ltemp);
    ltemp = 0;
    if (ltemp2)
    {
      ltemp = (struct laurentpoly *) malloc (POLYSIZE (1));
      ltemp->indets = 3;
      ltemp->minexpon = aminexpon;
      ltemp->stemdegree = 0;
      ltemp->stem[0].lx = ltemp2;
    }
    bmod = laurent_add (bj, ltemp);
    free_laurentpoly (ltemp);
    bw[j] = laurent_divide_by_1minusw (bmod, &rest);
    assert (rest == 0);

    /* now check that c is cu(1-u) + cv(1-v) with cu = -aw and cv = -bw */
    for (k = 0; k < r->length; k++) lincomb1[k] = 0;
    lincomb1[r->length] = 1;
    foxderivative (r->length, r->var, lincomb1, lincomb2, p->gennum);
    cj = map_to_abelian (r->var, lincomb2, r->length, p->gennum - rank, rank);
    ltemp = laurent_dup (aw[j]);
    laurent_mulu (ltemp);
    ltemp2 = laurent_add_scal (aw[j], ltemp, -1);
    free_laurentpoly (ltemp);
    ltemp = laurent_add (cj, ltemp2);
    free_laurentpoly (cj);
    cj = ltemp;   /* this is C - Cu(1-u) */
    ltemp = laurent_dup (bw[j]);
    laurent_mulv (ltemp);
    ltemp2 = laurent_add_scal (bw[j], ltemp, -1);
    free_laurentpoly (ltemp);
    ltemp = laurent_add (cj, ltemp2);
    assert (ltemp == 0);
    free_laurentpoly (cj);

    /* transform av into a 3 indets polynomial */
    if (av[j])
    {
      ltemp = (struct laurentpoly *) malloc (POLYSIZE(1));
      ltemp->indets = 3;
      ltemp->minexpon = aminexpon;
      ltemp->stemdegree = 0;
      ltemp->stem[0].lx = av[j];
      av[j] = ltemp;
    }
  }
  free (lincomb1);
  free (lincomb2);

  jacobian->numcols = numrelators;
  jacobian->columns[numrelators - 2] = av;
  jacobian->columns[numrelators - 1] = aw;
  detu = laurent_compute_determinant (jacobian->columns, numrelators, 3);

  jacobian->columns[numrelators - 1] = bw;
  detv = laurent_compute_determinant (jacobian->columns, numrelators, 3);

  jacobian->columns[numrelators - 2] = aw;
  detw = laurent_compute_determinant (jacobian->columns, numrelators, 3);

  jacobian->columns[numrelators] = av;
  jacobian->numcols = numrelators + 1;
  
  laurent_free_matrix (jacobian);

  res = laurent_dup (detu);
  res = laurent_addto_scal (res, detv, 1);
  res = laurent_addto_scal (res, detw, 1);

  laurent_mulu (detu);
  laurent_mulv (detv);
  detw->minexpon++;

  res = laurent_addto_scal (res, detu, -1);
  res = laurent_addto_scal (res, detv, -1);
  res = laurent_addto_scal (res, detw, -1);

  free_laurentpoly (detu);
  free_laurentpoly (detv);
  free_laurentpoly (detw);

  laurent_canonify_exponents (res);

  ai = (struct alexanderideal *) malloc (AI_DIM(3));  /* to be sure there is space in fl2 */
  ai->max_generators_num = 3;
  ai->indets = 3;
  ai->fl2offset = ai->max_generators_num/2;
  ai->spread = 1;
  ai->l2num = 0;
  ai->fl2num = 1;
  ai->l[ai->fl2offset] = res;
  ai->gcd = 0;

  return(ai);
}

struct alexanderideal *
generic_ideal_computation (struct presentation *p, int indets, int minordim)
{
  extern int verbose, quiet;
  struct alexanderideal *ai;
  struct laurentpoly **column;
  struct presentationrule *r;
  struct laurentmatrix *jacobian;
  int i, j, k;
  int *lincomb1, *lincomb2;
  int numrelators = 0, maxlen = 0;

  /* compute the columns of the jacobian */

  assert (minordim >= 1);
  for (r = p->rules; r; r = r->next) numrelators++;

  /* the presentation is assumed to be preabelian, we compute the fox derivatives
   * of the relators with respect to the generators
   */

  for (i = 0, r = p->rules; r; i++, r = r->next) if (r->length > maxlen) maxlen = r->length;
  lincomb1 = (int *) malloc ((maxlen+1)*sizeof (int));
  lincomb2 = (int *) malloc ((maxlen+1)*sizeof (int));
  jacobian = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
  jacobian->numrows = numrelators;
  jacobian->numcols = p->gennum;
  jacobian->columns = (struct laurentpoly ***) malloc (p->gennum * sizeof (struct laurentpoly **));
  for (i = 0; i < p->gennum; i++)
  {
    jacobian->columns[i] = column = (struct laurentpoly **) malloc (numrelators * sizeof (struct laurentpoly *));
    for (j = 0, r = p->rules; r; r = r->next, j++)
    {
      for (k = 0; k < r->length; k++) lincomb1[k] = 0;
      lincomb1[r->length] = 1;
      foxderivative (r->length, r->var, lincomb1, lincomb2, i+1);
      column[j] = map_to_abelian (r->var, lincomb2, r->length, p->gennum - indets, indets);
    }
  }
  free (lincomb1);
  free (lincomb2);

  if (verbose) print_matrix (jacobian, indets);
  if (!quiet) printf ("# generic ideal computation...\n");

  ai = compute_invariant_factor (jacobian->columns, jacobian->numrows, jacobian->numcols, minordim, indets);
  if (ai)
  {
    for (k = 0; k < ai->l2num; k++)
    {
      laurent_canonify_exponents (ai->l[k]);
    }
  }

  if (ai == 0)
  {
    printf ("Computation of Alexander ideal is not yet implemented for %d by %d matrices in %d indeterminates, ",
            jacobian->numrows, jacobian->numcols, indets);
    printf ("invariant factor: %d\n", minordim);
  }
  laurent_free_matrix (jacobian);
  return (ai);
}

