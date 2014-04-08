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
  extern int verbose, quiet;
  struct presentationrule *r;
  struct laurentpoly *l, *determinant;
  int i, j, sum, rank, matrixrank;
  int numcols = 0;
  int gconj;
  struct laurentpoly **matrixcolumn;
  struct laurentpoly ***matrix;

  topreabelian (p);

  for (r = p->rules; r; r = r->next) numcols++;
  rank = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    assert (sum >= 0);
    if (sum) matrixrank = i;
      else gconj = i;
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
  gconj = p->gennum;

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank != numcols)
  {
    printf ("Cannot compute Alexander polynomial because there are too many relators\n");
    return (0);
  }

  matrix = (struct laurentpoly ***) malloc (numcols*sizeof(struct laurentpoly **));
  for (j = 0; j < numcols; j++)
  {
    matrix[j] = (struct laurentpoly **) malloc (matrixrank*sizeof(struct laurentpoly *));
  }

  if (verbose) print_presentation (p);

  /*
   * fill matrix
   */

  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix[j];
    for (i = 1; i <= matrixrank; i++)
    {
      matrixcolumn[i-1] = laurent_get_exp_sum (r, i, gconj);
    }
  }

  /*
   * stampa della matrice
   */

  if (verbose)
  {
    for (i = 0; i < matrixrank; i++)
    {
      for (j = 0; j < numcols; j++)
      {
        matrixcolumn = matrix[j];
        l = matrixcolumn[i];
        print_laurentpoly (l);
        printf ("; \t");
      }
      printf ("\n");
    }
  }

  determinant = laurent_compute_determinant (matrix, matrixrank);

  laurent_canonify (determinant);
  if (!quiet) printf ("Alexander polynomial (up to t -> 1/t):\n");
  print_laurentpoly (determinant);
  printf ("\n");

  /*
   * free allocated space
   */

  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    for (i = 0; i < matrixrank; i++)
    {
      l = matrixcolumn[i];
      if (l) free (l);
    }
    free (matrixcolumn);
  }
  free (matrix);

  return (1);
}

/*
 *
 */

struct laurentpoly *
laurent_get_exp_sum (struct presentationrule *r, int n, int gconj)
{
  int k, d;
  int stemdegree;
  int runningexp, minexp, maxexp;
  struct laurentpoly *l;

  runningexp = 0;
  minexp = INT_MAX;
  maxexp = INT_MIN;
  for (k = 0; k < r->length; k++)
  {
    if (abs(r->var[k]) == gconj)
    {
      if (r->var[k] > 0) runningexp++;
        else runningexp--;
    }
    if (abs(r->var[k]) == n)
    {
      if (runningexp < minexp) minexp = runningexp;
      if (runningexp > maxexp) maxexp = runningexp;
    }
  }

  if (minexp > maxexp) return (0);
  stemdegree = maxexp - minexp;
  l = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + (stemdegree + 1)*sizeof(int));
  l->denom = 1;
  l->stemdegree = stemdegree;
  l->minexpon = minexp;
  for (d = 0; d <= stemdegree; d++) l->stem[d] = 0;

  runningexp = 0;
  for (k = 0; k < r->length; k++)
  {
    if (abs(r->var[k]) == gconj)
    {
      if (r->var[k] > 0) runningexp++;
        else runningexp--;
    }
    if (abs(r->var[k]) == n)
    {
      d = runningexp - minexp;
      if (r->var[k] > 0) l->stem[d]++;
        else l->stem[d]--;
    }
  }

  return (l);
}

/*
 * compute the determinant of the matrix
 */

struct laurentpoly *
laurent_compute_determinant (struct laurentpoly ***matrix, int n)
{
  int i, ii, jj, i1;
  int sign;
  struct laurentpoly *determinant = 0, *subdeterminant, *product, *sum;
  struct laurentpoly ***submatrix;
  struct laurentpoly **matrixcol, **submatrixcol, **firstcolumn;

  assert (n >= 0);
  if (n == 0)
  {
    determinant = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + sizeof (int));
    determinant->denom = 1;
    determinant->minexpon = 0;
    determinant->stemdegree = 0;
    determinant->stem[0] = 1;
    return (determinant);
  }
  if (n == 1) return (laurent_dup(matrix[0][0]));

  /*
   * developping the determinant about the first column
   */

  submatrix = (struct laurentpoly ***) malloc ( (n-1)*sizeof (struct laurentpoly **) );
  for (jj = 0; jj < n-1; jj++)
    submatrix[jj] = (struct laurentpoly **) malloc ( (n-1)*sizeof (struct laurentpoly *) );

  firstcolumn = matrix[0];
  sign = 1;
  for (i = 0, sign = 1; i < n; i++, sign = -sign)
  {
    if (firstcolumn[i] == 0) continue;
    for (jj = 0; jj < n - 1; jj++)
    {
      submatrixcol = submatrix[jj];
      matrixcol = matrix[jj + 1];
      for (i1 = 0, ii = 0; i1 < n; i1++)
      {
        if (i1 == i) continue;
        submatrixcol[ii++] = matrixcol[i1];
      }
    }

    subdeterminant = laurent_compute_determinant (submatrix, n-1);
    product = laurent_mul (subdeterminant, firstcolumn[i]);
    if (subdeterminant) free (subdeterminant);
    if (sign < 0) laurent_negate (product);
    sum = laurent_add (determinant, product);
    if (product && product != sum) free (product);
    if (determinant && sum != determinant) free (determinant);
    determinant = sum;
  }

  for (jj = 0; jj < n-1; jj++) free (submatrix[jj]);
  free (submatrix);

  return (determinant);
}

/*
 * negate polynomial
 */

void
laurent_negate (struct laurentpoly *term)
{
  int i;

  if (term == 0) return;
  for (i = 0; i <= term->stemdegree; i++) term->stem[i] = -term->stem[i];

  return;
}

/*
 * canonify with respect to t->1/t
 */

void
laurent_canonify (struct laurentpoly *l)
{
  int k, kk, minexp, maxexp;

printf ("in canonify\n");
  if (l == 0) return;

  minexp = l->minexpon;
  maxexp = minexp + l->stemdegree;

  if (minexp + maxexp < 0)
  {
    laurent_t_to_oneovert (l);
    return;
  }

  /* The degrees are balanced, making lexicographic comparison */

  for (k = 0, kk = l->stemdegree; k < kk; k++, kk--)
  {
    if (l->stem[k] > l->stem[kk])
    {
      laurent_t_to_oneovert (l);
      return;
    }
  }
  return;
}

/*
 * transform t->1/t
 */

void
laurent_t_to_oneovert (struct laurentpoly *l)
{
  int k, kk, saved;

  if (l == 0) return;

  l->minexpon = -(l->minexpon + l->stemdegree);

  for (k = 0, kk = l->stemdegree; k < kk; k++, kk--)
  {
    saved = l->stem[k];
    l->stem[k] = l->stem[kk];
    l->stem[kk] = saved;
  }
  return;
}

/*
 * duplicate a laurent polynomial
 */

struct laurentpoly *
laurent_dup (struct laurentpoly *l)
{
  int k;
  struct laurentpoly *res;

  if (l == 0) return (0);

  res = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) +
            (l->stemdegree + 1)*sizeof(int));

  res->denom = l->denom;
  res->minexpon = l->minexpon;
  res->stemdegree = l->stemdegree;

  for (k = 0; k <= l->stemdegree; k++) res->stem[k] = l->stem[k];

  return (res);
}

/*
 * add two laurent polynomials
 */

struct laurentpoly *
laurent_add (struct laurentpoly *a1, struct laurentpoly *a2)
{
  int minexp, maxexp, k;
  struct laurentpoly *res;

  if (a1 == 0) return (laurent_dup(a2));
  if (a2 == 0) return (laurent_dup(a1));

  minexp = a1->minexpon;
  if (a2->minexpon < minexp) minexp = a2->minexpon;

  maxexp = a1->minexpon + a1->stemdegree;
  if (a2->minexpon + a2->stemdegree > maxexp) maxexp = a2->minexpon + a2->stemdegree;

  res = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) +
           (maxexp - minexp + 1)*sizeof(int));

  assert (a1->denom == a2->denom);
  res->denom = a1->denom;
  res->minexpon = minexp;
  res->stemdegree = maxexp - minexp;

  for (k = 0; k <= res->stemdegree; k++) res->stem[k] = 0;

  for (k = 0; k <= a1->stemdegree; k++) res->stem[k + a1->minexpon - minexp] = a1->stem[k];
  for (k = 0; k <= a2->stemdegree; k++) res->stem[k + a2->minexpon - minexp] += a2->stem[k];

  /* normalize */

  while (res->stem[0] == 0)
  {
    if (res->stemdegree == 0)
    {
      free (res);
      return (0);
    }
    for (k = 0; k < res->stemdegree; k++) res->stem[k] = res->stem[k+1];
    res->stemdegree--;
  }

  while (res->stem[res->stemdegree] == 0)
  {
    assert (res->stemdegree > 0);
    res->stemdegree--;
  }

  return (res);
}

/*
 * multiply two laurent polynomials
 */

struct laurentpoly *
laurent_mul (struct laurentpoly *f1, struct laurentpoly *f2)
{
  int i, j;
  int resstemdegree;
  struct laurentpoly *res;

  if (f1 == 0 || f2 == 0) return (0);

  resstemdegree = f1->stemdegree + f2->stemdegree;

  res = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) +
          (resstemdegree + 1)*sizeof (int));
  res->denom = f1->denom * f2->denom;
  res->minexpon = f1->minexpon + f2->minexpon;
  res->stemdegree = resstemdegree;

  for (i = 0; i <= resstemdegree; i++) res->stem[i] = 0;

  for (i = 0; i <= f1->stemdegree; i++)
  {
    for (j = 0; j <= f2->stemdegree; j++)
    {
      res->stem[i + j] += f1->stem[i]*f2->stem[j];
    }
  }

  return (res);
}

/*
 * print a Laurent polynomial
 */

void
print_laurentpoly (struct laurentpoly *l)
{
  int i, expon;

  if (l == 0) {printf ("0"); return;}
  assert (l->denom == 1);
  for (i = 0; i <= l->stemdegree; i++)
  {
    expon = i + l->minexpon;
    if (l->stem[i])
    {
      if (abs(l->stem[i]) != 1 || expon == 0)
        printf ("%+d", l->stem[i]);
       else
      {
        if (l->stem[i] > 0) printf ("+");
         else printf ("-");
      }
      if (expon != 0)
      {
        printf ("t");
        if (expon != 1)
        {
          if (expon > 0) printf ("^%d", expon);
            else printf ("^(%d)", expon);
        }
      }
    }
  }
}
