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
  struct laurentpoly *determinant;
  struct laurentpoly2 *determinant2;
  int i, sum, rank, matrixrank;
  int numcols = 0;
  int gconj, gconj2;

  topreabelian (p);

  if (verbose) print_presentation (p);

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
  if (rank > 2)
  {
    printf ("Cannot compute Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  if (rank < 1)
  {
    printf ("Cannot compute Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  gconj = p->gennum;
  if (rank >= 2) gconj2 = p->gennum - 1;

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols)
  {
    printf ("Cannot compute Alexander polynomial because there are too many relators\n");
    return (0);
  }

  if (matrixrank > numcols)
  {
    printf ("Cannot compute Alexander polynomial because there are too few relators\n");
    return (0);
  }

  switch (rank)
  {
    case 1:
    determinant = laurent_eliminate_one_indeterminate (p, gconj);
    laurent_canonify (determinant);
    if (!quiet) printf ("Alexander polynomial (up to t -> 1/t):\n");
    print_laurentpoly (determinant, 't');
    printf ("\n");
    break;

    case 2:
    determinant2 = laurent_eliminate_two_indeterminates (p, gconj2, gconj);
    laurent_canonify2 (determinant2);
    printf ("*** Warning: canonization procedure not yet implemented ***\n");
    if (!quiet) printf ("Alexander polynomial (CAVEAT):\n");
    print_laurentpoly2 (determinant2, 'u', 'v');
    printf ("\n");
    break;

    default:
    printf ("Software unable to compute Alexander polynomial with these many indeterminates.\n");
    return (0);
  }

  return (1);
}

/*
 * compute the linking number of a two-component link
 */

int
linkingnumber (struct presentation *p)
{
  extern int verbose, quiet;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  int i, sum, rank, matrixrank;
  int deriv, val, e;
  int numcols = 0;
  int gconj;

  topreabelian (p);

  if (verbose) print_presentation (p);

  for (r = p->rules; r; r = r->next) numcols++;
  rank = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    assert (sum >= 0);
    if (sum) matrixrank = i;
    if (sum && sum != 1)
    {
      printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
      return (0);
    }
  }
  rank = p->gennum - matrixrank;
  if (rank <= 1)
  {
    printf ("Cannot compute linking number for groups with rank %d\n", rank);
    return (0);
  }
  if (rank > 2)
  {
    printf ("Cannot compute linking number for groups with rank %d\n", rank);
    return (0);
  }
  rank = 1;
  matrixrank = p->gennum - rank;
  gconj = p->gennum;

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols)
  {
    printf ("Cannot compute linking number because there are too many relators\n");
    return (0);
  }

  determinant = laurent_eliminate_one_indeterminate (p, gconj);

  //laurent_canonify (determinant);

  if (verbose)
  {
    printf ("Determinant of matrix: ");
    print_laurentpoly (determinant, 't');
    printf ("\n");
  }

  /* linking number is the derivative in t=1 */

  deriv = 0;
  val = 0;
  if (determinant != 0)
  {
    for (i = 0; i <= determinant->stemdegree; i++)
    {
      e = i + determinant->minexpon;
      deriv += e*determinant->stem[i];
      val += determinant->stem[i];
    }
  }
  if (val != 0) printf ("Warning: Value in 1 is %d, should be zero.\n", val);
  deriv = abs(deriv);
  if (quiet) printf ("%d\n", deriv);
   else printf ("Linking number is %d\n", deriv);
  return (1);
}

/*
 * corank one (experimental)
 */

int
corank_one_alexander (struct presentation *p)
{
  extern int verbose, quiet;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  int i, sum, rank, matrixrank;
  int deriv, val, e;
  int numcols = 0;
  int gconj;

  topreabelian (p);

  if (verbose) print_presentation (p);

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
      printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
      return (0);
    }
  }
  rank = p->gennum - matrixrank;
  if (rank <= 1)
  {
    printf ("Cannot compute corank one Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  if (rank > 2)
  {
    printf ("Cannot compute corank one Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  rank = 1;
  matrixrank = p->gennum - rank;
  gconj = p->gennum;

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols)
  {
    printf ("Cannot compute Alexander polynomial because there are too many relators\n");
    return (0);
  }

  determinant = laurent_eliminate_one_indeterminate (p, gconj);

  laurent_canonify (determinant);

  if (verbose)
  {
    printf ("Determinant of matrix (before canonization): ");
    print_laurentpoly (determinant, 't');
    printf ("\n");
  }

  /* calcolo della derivata in t=1 */

  deriv = 0;
  val = 0;
  if (determinant != 0)
  {
    for (i = 0; i <= determinant->stemdegree; i++)
    {
      e = i + determinant->minexpon;
      deriv += e*determinant->stem[i];
      val += determinant->stem[i];
    }
  }
  printf ("Value in 1: %d, deriv = %d\n", val, deriv);
  if (!quiet) printf ("Alexander polynomial (up to t -> 1/t):\n");
  print_laurentpoly (determinant, 't');
  printf ("\n");
  return (1);
}

/*
 * eliminate one indeterminate by conjugacy
 */

struct laurentpoly *
laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate)
{
  extern int verbose;
  struct presentationrule *r;
  struct laurentpoly *l, *determinant;
  struct laurentpoly **matrixcolumn;
  struct laurentpoly ***matrix;
  int numcols = 0;
  int numrows;
  int i, ii, j;

  assert (eliminate >= 1);
  assert (eliminate <= p->gennum);

  for (r = p->rules; r; r = r->next) numcols++;
  numrows = p->gennum - 1;
  assert (numcols <= numrows);

  if (numcols < numrows) return (0);

  matrix = (struct laurentpoly ***) malloc (numcols*sizeof(struct laurentpoly **));
  for (j = 0; j < numcols; j++)
  {
    matrix[j] = (struct laurentpoly **) malloc (numrows*sizeof(struct laurentpoly *));
  }

  /*
   * fill matrix
   */

  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix[j];
    for (i = 1, ii = 0; i <= p->gennum; i++)
    {
      if (i == eliminate) continue;
      matrixcolumn[ii++] = laurent_get_exp_sum (r, i, eliminate);
    }
  }

  /*
   * stampa della matrice
   */

  if (verbose)
  {
    printf ("Matrix entries:\n");
    for (i = 0; i < numrows; i++)
    {
      for (j = 0; j < numcols; j++)
      {
        matrixcolumn = matrix[j];
        l = matrixcolumn[i];
        print_laurentpoly (l, 't');
        printf ("; \t");
      }
      printf ("\n");
    }
  }

  determinant = laurent_compute_determinant (matrix, numcols);

  /*
   * free allocated space
   */

  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    for (i = 0; i < numrows; i++)
    {
      l = matrixcolumn[i];
      if (l) free (l);
    }
    free (matrixcolumn);
  }
  free (matrix);

  return (determinant);
}

/*
 * eliminate two indeterminates by conjugacy
 */

struct laurentpoly2 *
laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2)
{
  extern int verbose;
  struct presentationrule *r;
  struct laurentpoly2 *l, *determinant;
  struct laurentpoly2 **matrixcolumn;
  struct laurentpoly2 ***matrix;
  int numcols = 0;
  int numrows;
  int i, ii, j;

  assert (e1 >= 1 && e2 >= 1);
  assert (e1 <= p->gennum);
  assert (e2 <= p->gennum);
  assert (e1 != e2);

  for (r = p->rules; r; r = r->next) numcols++;
  numrows = p->gennum - 2;
  assert (numcols <= numrows);

  if (numcols < numrows) return (0);

  matrix = (struct laurentpoly2 ***) malloc (numcols*sizeof(struct laurentpoly2 **));
  for (j = 0; j < numcols; j++)
  {
    matrix[j] = (struct laurentpoly2 **) malloc (numrows*sizeof(struct laurentpoly2 *));
  }

  /*
   * fill matrix
   */

  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix[j];
    for (i = 1, ii = 0; i <= p->gennum; i++)
    {
      if (i == e1 || i == e2) continue;
      matrixcolumn[ii++] = laurent_get_exp_sum2 (r, i, e1, e2);
    }
  }

  /*
   * stampa della matrice
   */

  if (verbose)
  {
    printf ("Matrix entries:\n");
    for (i = 0; i < numrows; i++)
    {
      for (j = 0; j < numcols; j++)
      {
        matrixcolumn = matrix[j];
        l = matrixcolumn[i];
        print_laurentpoly2 (l, 'u', 'v');
        printf ("; \t");
      }
      printf ("\n");
    }
  }

  determinant = laurent_compute_determinant2 (matrix, numcols);

  /*
   * free allocated space
   */

  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    for (i = 0; i < numrows; i++)
    {
      l = matrixcolumn[i];
      if (l) free_laurentpoly2 (l);
    }
    free (matrixcolumn);
  }
  free (matrix);

  return (determinant);
}

/*
 * Compute matrix entries for Alexander polynomial computation
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

  assert (l);
  if (laurent_normalize(l) == 0)
  {
    free (l);
    return (0);
  }
  assert (l && l->stem[0]);
  return (l);
}

/*
 * Compute matrix entries for Alexander polynomial computation (two indeterminates)
 */

struct laurentpoly2 *
laurent_get_exp_sum2 (struct presentationrule *r, int n, int e1, int e2)
{
  int k, d;
  int stemdegree;
  int re1, re2, minexp1, minexp2, maxexp1, maxexp2;
  struct laurentpoly2 *l;

  re1 = re2 = 0;
  minexp1 = minexp2 = INT_MAX;
  maxexp1 = maxexp2 = INT_MIN;
  for (k = 0; k < r->length; k++)
  {
    if (abs(r->var[k]) == e1)
    {
      if (r->var[k] > 0) re1++;
        else re1--;
    }
    if (abs(r->var[k]) == e2)
    {
      if (r->var[k] > 0) re2++;
        else re2--;
    }
    if (abs(r->var[k]) == n)
    {
      if (re1 < minexp1) minexp1 = re1;
      if (re2 < minexp2) minexp2 = re2;
      if (re1 > maxexp1) maxexp1 = re1;
      if (re2 > maxexp2) maxexp2 = re2;
    }
  }

  if (minexp1 > maxexp1 || minexp2 > maxexp2) return (0);
  stemdegree = maxexp2 - minexp2;
  l = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2)
            + (stemdegree + 1)*sizeof(struct laurentpoly *));
  l->denom = 1;
  l->stemdegree = stemdegree;
  l->minexpon = minexp2;
  for (d = 0; d <= stemdegree; d++) l->stem[d] = 0;

  re1 = re2 = 0;
  for (k = 0; k < r->length; k++)
  {
    if (abs(r->var[k]) == e1)
    {
      if (r->var[k] > 0) re1++;
        else re1--;
    }
    if (abs(r->var[k]) == e2)
    {
      if (r->var[k] > 0) re2++;
        else re2--;
    }
    if (abs(r->var[k]) == n)
    {
      d = re2 - minexp2;
      if (r->var[k] > 0) l->stem[d] = laurentpoly_addmonom (l->stem[d], re1, 1);
        else l->stem[d] = laurentpoly_addmonom (l->stem[d], re1, -1);
    }
  }

  assert (l);
  if (laurent_normalize2(l) == 0)
  {
    free (l);
    return (0);
  }
  assert (l && l->stem[0]);
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
 * compute the determinant of the matrix (two indeterminates)
 */

struct laurentpoly2 *
laurent_compute_determinant2 (struct laurentpoly2 ***matrix, int n)
{
  int i, ii, jj, i1;
  int sign;
  struct laurentpoly2 *determinant = 0, *subdeterminant, *product, *sum;
  struct laurentpoly2 ***submatrix;
  struct laurentpoly2 **matrixcol, **submatrixcol, **firstcolumn;
  struct laurentpoly *lp1;

  assert (n >= 0);
  if (n == 0)
  {
    determinant = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2) + sizeof (struct laurentpoly *));
    determinant->denom = 1;
    determinant->minexpon = 0;
    determinant->stemdegree = 0;
    lp1 = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + sizeof (int));
    determinant->stem[0] = lp1;
    lp1->denom = 1;
    lp1->minexpon = 0;
    lp1->stemdegree = 0;
    lp1->stem[0] = 1;
    return (determinant);
  }
  if (n == 1) return (laurent_dup2(matrix[0][0]));

  /*
   * developping the determinant about the first column
   */

  submatrix = (struct laurentpoly2 ***) malloc ( (n-1)*sizeof (struct laurentpoly2 **) );
  for (jj = 0; jj < n-1; jj++)
    submatrix[jj] = (struct laurentpoly2 **) malloc ( (n-1)*sizeof (struct laurentpoly2 *) );

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

    subdeterminant = laurent_compute_determinant2 (submatrix, n-1);
    product = laurent_mul2 (subdeterminant, firstcolumn[i]);
    if (subdeterminant) free_laurentpoly2 (subdeterminant);
    if (sign < 0) laurent_negate2 (product);
    sum = laurent_add2 (determinant, product);
    if (product && product != sum) free_laurentpoly2 (product);
    if (determinant && sum != determinant) free_laurentpoly2 (determinant);
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
 * negate polynomial (two indeterminates)
 */

void
laurent_negate2 (struct laurentpoly2 *term)
{
  int i;

  if (term == 0) return;
  for (i = 0; i <= term->stemdegree; i++) laurent_negate (term->stem[i]);

  return;
}

/*
 * canonify with respect to t->1/t
 */

void
laurent_canonify (struct laurentpoly *l)
{
  int k, kk, sign;

  if (l == 0) return;

  assert (l->stem[0]);
  assert (l->stem[l->stemdegree]);
  //if (l->stem[0] < 0) laurent_negate (l);

  sign = 1;
  if (l->stem[l->stemdegree] < 0) sign = -1;

  /* Making lexicographic comparison */

  for (k = 0, kk = l->stemdegree; k < kk; k++, kk--)
  {
    if (l->stem[k] > sign*l->stem[kk])
    {
      laurent_t_to_oneovert (l);
      if (l->stem[0] < 0) laurent_negate (l);
      l->minexpon = 0;
      return;
    }
  }
  l->minexpon = 0; // we are interested only to the stem
  return;
}

/*
 * laurent_canonify2
 * for now simply shift exponents to minimal nonnegative
 */

void
laurent_canonify2 (struct laurentpoly2 *l)
{
  int i, minexpon1;
  struct laurentpoly *l1;

  if (l == 0) return;

  assert (l->stem[0]);

  l->minexpon = 0;

  minexpon1 = (l->stem[0])->minexpon;
  for (i = 0; i <= l->stemdegree; i++)
  {
    l1 = l->stem[i];
    if (l1 == 0) continue;
    if (l1->minexpon < minexpon1) minexpon1 = l1->minexpon;
  }
  for (i = 0; i <= l->stemdegree; i++)
  {
    l1 = l->stem[i];
    if (l1 == 0) continue;
    l1->minexpon -= minexpon1;
  }
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
 * duplicate a laurent polynomial (two indeterminates)
 */

struct laurentpoly2 *
laurent_dup2 (struct laurentpoly2 *l)
{
  int k;
  struct laurentpoly2 *res;

  if (l == 0) return (0);

  res = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2) +
            (l->stemdegree + 1)*sizeof(struct laurentpoly *));

  res->denom = l->denom;
  res->minexpon = l->minexpon;
  res->stemdegree = l->stemdegree;

  for (k = 0; k <= l->stemdegree; k++) res->stem[k] = laurent_dup (l->stem[k]);

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

  if (res && laurent_normalize (res) == 0)
  {
    free (res);
    res = 0;
  }

  return (res);
}

/*
 * add two laurent polynomials (two indeterminates)
 */

struct laurentpoly2 *
laurent_add2 (struct laurentpoly2 *a1, struct laurentpoly2 *a2)
{
  int minexp, maxexp, k;
  struct laurentpoly2 *res;
  struct laurentpoly *addres;

  if (a1 == 0) return (laurent_dup2(a2));
  if (a2 == 0) return (laurent_dup2(a1));

  minexp = a1->minexpon;
  if (a2->minexpon < minexp) minexp = a2->minexpon;

  maxexp = a1->minexpon + a1->stemdegree;
  if (a2->minexpon + a2->stemdegree > maxexp) maxexp = a2->minexpon + a2->stemdegree;

  res = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2) +
           (maxexp - minexp + 1)*sizeof(struct laurentpoly *));

  assert (a1->denom == a2->denom);
  res->denom = a1->denom;
  res->minexpon = minexp;
  res->stemdegree = maxexp - minexp;

  for (k = 0; k <= res->stemdegree; k++) res->stem[k] = 0;

  for (k = 0; k <= a1->stemdegree; k++) res->stem[k + a1->minexpon - minexp] = laurent_dup (a1->stem[k]);
  for (k = 0; k <= a2->stemdegree; k++)
  {
    if (a2->stem[k] == 0) continue;
    addres = laurent_add (res->stem[k + a2->minexpon - minexp], a2->stem[k]);
    if (res->stem[k + a2->minexpon - minexp]) free (res->stem[k + a2->minexpon - minexp]);
    res->stem[k + a2->minexpon - minexp] = addres;
  }

  /* normalize */

  if (res && laurent_normalize2 (res) == 0)
  {
    free (res);
    res = 0;
  }

  return (res);
}

/*
 * normalize a laurent polynomial
 * such that first and last element in
 * stem are nonzero
 * returns 0 if the laurent polynomial
 * is identically zero
 */

struct laurentpoly *
laurent_normalize (struct laurentpoly *l)
{
  int k;

  if (l == 0) return (0);

  while (l->stem[0] == 0)
  {
    if (l->stemdegree == 0)
    {
      return (0);
    }
    for (k = 0; k < l->stemdegree; k++) l->stem[k] = l->stem[k+1];
    l->stemdegree--;
    l->minexpon++;
  }

  while (l->stem[l->stemdegree] == 0)
  {
    assert (l->stemdegree > 0);
    l->stemdegree--;
  }
  return (l);
}

struct laurentpoly2 *
laurent_normalize2 (struct laurentpoly2 *l)
{
  int k;

  if (l == 0) return (0);

  while (l->stem[0] == 0)
  {
    if (l->stemdegree == 0)
    {
      return (0);
    }
    for (k = 0; k < l->stemdegree; k++) l->stem[k] = l->stem[k+1];
    l->stemdegree--;
    l->minexpon++;
  }

  while (l->stem[l->stemdegree] == 0)
  {
    assert (l->stemdegree > 0);
    l->stemdegree--;
  }
  return (l);
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
 * multiply two laurent polynomials (two indeterminates)
 */

struct laurentpoly2 *
laurent_mul2 (struct laurentpoly2 *f1, struct laurentpoly2 *f2)
{
  int i, j;
  int resstemdegree;
  struct laurentpoly2 *res;
  struct laurentpoly *mulres, *addres;

  if (f1 == 0 || f2 == 0) return (0);

  resstemdegree = f1->stemdegree + f2->stemdegree;

  res = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2) +
          (resstemdegree + 1)*sizeof (struct laurentpoly *));
  res->denom = f1->denom * f2->denom;
  res->minexpon = f1->minexpon + f2->minexpon;
  res->stemdegree = resstemdegree;

  for (i = 0; i <= resstemdegree; i++) res->stem[i] = 0;

  for (i = 0; i <= f1->stemdegree; i++)
  {
    for (j = 0; j <= f2->stemdegree; j++)
    {
      mulres = laurent_mul (f1->stem[i], f2->stem[j]);
      addres = laurent_add (res->stem[i + j], mulres);
      if (mulres) free (mulres);
      if (res->stem[i + j]) free (res->stem[i + j]);
      res->stem[i + j] = addres;
    }
  }

  return (res);
}

/*
 * print a Laurent polynomial
 */

void
print_laurentpoly (struct laurentpoly *l, char indet)
{
  int i, expon;

  if (l == 0) {printf ("0"); return;}
  assert (l->denom == 1);
  assert (l->stem[0]);
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
        printf ("%c", indet);
        if (expon != 1)
        {
          if (expon > 0) printf ("^%d", expon);
            else printf ("^(%d)", expon);
        }
      }
    }
  }
}

/*
 * print a Laurent polynomial with two indeterminates
 */

void
print_laurentpoly2 (struct laurentpoly2 *l, char indet1, char indet2)
{
  int i, j, expon, expon1;
  struct laurentpoly *l1;

  if (l == 0) {printf ("0"); return;}
  assert (l->denom == 1);
  assert (l->stem[0]);
  for (i = 0; i <= l->stemdegree; i++)
  {
    expon = i + l->minexpon;
    if ((l1 = l->stem[i]) != 0)
    {
      for (j = 0; j <= l1->stemdegree; j++)
      {
        expon1 = j + l1->minexpon;
        if (l1->stem[j] == 0) continue;
        if (abs(l1->stem[j]) != 1 || ((expon1 == 0)&&(expon == 0)))
        {
          printf ("%+d", l1->stem[j]);
        } else {
          if (l1->stem[j] > 0) printf ("+"); else printf ("-");
        }
        if (expon1 != 0)
        {
          printf ("%c", indet1);
          if (expon1 != 1)
          {
            if (expon1 > 0) printf ("^%d", expon1);
             else printf ("^(%d)", expon1);
          }
        }
        if (expon != 0)
        {
          printf ("%c", indet2);
          if (expon != 1)
          {
            if (expon > 0) printf ("^%d", expon);
             else printf ("^(%d)", expon);
          }
        }
      }
    }
  }
}

/*
 * free data associated to a laurent polynomial with 2 indeterminates
 */

void
free_laurentpoly2 (struct laurentpoly2 *l)
{
  int i;

  if (l == 0) return;

  for (i = 0; i <= l->stemdegree; i++)
  {
    if (l->stem[i]) free (l->stem[i]);
  }
  free (l);
}

/*
 * add a monomial to a laurent_poly
 */

struct laurentpoly *
laurentpoly_addmonom (struct laurentpoly *l, int expon, int coef)
{
  struct laurentpoly *m = 0;
  struct laurentpoly *res;

  if (coef == 0) return (l);

  if (l && expon >= l->minexpon && expon <= l->minexpon + l->stemdegree)
  {
    l->stem[expon - l->minexpon] += coef;
    if (laurent_normalize(l) == 0)
    {
      free (l);
      return (0);
    }
    return (l);
  }

  m = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + sizeof (int));
  m->denom = 1;
  m->minexpon = expon;
  m->stemdegree = 0;
  m->stem[0] = coef;

  if (l == 0) return (m);

  res = laurent_add (l, m);
  free (m);
  free (l);

  return (res);
}
