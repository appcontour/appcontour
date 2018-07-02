#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "contour.h"
#include "representations.h"

extern int quiet, verbose;
extern struct global_data globals;

int
cccountsl2zp (struct presentation *pst)
{
  int count, p;

  p = globals.p;
  count = count_sl2zp_cclasses (pst, p);

  if (quiet) printf ("%d\n", count);
   else printf ("Result: %d\n", count);

  return (1);
}

/*
 * computing the number of conjugate classes of homomorphisms from the
 * fundamental group into SL(2,Z/pZ)
 */

int
count_sl2zp_cclasses (struct presentation *pst, int p)
{
  struct sl2elem sl2vec[MAXGENNUM];
  struct sl2elem sl2vecinv[MAXGENNUM];
  int i, gennum, rem, count;
  int clevel;

  gennum = pst->gennum;
  assert (gennum <= MAXGENNUM);

  /* we need space for a vector of gennum matrices and their inverses */

  for (i = 0; i < gennum; i++)
  {
    sl2_clear (sl2vec[i].a);
    rem = sl2_next_det1 (sl2vec[i].a, p);
    assert (rem == 0);
    sl2_invert (sl2vec[i].a, sl2vecinv[i].a, p);
  }

  count = 0;

  do {
    if ((clevel = sl2_isnotcanon(sl2vec, gennum, p)))
    {
      if (clevel < gennum)
      {
        /*
         * we can prune away a whole branch of homomorphisms
         * by setting the remaining matrices to the lexicographically "last" value
         */
        for (i = clevel; i < gennum; i++)
        {
          sl2_set (sl2vec[i].a, p);
        }
      }
      continue;
    }
    if (sl2_checkrelators(sl2vec, sl2vecinv, pst, p) == 0) continue;
    count++;
    if (verbose)
    {
      printf ("Homomorphism #%d defined by the matrices:\n", count);
      for (i = 0; i < gennum; i++)
      {
        sl2_print (sl2vec[i].a);
        printf ("---\n");
      }
    }
  } while (sl2_nextmap (sl2vec, sl2vecinv, gennum, p) == 0);

  return (count);
}

/*
 * invert a 2x2 matrix mod p
 */

void
sl2_invert (int m[2][2], int minv[2][2], int p)
{
  int det, detinv;

  det = sl2_det (m, p);
  assert (det != 0);

  detinv = 1;
  if (det != 1) detinv = inv_modp (det, p);

  minv[0][0] = (m[1][1]*detinv) % p;
  minv[1][1] = (m[0][0]*detinv) % p;
  minv[0][1] = ((p - m[0][1])*detinv) % p;
  minv[1][0] = ((p - m[1][0])*detinv) % p;
}

void
sl2_clear (int m[2][2])
{
  int i, j;

  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      m[i][j] = 0;
}

void
sl2_set (int m[2][2], int p)
{
  int i, j;

  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      m[i][j] = p - 1;
}

/*
 * lexicographically next 2x2 matrix mod p with det=1
 */

int
sl2_next_det1 (int m[2][2], int p)
{
  while (sl2_next (m, p) == 0)
  {
    if (sl2_det (m, p) == 1) return (0);
  }
  while (sl2_next (m, p) == 0)
  {
    if (sl2_det (m, p) == 1) return (1);
  }
  printf ("FATAL: Cannot find a nondegenerate matrix\n");
  exit (50);
}

/*
 * lexicographically next 2x2 matrix mod p
 */

int
sl2_next (int m[2][2], int p)
{
  m[1][1]++;
  if (m[1][1] < p) return (0);
  m[1][1] = 0;

  m[1][0]++;

  if (m[1][0] < p) return (0);
  m[1][0] = 0;

  m[0][1]++;

  if (m[0][1] < p) return (0);
  m[0][1] = 0;

  m[0][0]++;

  if (m[0][0] < p) return (0);
  m[0][0] = 0;

  return (1);
}

int
sl2_nextmap (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, int gennum, int p)
{
  int rem;

  if (gennum <= 0) return (1);       /* cycling wrapped */
  rem = sl2_nextmap (sl2vec + 1, sl2vecinv + 1, gennum - 1, p);
  if (rem == 0) return (0);          /* found next configuration */
  /* rem == 1 means that the subsequent elements cycled */
  rem = sl2_next_det1 (sl2vec[0].a, p);
  sl2_invert (sl2vec[0].a, sl2vecinv[0].a, p);
  return (rem);
}

int
sl2_det (int m[2][2], int p)
{
  int ldet = m[0][0]*m[1][1] - m[0][1]*m[1][0];

  ldet = ldet % p;
  if (ldet >= 0) return (ldet);

  return ((ldet + p) % p);
}

int
inv_modp (int n, int p)
{
  int i, x, y;

  assert (p > 0);
  assert (n > 0);

  x = n % p;
  for (i = 0; i < p; i++)
  {
    y = (x*n) % p;
    if (y == 1) return (x);
    x = y;
  }
  printf ("FATAL: Cannot find inverse of %d mod %d\n", n, p);
  exit (55);
}

void
sl2_print (int m[2][2])
{
  printf ("[%d %d]\n", m[0][0], m[0][1]);
  printf ("[%d %d]\n", m[1][0], m[1][1]);
}

int
sl2_checkrelators (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, struct presentation *pst, int p)
{
  struct presentationrule *rule;

  for (rule = pst->rules; rule; rule = rule->next)
    if (sl2_checkrelator (sl2vec, sl2vecinv, pst->gennum, rule, p) == 0) return (0);
  return (1);
}

int
sl2_checkrelator (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, int gennum,
                  struct presentationrule *rule, int p)
{
  int i, var;
  int partial[2][2];

  partial[0][0] = partial[1][1] = 1;
  partial[0][1] = partial[1][0] = 0;

  for (i = 0; i < rule->length; i++)
  {
    var = rule->var[i];
    assert (var);
    if (var > 0)
      sl2_matmul (partial, sl2vec[var-1].a, p);
     else
      sl2_matmul (partial, sl2vecinv[-var-1].a, p);
  }
  if (partial[0][0] != 1) return (0);
  if (partial[1][1] != 1) return (0);
  if (partial[0][1] != 0) return (0);
  if (partial[1][0] != 0) return (0);
  return (1);
}

/*
 * copy matrix
 */

void
sl2_matcopy (int acc[2][2], int mat[2][2])
{
  acc[0][0] = mat[0][0];
  acc[0][1] = mat[0][1];
  acc[1][0] = mat[1][0];
  acc[1][1] = mat[1][1];
}

/*
 * lexicographical comparison of matrices
 */

int
sl2_matcmp (int m1[2][2], int m2[2][2])
{
  int a1, a2;

  a1 = m1[0][0];
  a2 = m2[0][0];
  if (a1 != a2) return ( (a1 < a2)?(-1):1 );

  a1 = m1[0][1];
  a2 = m2[0][1];
  if (a1 != a2) return ( (a1 < a2)?(-1):1 );

  a1 = m1[1][0];
  a2 = m2[1][0];
  if (a1 != a2) return ( (a1 < a2)?(-1):1 );

  a1 = m1[1][1];
  a2 = m2[1][1];
  if (a1 != a2) return ( (a1 < a2)?(-1):1 );

  return (0);
}

/*
 * multiply matrix acc times matrix mat2; result in acc
 */

void
sl2_matmul (int acc[2][2], int mat2[2][2], int p)
{
  int r00, r01, r10, r11;

  r00 = acc[0][0]*mat2[0][0] + acc[0][1]*mat2[1][0];
  r00 = r00 % p;

  r01 = acc[0][0]*mat2[0][1] + acc[0][1]*mat2[1][1];
  r01 = r01 % p;

  r10 = acc[1][0]*mat2[0][0] + acc[1][1]*mat2[1][0];
  r10 = r10 % p;

  r11 = acc[1][0]*mat2[0][1] + acc[1][1]*mat2[1][1];
  r11 = r11 % p;

  acc[0][0] = r00;
  acc[0][1] = r01;
  acc[1][0] = r10;
  acc[1][1] = r11;
}

/*
 * check if this set is not a canonical representative up to conjugation
 */

int
sl2_isnotcanon (struct sl2elem *sl2vec, int gennum, int p)
{
  int g[2][2];
  int ginv[2][2];
  int acc[2][2];
  int i, cmpres;

  sl2_clear (g);
  /* cycle on the elements of SL(2,Zp) */
  while (sl2_next_det1 (g, p) == 0)
  {
    sl2_invert (g, ginv, p);
    /* now cycle on the elements of the set and compute the conjugate */
    for (i = 0; i < gennum; i++)
    {
      sl2_matcopy (acc, ginv);
      sl2_matmul (acc, sl2vec[i].a, p);
      sl2_matmul (acc, g, p);
      cmpres = sl2_matcmp (acc, sl2vec[i].a);
      if (cmpres < 0) return (i+1); /* found a better representative: non canonical */
      if (cmpres > 0) break;      /* it is needless to conjugate the other elements */
    }
  }
  return (0);
}

/*
 * check if this set is a canonical representative up to conjugation
 */

int
sl2_iscanon (struct sl2elem *sl2vec, int gennum, int p)
{
  int g[2][2];
  int ginv[2][2];
  int acc[2][2];
  int i, cmpres;

  sl2_clear (g);
  /* cycle on the elements of SL(2,Zp) */
  while (sl2_next_det1 (g, p) == 0)
  {
    sl2_invert (g, ginv, p);
    /* now cycle on the elements of the set and compute the conjugate */
    for (i = 0; i < gennum; i++)
    {
      sl2_matcopy (acc, ginv);
      sl2_matmul (acc, sl2vec[i].a, p);
      sl2_matmul (acc, g, p);
      cmpres = sl2_matcmp (acc, sl2vec[i].a);
      if (cmpres < 0) return (0); /* found a better representative: non canonical */
      if (cmpres > 0) break;      /* it is needless to conjugate the other elements */
    }
  }
  return (1);
}

/*
 * ======= SIMMETRIC AND ALTERNATING GROUPS =========
 */

int
cccountsn (struct presentation *pst)
{
  int count, n;

  n = globals.n;
  count = count_sn_cclasses (pst, n);

  if (quiet) printf ("%d\n", count);
   else printf ("Result: %d\n", count);

  return (1);
}

/*
 * computing the number of conjugate classes of homomorphisms from the
 * fundamental group into Sn or An (based on globals.onlyeven
 */

int
count_sn_cclasses (struct presentation *pst, int n)
{
  struct snelem snvec[MAXGENNUM];
  struct snelem snvecinv[MAXGENNUM];
  int i, gennum, rem, count;
  int clevel;

  gennum = pst->gennum;
  assert (gennum <= MAXGENNUM);
  assert (n <= REPR_MAXN);

  /* we need space for a vector of gennum permutations  and their inverses */

  for (i = 0; i < gennum; i++)
  {
    sn_init (&snvec[i], n);
    rem = sn_next_cond (snvec + i);
    assert (rem == 0);
    sn_invert (snvec + i, snvecinv + i);
  }

  count = 0;

  do {
    if ((clevel = sn_isnotcanon(snvec, gennum, n)))
    {
      if (clevel < gennum)
      {
        /*
         * we can prune away a whole branch of homomorphisms
         * by setting the remaining permutations to the lexicographically "last" value
         */
        for (i = clevel; i < gennum; i++)
        {
          sn_setlast (snvec + i, n);
        }
      }
      continue;
    }
    if (sn_checkrelators(snvec, snvecinv, pst, n) == 0) continue;
    count++;
    if (verbose)
    {
      printf ("Homomorphism #%d defined by the permutations:\n", count);
      for (i = 0; i < gennum; i++)
      {
        sn_print (snvec + i);
        printf ("---\n");
      }
    }
  } while (sn_nextmap (snvec, snvecinv, gennum, n) == 0);

  return (count);
}

/*
 * initialize with the identity permutation
 */

void
sn_init (struct snelem *perm, int n)
{
  int i;

  assert (n <= REPR_MAXN);
  perm->n = n;
  for (i = 0; i < n; i++) perm->perm[i] = i;
}

int
sn_next_cond (struct snelem *perm)
{
  while (sn_next (&perm->perm[0], perm->n) == 0)
  {
    if (globals.onlyeven == 0) return (0);
    assert (0);  // add iseven check
  }
  while (sn_next (&perm->perm[0], perm->n) == 0)
  {
    if (globals.onlyeven == 0) return (1);
    assert (0);  // add iseven check
  }
  printf ("FATAL: Cannot find a suitable permutation\n");
  exit (50);
}

int
sn_next (int *perm, int n)
{
  assert (0);
}

void
sn_invert (struct snelem *perm, struct snelem *perminv)
{
  assert (0);
}

int
sn_isnotcanon(struct snelem *perms, int gennum, int n)
{
  assert (0);
}

void
sn_setlast (struct snelem *perm, int n)
{
  assert (0);
}

int
sn_checkrelators(struct snelem *perms, struct snelem *permsinv,
                      struct presentation *pst, int n)
{
  assert (0);
}

void
sn_print (struct snelem *perm)
{
  int i;

  printf ("[");
  for (i = 0; i < perm->n; i++)
  {
    if (i > 0) printf (" ");
    printf ("%d", perm->perm[i]);
  }
  printf ("] --> ");
  printf ("== cycle notation not implemented ==\n");
}

int
sn_nextmap (struct snelem *perms, struct snelem *permsinv, int gennum, int n)
{
  assert (0);
}

