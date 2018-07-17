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
  int p;
  struct presentationlist pstlist;
  int results[1];

  p = globals.p;
  pstlist.p = pst;
  pstlist.next = 0;
  count_sl2zp_cclasses (&pstlist, p, results);

  if (quiet) printf ("%d\n", *results);
   else printf ("Result: %d\n", *results);

  return (1);
}

int
cccountsl2zp_list (struct presentationlist *pstlist)
{
  struct presentationlist *pstl;
  int i, count, p;
  int *results;

  p = globals.p;
  for (count = 0, pstl = pstlist; pstl; pstl = pstl->next, count++);
  results = (int *) malloc (count*sizeof(int));

  count_sl2zp_cclasses (pstlist, p, results);

  for (i = 0; i < count; i++)
  {
    if (quiet) printf ("%d\n", results[i]);
     else printf ("Result: %d\n", results[i]);
  }

  return (1);
}

/*
 * computing the number of conjugate classes of homomorphisms from the
 * fundamental group into SL(2,Z/pZ)
 */

void
count_sl2zp_cclasses (struct presentationlist *pstlist, int p, int *results)
{
  struct sl2elem sl2vec[MAXGENNUM];
  struct sl2elem sl2vecinv[MAXGENNUM];
  struct presentationlist *pstl;
  int i, k, gennum, rem;
  int clevel;

  gennum = pstlist->p->gennum;
  assert (gennum <= MAXGENNUM);

  if (!isprime (p))
  {
    if (globals.insist)
    {
      fprintf (stderr, "Warning p = %d is not prime\n", p);
    } else {
      printf ("FATAL: p = %d must be a prime number\n", p);
      printf ("if you insist in doing this computation add option '--insist'\n");
      exit (111);
    }
  }

  /* we need space for a vector of gennum matrices and their inverses */

  for (i = 0; i < gennum; i++)
  {
    sl2_clear (sl2vec[i].a);
    rem = sl2_next_det1 (sl2vec[i].a, p);
    assert (rem == 0);
    sl2_invert (sl2vec[i].a, sl2vecinv[i].a, p);
  }

  for (k = 0, pstl = pstlist; pstl; pstl = pstl->next, k++)
  {
    results[k] = 0;
    if (pstl->p->gennum != gennum)
    {
      printf ("FATAL: all groups must have the same number of generators!\n");
      exit (123);
    }
  }

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
    for (k = 0, pstl = pstlist; pstl; pstl = pstl->next, k++)
    {
      if (sl2_checkrelators(sl2vec, sl2vecinv, pstl->p, p))
      {
        results[k]++;
        if (verbose)
        {
          if (k == 0) printf ("Homomorphism #%d defined by the matrices:\n", results[k]);
           else printf ("Homomorphism #%d for group #%d defined by the matrices:\n", results[k], k);
          for (i = 0; i < gennum; i++)
          {
            sl2_print (sl2vec[i].a);
            printf ("---\n");
          }
        }
      }
    }
  } while (sl2_nextmap (sl2vec, sl2vecinv, gennum, p) == 0);
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

  if (globals.dontidentify) return (0); /* user request is to not identify by conjugacy */
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
 * ======= PSL(2,q) =================================
 */

int
cccountpsl2q (struct presentation *pst)
{
  int q;
  struct presentationlist pstlist;
  int results[1];

  q = globals.q;
  pstlist.p = pst;
  pstlist.next = 0;
  count_psl2q_cclasses (&pstlist, q, results);

  if (quiet) printf ("%d\n", *results);
   else printf ("Result: %d\n", *results);

  return (1);
}

int
cccountpsl2q_list (struct presentationlist *pstlist)
{
  struct presentationlist *pstl;
  int i, count, q;
  int *results;

  q = globals.q;
  for (count = 0, pstl = pstlist; pstl; pstl = pstl->next, count++);
  results = (int *) malloc (count*sizeof(int));

  count_psl2q_cclasses (pstlist, q, results);

  for (i = 0; i < count; i++)
  {
    if (quiet) printf ("%d\n", results[i]);
     else printf ("Result: %d\n", results[i]);
  }

  return (1);
}

/*
 * computing the number of conjugate classes of homomorphisms from the
 * fundamental group into PSL(2,q) (based on globals.q)
 */

void
count_psl2q_cclasses (struct presentationlist *pstlist, int q, int *results)
{
  struct sl2elem sl2vec[MAXGENNUM];
  struct sl2elem sl2vecinv[MAXGENNUM];
  struct presentationlist *pstl;
  int i, k, gennum, rem;
  int clevel;

  gennum = pstlist->p->gennum;
  assert (gennum <= MAXGENNUM);

  if (q <= 2 || !isprime (q))
  {
    if (globals.insist)
    {
      fprintf (stderr, "Warning p = %d must be an odd prime\n", q);
    } else {
      printf ("FATAL: p = %d must be an odd prime number\n", q);
      printf ("if you insist in doing this computation add option '--insist'\n");
      exit (111);
    }
  }

  /* we need space for a vector of gennum matrices and their inverses */

  for (i = 0; i < gennum; i++)
  {
    sl2_clear (sl2vec[i].a);
    rem = psl2_next_det1 (sl2vec[i].a, q);
    assert (rem == 0);
    sl2_invert (sl2vec[i].a, sl2vecinv[i].a, q);
  }

  for (k = 0, pstl = pstlist; pstl; pstl = pstl->next, k++)
  {
    results[k] = 0;
    if (pstl->p->gennum != gennum)
    {
      printf ("FATAL: all groups must have the same number of generators!\n");
      exit (123);
    }
  }

  do {
    if ((clevel = psl2_isnotcanon(sl2vec, gennum, q)))
    {
      if (clevel < gennum)
      {
        /*
         * we can prune away a whole branch of homomorphisms
         * by setting the remaining matrices to the lexicographically "last" value
         */
        for (i = clevel; i < gennum; i++)
        {
          sl2_set (sl2vec[i].a, q);
        }
      }
      continue;
    }
    for (k = 0, pstl = pstlist; pstl; pstl = pstl->next, k++)
    {
      if (psl2_checkrelators(sl2vec, sl2vecinv, pstl->p, q))
      {
        results[k]++;
        if (verbose)
        {
          if (k == 0) printf ("Homomorphism #%d defined by the matrices:\n", results[k]);
           else printf ("Homomorphism #%d for group #%d defined by the matrices:\n", results[k], k);
          for (i = 0; i < gennum; i++)
          {
            sl2_print (sl2vec[i].a);
            printf ("---\n");
          }
        }
      }
    }
  } while (psl2_nextmap (sl2vec, sl2vecinv, gennum, q) == 0);
}

/*
 * check if this set is not a canonical representative up to conjugation
 */

int
psl2_isnotcanon (struct sl2elem *sl2vec, int gennum, int q)
{
  int g[2][2];
  int ginv[2][2];
  int acc[2][2];
  int i, cmpres;

  if (globals.dontidentify) return (0); /* user request is to not identify by conjugacy */
  sl2_clear (g);
  /* cycle on the elements of PSL(2,q) */
  while (psl2_next_det1 (g, q) == 0)
  {
    sl2_invert (g, ginv, q);
    /* now cycle on the elements of the set and compute the conjugate */
    for (i = 0; i < gennum; i++)
    {
      sl2_matcopy (acc, ginv);
      sl2_matmul (acc, sl2vec[i].a, q);
      sl2_matmul (acc, g, q);
      psl2_canon (acc, q);
      cmpres = sl2_matcmp (acc, sl2vec[i].a);
      if (cmpres < 0) return (i+1); /* found a better representative: non canonical */
      if (cmpres > 0) break;      /* it is needless to conjugate the other elements */
    }
  }
  return (0);
}

/*
 * lexicographically next 2x2 matrix mod q with det=1, modulo its sign
 */

int
psl2_next_det1 (int m[2][2], int q)
{
  while (psl_next (2, m, q) == 0)
  {
    if (sl2_det (m, q) == 1) return (0);
  }
  do {
    if (sl2_det (m, q) == 1) return (1);
  } while (psl_next (2, m, q) == 0);
  printf ("FATAL: Cannot find a nondegenerate matrix\n");
  exit (50);
}

/*
 * lexicographically next 2x2 matrix mod p, modulo its sign
 */

int
psl_next (int n, int m[2][2], int q)
{
  int i, j, fi, fj;
  int firstnzero = 0;

  assert (n == 2);  /* otherwise the final if would be different */
  assert (q > 2 && (q % 2) == 1);
  fi = -1;
  for (i = 0; i < n && firstnzero == 0; i++)
  {
    for (j = 0; j < n && firstnzero == 0; j++)
    {
      firstnzero = m[i][j];
      if (firstnzero == 0) continue;
      fi = i;
      fj = j;
    }
  }
  for (i = n-1; i >= 0; i--)
  {
    for (j = n-1; j >= 0; j--)
    {
      m[i][j]++;
      if (m[i][j] < (q+1)/2) return (0);
      if ((i != fi || j != fj) && m[i][j] < q) return (0);
      m[i][j] = 0;
    }
  }
  return (1);
}

/*
 * change sign of a 2x2 matrix if it is "negative"
 */

void
psl2_canon (int m[2][2], int q)
{
  int i, j;
  int n = 2;
  int firstnzero = 0;

  assert (q > 2 && (q % 2) == 1);

  for (i = 0; i < n && firstnzero == 0; i++)
  {
    for (j = 0; j < n && firstnzero == 0; j++)
    {
      firstnzero = m[i][j];
    }
  }

  if (firstnzero == 0)
  {
    fprintf (stderr, "Zero matrix! This should not happen\n");
    return;
  }
  if (firstnzero < (q+1)/2) return;
  /* change sign */
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (m[i][j]) m[i][j] = q - m[i][j];
    }
  }
}

int
psl2_nextmap (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, int gennum, int q)
{
  int rem;

  if (gennum <= 0) return (1);       /* cycling wrapped */
  rem = psl2_nextmap (sl2vec + 1, sl2vecinv + 1, gennum - 1, q);
  if (rem == 0) return (0);          /* found next configuration */
  /* rem == 1 means that the subsequent elements cycled */
  rem = psl2_next_det1 (sl2vec[0].a, q);
  sl2_invert (sl2vec[0].a, sl2vecinv[0].a, q);
  return (rem);
}

int
psl2_checkrelators (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, struct presentation *pst, int q)
{
  struct presentationrule *rule;

  for (rule = pst->rules; rule; rule = rule->next)
    if (psl2_checkrelator (sl2vec, sl2vecinv, pst->gennum, rule, q) == 0) return (0);
  return (1);
}

int
psl2_checkrelator (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, int gennum,
                  struct presentationrule *rule, int q)
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
      sl2_matmul (partial, sl2vec[var-1].a, q);
     else
      sl2_matmul (partial, sl2vecinv[-var-1].a, q);
  }
  if (partial[0][1] != 0) return (0);
  if (partial[1][0] != 0) return (0);
  if (partial[0][0] != partial[1][1]) return (0);
  if (partial[0][0] == 1) return (1);
  if (q - partial[0][0] != 1) return (0);
  return (1);
}

/*
 * ======= SYMMETRIC AND ALTERNATING GROUPS =========
 */

int
cccountsn (struct presentation *pst)
{
  int n;
  struct presentationlist pstlist;
  int results[1];

  n = globals.n;
  pstlist.p = pst;
  pstlist.next = 0;
  count_sn_cclasses (&pstlist, n, results);

  if (quiet) printf ("%d\n", *results);
   else printf ("Result: %d\n", *results);

  return (1);
}

int
cccountsn_list (struct presentationlist *pstlist)
{
  struct presentationlist *pstl;
  int i, count, n;
  int *results;

  n = globals.n;
  for (count = 0, pstl = pstlist; pstl; pstl = pstl->next, count++);
  results = (int *) malloc (count*sizeof(int));

  count_sn_cclasses (pstlist, n, results);

  for (i = 0; i < count; i++)
  {
    if (quiet) printf ("%d\n", results[i]);
     else printf ("Result: %d\n", results[i]);
  }

  return (1);
}

/*
 * computing the number of conjugate classes of homomorphisms from the
 * fundamental group into Sn or An (based on globals.onlyeven
 */

void
count_sn_cclasses (struct presentationlist *pstlist, int n, int *results)
{
  struct snelem snvec[MAXGENNUM];
  struct snelem snvecinv[MAXGENNUM];
  struct presentationlist *pstl;
  int i, k, gennum, rem;
  int clevel;

  gennum = pstlist->p->gennum;
  assert (gennum <= MAXGENNUM);
  assert (n <= REPR_MAXN);

  /* we need space for a vector of gennum permutations  and their inverses */

  for (i = 0; i < gennum; i++)
  {
    sn_init (snvec + i, n);
    rem = sn_next_cond (snvec[i].perm, snvec[i].n);
    assert (rem == 1);
    snvecinv[i].n = n;
    sn_invert (snvec[i].perm, snvecinv[i].perm, n);
  }

  for (k = 0, pstl = pstlist; pstl; pstl = pstl->next, k++)
  {
    results[k] = 0;
    if (pstl->p->gennum != gennum)
    {
      printf ("FATAL: all groups must have the same number of generators!\n");
      exit (123);
    }
  }

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
          sn_setlast (snvec[i].perm, n);
        }
      }
      continue;
    }
    for (k = 0, pstl = pstlist; pstl; pstl = pstl->next, k++)
    {
      if (sn_checkrelators(snvec, snvecinv, pstl->p, n))
      {
        results[k]++;
        if (verbose)
        {
          if (k == 0) printf ("Homomorphism #%d defined by the permutations:\n", results[k]);
           else printf ("Homomorphism #%d for group #%d defined by the permutations:\n", results[k], k);
          for (i = 0; i < gennum; i++)
          {
            sn_print (snvec[i].perm, n);
          }
        }
      }
    }
  } while (sn_nextmap (snvec, snvecinv, gennum) == 0);
}

/*
 * check if a permutation is even
 */

int
sn_iseven (int *perm, int n)
{
  int i, j, iseven;

  iseven = 1;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < i; j++)
    if (perm[j] > perm[i]) iseven = 1 - iseven;
  }
  return (iseven);
}

/*
 * initialize with the permutation that reverts the order.
 * This is the last permutation in lexicographic ordering, so that
 * when we ask for the "next" permutation we should get the
 * identity permutation (the first in lexicographic order) with
 * a "rem" = 1
 */

void
sn_init (struct snelem *perm, int n)
{
  int i;

  assert (n <= REPR_MAXN);
  perm->n = n;
  for (i = 0; i < n; i++) perm->perm[i] = n - i - 1;
}

int
sn_next_cond (int *perm, int n)
{
  while (sn_next (perm, n) == 0)
  {
    if (globals.onlyeven == 0) return (0);
    if (sn_iseven (perm, n)) return (0);
  }
  do {
    if (globals.onlyeven == 0) return (1);
    if (sn_iseven (perm, n)) return (1);
  } while (sn_next (perm, n) == 0);
  printf ("FATAL: Cannot find a suitable permutation\n");
  exit (50);
}

int
sn_next (int *perm, int n)
{
  int rem, i, first;

  if (n <= 1) return (1);
  rem = sn_next (perm + 1, n-1);
  if (rem == 0) return (0);
  /*
   * rem == 1 means that the remainder of the permutation is
   * now ordered in increasing order; I have something like
   * in the example
   *   6 1 3 4 7 9
   * to be transformed in
   *   7 1 3 4 6 9
   * returning 0
   * or
   *   9 1 3 4 6 7
   * to be transformed in
   *   1 3 4 6 7 9
   * returning 1
   */
  first = perm[0];
  for (i = 1; i < n; i++)
  {
    if (perm[i] > first)
    {
      perm[0] = perm[i];
      perm[i] = first;
      return (0);
    }
  }
  for (i = 0; i < n-1; i++)
  {
    perm[i] = perm[i+1];
  }
  perm[n-1] = first;
  return (1);
}

/*
 * of course this works under the assumption that
 * perm->perm[i] is a permutation
 */

void
sn_invert (int *perm, int *perminv, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    perminv[perm[i]] = i;
  }
}

/*
 * multiplication of permutations (left to right)
 * result copied in first argument
 */

void
sn_permmul (int *acc, int *perm, int n)
{
  int i;
  int result[REPR_MAXN];

  assert (n <= REPR_MAXN);

  for (i = 0; i < n; i++) result[i] = perm[acc[i]];
  for (i = 0; i < n; i++) acc[i] = result[i];
}

/*
 * check if this is the representative (minimal with respect to lexicographic order)
 * up to conjugation
 */

int
sn_isnotcanon (struct snelem *perms, int gennum, int n)
{
  int i, j, cmpres;
  int g[REPR_MAXN];
  int ginv[REPR_MAXN];
  int acc[REPR_MAXN];

  if (globals.dontidentify) return (0); /* user request is to not identify by conjugacy */

  assert (n < REPR_MAXN);
  for (j = 0; j < n; j++) g[j] = j;

  /* cycle on the elements of S_n */
  do {
    sn_invert (g, ginv, n);  // WARNING: should we consider only even permutations?

    /* now cycle on the elements of the set and compute the conjugate */
    for (i = 0; i < gennum; i++)
    {
      for (j = 0; j < n; j++) acc[j] = ginv[j];
      sn_permmul (acc, perms[i].perm, n);
      sn_permmul (acc, g, n);
      cmpres = sn_permcmp (acc, perms[i].perm, n);
      if (cmpres < 0) return (i+1); /* found a better representative: non canonical */
      if (cmpres > 0) break;      /* it is needless to conjugate the other elements */
    }
  /* chose whether conjugacy should be done in the larger group S_n or in the smaller A_n */
  } while (sn_next_conj (g, n) == 0);
  return (0);
}

/*
 */

int
sn_next_conj (int *perm, int n)
{
  if (globals.onlyeven == 0) return (sn_next (perm, n));
  if (globals.inner) return (sn_next_cond (perm, n));
  return (sn_next (perm, n));
}

/*
 * lexicographic comparison between two permutations
 */

int
sn_permcmp (int *p1, int *p2, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    if (p1[i] < p2[i]) return (-1);
    if (p1[i] > p2[i]) return (1);
  }
  return (0);
}

void
sn_setlast (int *perm, int n)
{
  int i;

  for (i = 0; i < n; i++) perm[i] = n - i - 1;
}

int
sn_checkrelators (struct snelem *perms, struct snelem *permsinv,
                  struct presentation *pst, int n)
{
  struct presentationrule *rule;

  for (rule = pst->rules; rule; rule = rule->next)
    if (sn_checkrelator (perms, permsinv, pst->gennum, rule, n) == 0) return (0);
  return (1);
}

int
sn_checkrelator (struct snelem *perms, struct snelem *permsinv, int gennum,
                 struct presentationrule *rule, int n)
{
  int i, var;
  int partial[REPR_MAXN];

  for (i = 0; i < n; i++) partial[i] = i;

  for (i = 0; i < rule->length; i++)
  {
    var = rule->var[i];
    assert (var);
    if (var > 0)
      sn_permmul (partial, perms[var-1].perm, n);
     else
      sn_permmul (partial, permsinv[-var-1].perm, n);
  }
  for (i = 0; i < n; i++)
  {
    if (partial[i] != i) return (0);
  }
  return (1);
}

void
sn_print (int *perm, int n)
{
  int i, j, isidentity;
  int mark[REPR_MAXN];

  assert (n <= REPR_MAXN);
  printf ("[");
  for (i = 0; i < n; i++)
  {
    if (i > 0) printf (" ");
    printf ("%d", perm[i] + 1);
  }
  printf ("] --> ");

  for (i = 0; i < n; i++) mark[i] = 0;

  isidentity = 1;
  for (i = 0; i < n; i++)
  {
    if (mark[i]) continue;
    if (perm[i] == i) continue;
    printf ("(");
    isidentity = 0;
    j = i;
    do {
      if (i != j) printf (" ");
      mark[j] = 1;
      printf ("%d", j + 1);
      j = perm[j];
    } while (perm[j] != perm[i]);
    printf (")");
  }
  if (isidentity) printf ("()");
  printf ("\n");
}

int
sn_nextmap (struct snelem *perms, struct snelem *permsinv, int gennum)
{
  int rem;

  if (gennum <= 0) return (1);       /* cycling wrapped */
  rem = sn_nextmap (perms + 1, permsinv + 1, gennum - 1);
  if (rem == 0) return (0);          /* found next configuration */
  /* rem == 1 means that the subsequent elements cycled */
  rem = sn_next_cond (perms[0].perm, perms[0].n);
  sn_invert (perms[0].perm, permsinv[0].perm, perms[0].n);
  return (rem);
}

/*
 * various utilities
 */

/*
 * this is crude!
 * however we do not expect to use this for large values of p
 */

int
isprime (int p)
{
  int a;

  if (p < 2) return (0);

  for (a = 2; a*a <= p; a++)
  {
    if ((p % a) == 0) return (0);
  }

  return (1);
}

int
isprimepower (int q)
{
  printf ("NOT IMPL\n");
  exit (100);
}
