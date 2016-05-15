/*
 * computations with laurent polynomials
 * in one or two indeterminates
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "laurent.h"
#include "parser.h"

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
 * print a Laurent polynomial
 */

void
print_laurentpoly (struct laurentpoly *l, char indet)
{
  int i, expon;

  if (l == 0) {printf ("0"); return;}
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
 * read a laurent polynomial in one indet
 */

struct laurentpoly *
read_laurentpoly (FILE *file, char indet_names[2])
{
  char ch;
  int sign, coef, exp1, exp2;
  struct laurentpoly *l1 = 0;

  ch = mygetchar (file);
  if (ch == ';') return (l1);

  sign = 1;
  if (ch == '-') sign = -1;
  if (ch != '+' && ch != '-') ungetc (ch, file);
  indet_names[1]='*';   // fake a second impossible indet
  while (get_unsignedmonomial2 (file, indet_names, &coef, &exp1, &exp2))
  {
    assert (exp2 == 0);
    coef = sign*coef;
    l1 = laurentpoly_addmonom (l1, exp1, coef);
    ch = mygetchar (file);
    if (ch == ';') return (l1);
    assert (ch == '+' || ch == '-');
    sign = 1;
    if (ch == '-') sign = -1;
  }
  assert (0);
  return (l1);
}

/*
 * read a laurent polynomial in two indets
 */

struct laurentpoly2 *
read_laurentpoly2 (FILE *file, char indet_names[2])
{
  char ch;
  int sign, coef, exp1, exp2;
  struct laurentpoly2 *l2 = 0;

  ch = mygetchar (file);
  if (ch == ';') return (l2);

  sign = 1;
  if (ch == '-') sign = -1;
  if (ch != '+' && ch != '-') ungetc (ch, file);
  while (get_unsignedmonomial2 (file, indet_names, &coef, &exp1, &exp2))
  {
    coef = sign*coef;
    l2 = laurentpoly2_addmonom (l2, exp1, exp2, coef);
    ch = mygetchar (file);
    if (ch == ';') return (l2);
    assert (ch == '+' || ch == '-');
    sign = 1;
    if (ch == '-') sign = -1;
  }
  assert (0);
  return (l2);
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
    free_laurentpoly2 (res);
    res = 0;
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

  res->minexpon = l->minexpon;
  res->stemdegree = l->stemdegree;

  for (k = 0; k <= l->stemdegree; k++) res->stem[k] = laurent_dup (l->stem[k]);

  return (res);
}

/*
 * canonify with respect to t->1/t
 *
 * note that the sign is such that evaluation in 1 gives
 * 1 (should be careful in general)
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
  m->minexpon = expon;
  m->stemdegree = 0;
  m->stem[0] = coef;

  if (l == 0) return (m);

  res = laurent_add (l, m);
  free (m);
  free (l);

  return (res);
}

/*
 * add a monomial to a laurentpoly2
 */

struct laurentpoly2 *
laurentpoly2_addmonom (struct laurentpoly2 *l, int degu, int degv, int coef)
{
  struct laurentpoly2 *m = 0;
  struct laurentpoly *mu = 0;
  struct laurentpoly2 *res;
  int iv = 0;

  if (coef == 0) return (l);

  if (l) iv = degv - l->minexpon;
  if (l && iv >= 0 && iv <= l->stemdegree)
  {
    l->stem[iv] = laurentpoly_addmonom (l->stem[iv], degu, coef);
    if (laurent_normalize2(l) == 0)
    {
      free_laurentpoly2 (l);
      return (0);
    }
    return (l);
  }

  m = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2) + sizeof (struct laurentpoly *));
  mu = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + sizeof (int));
  m->minexpon = degv;
  m->stemdegree = 0;
  m->stem[0] = mu;
  mu->minexpon = degu;
  mu->stemdegree = 0;
  mu->stem[0] = coef;

  if (l == 0) return (m);

  res = laurent_add2 (l, m);
  free_laurentpoly2 (m);
  free_laurentpoly2 (l);

  return (res);
}

/*
 * add up all coefficients of a laurentpoly
 * (equivalent to evaluation in t=1)
 */

int
laurent_sum_coefficients (struct laurentpoly *l)
{
  int k;
  int res = 0;

  if (l == 0) return (0);

  for (k = 0; k <= l->stemdegree; k++) res += l->stem[k];
  return (res);
}

/*
 * add up all coefficients of a laurentpoly2
 * (equivalent to evaluation in v=1)
 */

struct laurentpoly *
laurent_sum_coefficients2 (struct laurentpoly2 *l)
{
  int k;
  struct laurentpoly *res, *addres;

  if (l == 0) return (0);
  res = laurent_dup(l->stem[0]);

  for (k = 1; k <= l->stemdegree; k++)
  {
    addres = laurent_add (res, l->stem[k]);
    if (res) free (res);
    res = addres;
  }
  return (res);
}

/*
 * add up all coefficients of each coefficient or a laurentpoly2
 * (equivalent to evaluation in u=1)
 */

struct laurentpoly *
laurent_sum_each_coefficient2 (struct laurentpoly2 *l)
{
  int k;
  struct laurentpoly *res;

  if (l == 0) return (0);

  res = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + (l->stemdegree + 1)*sizeof(int));
  res->minexpon = l->minexpon;
  res->stemdegree = l->stemdegree;

  for (k = 0; k <= l->stemdegree; k++)
  {
    res->stem[k] = laurent_sum_coefficients (l->stem[k]);
  }
  if (laurent_normalize (res) == 0)
  {
    free (res);
    return (0);
  }

  return (res);
}

/*
 *
 */

struct laurentpoly *
laurent_gcd (struct laurentpoly *p1, struct laurentpoly *p2)
{
  struct laurentpoly *resgcd;
  int i, cont_p1, cont_p2, cont_gcd;

  /* special cases */
  if (p1 == 0)
  {
    return (laurent_dup (p2));
  }
  if (p2 == 0)
  {
    return (laurent_dup (p1));
  }

  /*
   * first: factor out content of the two polynomials
   */

  cont_p1 = laurent_factor_content (p1);
  cont_p2 = laurent_factor_content (p2);

  resgcd = laurent_euclid (p1, p2);
  laurent_factor_content (resgcd);
  cont_gcd = gcd (cont_p1, cont_p2);

  for (i = 0; i <= p1->stemdegree; i++) p1->stem[i] *= cont_p1;
  for (i = 0; i <= p2->stemdegree; i++) p2->stem[i] *= cont_p2;
  for (i = 0; i <= resgcd->stemdegree; i++) resgcd->stem[i] *= cont_gcd;

  return (resgcd);
}

/*
 * minexpon has no effect
 *
 * the resulting gcd is not guaranteed to be primitive
 */

struct laurentpoly *
laurent_euclid (struct laurentpoly *p1, struct laurentpoly *p2)
{
  struct laurentpoly *resgcd;
  int *a, *b, *c;
  int dega, degb, degc, vecsize;
  int i, g, fa, fb, temp;
  int c_is_null;

  //assert (p1->stemdegree >= p2->stemdegree);

  dega = p1->stemdegree;
  degb = p2->stemdegree;
  vecsize = dega;
  if (degb > dega) vecsize = degb;
  vecsize++;
  a = (int *) malloc ( vecsize*sizeof(int) );
  b = (int *) malloc ( vecsize*sizeof(int) );
  c = (int *) malloc ( vecsize*sizeof(int) );
  for (i = 0; i <= dega; i++) a[i] = p1->stem[dega - i];
  for (i = 0; i <= degb; i++) b[i] = p2->stem[degb - i];

  /*
   * the three buffers a, b, c contain the integral coefficients
   * from the higher degree term of the divisor, dividend, rest
   * possibly exchanging a <-> b we assume dega >= degb
   * subtracting from a multiple of 'a' a monomial multiple of 'b'
   * we decrease the degree of 'a' by one
   */
  while (1)
  {
    if (degb > dega)
    {
      for (i = 0; i <= degb; i++)
      {
        temp = a[i];
        a[i] = b[i];
        b[i] = temp;
      }
      temp = dega;
      dega = degb;
      degb = temp;
    }
    for (i = degb + 1; i <= dega; i++) b[i] = 0;
    g = gcd (a[0], b[0]);
    fa = b[0]/g;
    fb = a[0]/g;
    c_is_null = 1;
    for (i = 1; i <= dega; i++)
    {
      c[i-1] = fa*a[i] - fb*b[i];
      if (c[i-1] != 0) c_is_null = 0;
    }
    degc = dega - 1;
    if (c_is_null) break;
    while (c[0] == 0)
    {
      for (i = 0; i < degc; i++) c[i] = c[i+1];
      degc--;
    }
    /*
     * b becomes the new a, c becomes the new b
     */
    dega = degb;
    for (i = 0; i <= dega; i++) a[i] = b[i];
    degb = degc;
    for (i = 0; i <= degb; i++) b[i] = c[i];
  }

  /*
   * c is identically zero, b is a multiple of the gcd
   */

  resgcd = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) +
                                         (degb+1)*sizeof(int) );
  for (i = 0; i <= degb; i++)
  {
    resgcd->stem[i] = b[degb - i];
  }
  resgcd->stemdegree = degb;
  resgcd->minexpon = 0;
  free (a);
  free (b);
  free (c);
  return (resgcd);
}

/*
 *
 */

int
laurent_factor_content (struct laurentpoly *p)
{
  int i, cont;

  assert (p);
  cont = p->stem[0];
  for (i = 1; i <= p->stemdegree; i++) cont = gcd (cont, p->stem[i]);
  for (i = 0; i <= p->stemdegree; i++) p->stem[i] /= cont;

  return (cont);
}

/*
 * compute the GCD of two integers
 */

int
gcd (int a, int b)
{
  int asaved;

  b = abs(b);
  if (a == 0) return (b);
  a = abs(a);
  if (b == 0) return (a);

  if (a < b)
  {
    asaved = a;
    a = b;
    b = asaved;
  }

  return (gcd (b, a%b));
}
