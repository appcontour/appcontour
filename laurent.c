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
 * implement extended Euclidean algorithm for two Laurent polynomials
 */

struct lpol {
    int start;
    int end;
    int buf[];
  };

void euclid_alfabetatk (struct lpol *a, struct lpol *b, int alfa, int beta, int k);
void euclid_reduce3 (struct lpol *a, struct lpol *b, struct lpol *c);
void euclid_printlpol (struct lpol *a);

static int bufsize;

/*
 * test fatti con:
 * ./utils/kanenobu.sh -4 2 | ./contour --out alexander --foxd 2 --principal -v
 *
 * che risulta in una matrice 2x2, e un generatore dell'ideale secondo e' 3
 */

struct laurentpoly *
laurent_extended_euclid (struct laurentpoly *p1, struct laurentpoly *p2)
{
  struct laurentpoly *res;
  struct lpol *r_i, *r_im1, *s_i, *s_im1, *t_i, *t_im1, *temp;
  int i, j, k_i, ri0, rim10, c, alpha_i, beta_i;
  extern int verbose;

  bufsize=100;  //TODO  dare un valore decente!!
  /*
   *
   */

  if (verbose)
  {
    printf ("===> Extended Euclid: "); print_laurentpoly (p1, 'x'); printf (" --- ");
    print_laurentpoly (p2, 'x'); printf ("\n");
  }
  if (p1 == 0) return (laurent_dup (p2));
  if (p2 == 0) return (laurent_dup (p1));

  r_im1 = (struct lpol *) malloc (sizeof (struct lpol) + bufsize*sizeof (int));
  r_i = (struct lpol *) malloc (sizeof (struct lpol) + bufsize*sizeof (int));
  s_im1 = (struct lpol *) malloc (sizeof (struct lpol) + bufsize*sizeof (int));
  s_i = (struct lpol *) malloc (sizeof (struct lpol) + bufsize*sizeof (int));
  t_im1 = (struct lpol *) malloc (sizeof (struct lpol) + bufsize*sizeof (int));
  t_i = (struct lpol *) malloc (sizeof (struct lpol) + bufsize*sizeof (int));

  for (i = 0; i < bufsize; i++)
  {
    r_im1->buf[i] = 0;
    s_im1->buf[i] = 0;
    t_im1->buf[i] = 0;
    r_i->buf[i] = 0;
    s_i->buf[i] = 0;
    t_i->buf[i] = 0;
  }

  r_im1->start = r_i->start = r_im1->end = r_i->end = bufsize/2;
  s_im1->start = s_i->start = s_im1->end = s_i->end = bufsize/2;
  t_im1->start = t_i->start = t_im1->end = t_i->end = bufsize/2;
  if (p1)
  {
    for (i = 0, j = r_im1->start; i <= p1->stemdegree; i++) r_im1->buf[j++] = p1->stem[i];
    r_im1->end += p1->stemdegree + 1;
  }
  if (p2)
  {
    for (i = 0, j = r_i->start; i <= p2->stemdegree; i++) r_i->buf[j++] = p2->stem[i];
    r_i->end += p2->stemdegree + 1;
  }
  s_im1->buf[s_im1->start] = 1;
  s_im1->end = s_im1->start + 1;
  t_im1->end = t_im1->start;
  s_i->end = s_i->start;
  t_i->buf[t_i->start] = 1;
  t_i->end = t_i->start + 1;

  while (1)
  {
    if ((r_im1->end - r_im1->start) < (r_i->end - r_i->start))
    {
      temp = r_im1; r_im1 = r_i; r_i = temp;
      temp = s_im1; s_im1 = s_i; s_i = temp;
      temp = t_im1; t_im1 = t_i; t_i = temp;
    }
    /* first: find alfa_i and q_i=beta_i t^k_i */
    //r_i_size = r_i->end - e_i->start;
    //if (r_i_size < 0) r_i_size += bufsize;
    k_i = r_im1->start - r_i->start;
    ri0 = r_i->buf[r_i->start];
    rim10 = r_im1->buf[r_im1->start];
    c = gcd (ri0, rim10);
    alpha_i = ri0/c;
    beta_i = rim10/c;
    euclid_alfabetatk (r_im1, r_i, alpha_i, beta_i, k_i);
    euclid_alfabetatk (s_im1, s_i, alpha_i, beta_i, k_i);
    euclid_alfabetatk (t_im1, t_i, alpha_i, beta_i, k_i);
    euclid_reduce3 (r_im1, s_im1, t_im1);

    if (r_im1->start == r_im1->end) break;
  }

  if (verbose)
  {
    printf ("gcd: "); euclid_printlpol (r_i);
    printf ("  s: "); euclid_printlpol (s_i);
    printf ("  t: "); euclid_printlpol (t_i);
  }
  assert (r_i->start != r_i->end);
  res = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + (r_i->end-r_i->start)*sizeof(int));
  res->minexpon = 0;
  res->stemdegree = r_i->end - r_i->start - 1;
  for (i = r_i->start, j = 0; i < r_i->end; i++) res->stem[j++] = r_i->buf[i];
  free (r_im1);
  free (s_im1);
  free (t_im1);
  free (r_i);
  free (s_i);
  free (t_i);
  return (res);
}

/*
 *
 */

void
euclid_printlpol (struct lpol *a)
{
  int i;

  printf ("start: %d,", a->start);
  for (i = a->start; i < a->end; i++) printf (" %d", a->buf[i]);
  printf ("\n");
}

void
euclid_alfabetatk (struct lpol *a, struct lpol *b, int alpha, int beta, int k)
{
  int i;

  assert (b->start + k >= 0);
  assert (b->end + k < bufsize);
  for (i = a->start; i < a->end; i++) a->buf[i] *= alpha;
  for (i = b->start; i < b->end; i++) a->buf[i+k] -= beta*b->buf[i];
  if (a->end < b->end + k) a->end = b->end + k;
  if (a->start > b->start + k) a->start = b->start + k;
  while (a->buf[a->start] == 0 && a->start < a->end) a->start++;
  if (a->start == a->end) return;
  while (a->buf[a->end - 1] == 0) a->end--;
}

void
euclid_reduce3 (struct lpol *a, struct lpol *b, struct lpol *c)
{
  int i;
  int cgcd = 0;

  for (i = a->start; i < a->end; i++) cgcd = gcd (cgcd, a->buf[i]);
  for (i = b->start; i < b->end; i++) cgcd = gcd (cgcd, b->buf[i]);
  for (i = c->start; i < c->end; i++) cgcd = gcd (cgcd, c->buf[i]);

  if (cgcd == 0) return;
  if (cgcd < 0) cgcd = -cgcd;
  if (cgcd == 1) return;

  for (i = a->start; i < a->end; i++) a->buf[i] /= cgcd;
  for (i = b->start; i < b->end; i++) b->buf[i] /= cgcd;
  for (i = c->start; i < c->end; i++) c->buf[i] /= cgcd;
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
 * if unsure, call with "inspread = 1"
 */

struct laurentpoly *
laurent_gcd (int inspread, struct laurentpoly *p1, struct laurentpoly *p2, int *spreadpt)
{
  struct laurentpoly *resgcd;
  struct laurentpoly *resgcd2;
  int i, alpha, cont_p1, cont_p2, cont_gcd;
  extern int verbose, internalcheck;

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
  if (p1) {for (i = 0; i <= p1->stemdegree; i++) p1->stem[i] *= inspread*cont_p1;}  /* restore p1 */
  resgcd = laurent_extended_euclid (p1, p2);
  cont_p2 = laurent_factor_content (p2);

  alpha = laurent_factor_content (resgcd);

  if (internalcheck)
  {
    resgcd2 = laurent_euclid (p1, p2);
    laurent_factor_content (resgcd2);
    assert (resgcd2->stemdegree == resgcd->stemdegree);
    for (i = 0; i <= resgcd2->stemdegree; i++) assert (resgcd2->stem[i] == resgcd->stem[i]);
    if (resgcd2) free (resgcd2);
  }
  cont_gcd = gcd (cont_p1, cont_p2);
  if (cont_gcd < 0) cont_gcd = -cont_gcd;
  if (alpha < 0) alpha = -alpha;
  *spreadpt = alpha/cont_gcd;
  assert (alpha == (*spreadpt)*cont_gcd);
  for (i = 0; i <= p1->stemdegree; i++) p1->stem[i] /= inspread;
  for (i = 0; i <= p2->stemdegree; i++) p2->stem[i] *= cont_p2;
  for (i = 0; i <= resgcd->stemdegree; i++) resgcd->stem[i] *= cont_gcd;  /* obtain real gcd */

  if (verbose) printf ("spread: %d\n", *spreadpt);
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
  int i, sign, cont;

  assert (p);
  assert (p->stem[0] != 0);
  sign = 1;
  if (p->stem[0] < 0) sign = -1;
  cont = p->stem[0];
  for (i = 1; i <= p->stemdegree; i++) cont = gcd (cont, p->stem[i]);
  if (cont*sign < 0) cont = -cont;
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
