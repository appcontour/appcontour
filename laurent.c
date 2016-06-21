/*
 * computations with laurent polynomials
 * in one or two indeterminates
 */

#include <assert.h>
#include <limits.h>
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
 * find the element of the ideal generated by p1 and p2 with smallest stemdegree and
 * smallest leading coefficient.
 * This should be the result of a kind of extended Euclidean algorithm, also giving
 * the coefficients of the linear combination.
 * The resulting polynomial should also be an integral multiple of the gcd of p1 and
 * p2.
 *
 * test e.g. with:
 * ./utils/kanenobu.sh -4 2 | ./contour --out alexander --foxd 2 --principal -v
 *
 * resulting in a 2x2 Alexander matrix and 3 should be the least generator of the ideal
 */

struct laurentpoly *
laurent_extended_euclid (struct laurentpoly *p1, struct laurentpoly *p2)
{
  struct laurentpoly *res;
  struct lpol *r_i, *r_im1, *s_i, *s_im1, *t_i, *t_im1, *temp;
  int i, j, k_i, ri0, rim10, c, alpha_i, beta_i;
  extern int verbose;

  if (verbose)
  {
    printf ("Extended Euclid between: "); print_laurentpoly (p1, "x"); printf (" --- ");
    print_laurentpoly (p2, "x"); printf ("\n");
  }
  if (p1 == 0) return (laurent_dup (p2));
  if (p2 == 0) return (laurent_dup (p1));

  bufsize = p1->stemdegree;
  if (p2->stemdegree > p1->stemdegree) bufsize = p2->stemdegree;
  if (bufsize == 0) /* special case of two nonzero monomials */
  {
    c = abs(gcd(p1->stem[0].l0, p2->stem[0].l0));
    res = (struct laurentpoly *) malloc (POLYSIZE(1));
    res->indets = 1;
    res->minexpon = 0;
    res->stemdegree = 0;
    res->stem[0].l0 = c;
    return (res);
  }
  bufsize += 2;
  bufsize *= 2;

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
    for (i = 0, j = r_im1->start; i <= p1->stemdegree; i++) r_im1->buf[j++] = p1->stem[i].l0;
    r_im1->end += p1->stemdegree + 1;
  }
  if (p2)
  {
    for (i = 0, j = r_i->start; i <= p2->stemdegree; i++) r_i->buf[j++] = p2->stem[i].l0;
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
  res = (struct laurentpoly *) malloc (POLYSIZE(r_i->end-r_i->start));
  res->indets = 1;
  res->minexpon = 0;
  res->stemdegree = r_i->end - r_i->start - 1;
  for (i = r_i->start, j = 0; i < r_i->end; i++) res->stem[j++].l0 = r_i->buf[i];
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
  int maxca = 0;
  int maxcb = 0;

  for (i = a->start; i < a->end; i++) if (maxca < abs(a->buf[i])) maxca = abs(a->buf[i]);
  for (i = b->start; i < b->end; i++) if (maxcb < abs(b->buf[i])) maxcb = abs(b->buf[i]);
  if (alpha) assert (maxca < (INT_MAX/abs(alpha)/2));
  if (beta) assert (maxcb < (INT_MAX/abs(beta)/2));
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
  if (term->indets <= 1)
  {
    for (i = 0; i <= term->stemdegree; i++) term->stem[i].l0 = -term->stem[i].l0;
    return;
  }
  for (i = 0; i <= term->stemdegree; i++) laurent_negate (term->stem[i].lx);

  return;
}

/*
 * multiply a poly by first indet
 */

void
laurent_mulu (struct laurentpoly *l)
{
  int i;

  if (l == 0) return;
  if (l->indets == 1)
  {
    l->minexpon++;
    return;
  }
  for (i = 0; i <= l->stemdegree; i++) laurent_mulu (l->stem[i].lx);
  return;
}

/*
 * multiply a poly by second indet
 */

void
laurent_mulv (struct laurentpoly *l)
{
  int i;

  if (l == 0) return;
  assert (l->indets >= 2);
  if (l->indets == 2)
  {
    l->minexpon++;
    return;
  }
  for (i = 0; i <= l->stemdegree; i++) laurent_mulv (l->stem[i].lx);
  return;
}

/*
 * print a laurent polynomial in a generic number of indeterminates
 */

void print_laurentpoly_map (struct laurentpoly *l, char *indetnames, char *mapappend);

void
print_laurentpoly (struct laurentpoly *l, char *indetnames)
{
  print_laurentpoly_map (l, indetnames, "");
}

void
print_laurentpoly_map (struct laurentpoly *l, char *indetnames, char *mapappend)
{
  struct laurentpoly *l1;
  int i, expon, indets;
  char appstring[80];

  if (l == 0) {printf ("0"); return;}
  indets = l->indets;
  assert (l->indets <= strlen (indetnames));

  if (indets <= 1)
  {
    assert (l->stem[0].l0);
    for (i = 0; i <= l->stemdegree; i++)
    {
      expon = i + l->minexpon;
      if (l->stem[i].l0)
      {
        if (abs(l->stem[i].l0) != 1 || expon == 0)
        {
          if (*mapappend == 0)
          {
            printf ("%+d", l->stem[i].l0);
          } else {
            if (abs(l->stem[i].l0) == 1)
            {
              printf ("%c", (l->stem[i].l0 > 0)?'+':'-');
            } else {
              printf ("%+d", l->stem[i].l0);
            }
            if (expon == 0) printf ("%s", mapappend);
          }
        } else {
          if (l->stem[i].l0 > 0) printf ("+");
           else printf ("-");
        }
        if (expon != 0)
        {
          printf ("%c", *indetnames);
          if (expon != 1)
          {
            if (expon > 0) printf ("^%d", expon);
              else printf ("^(%d)", expon);
          }
          if (*mapappend) printf ("%s", mapappend);
        }
      }
    }
    return;
  }

  assert (l->stem[0].lx);
  for (i = 0; i <= l->stemdegree; i++)
  {
    expon = i + l->minexpon;
    if ((l1 = l->stem[i].lx) != 0)
    {
      appstring[0] = 0;
      if (expon != 0)
      {
        appstring[0] = indetnames[indets - 1];
        appstring[1] = 0;
        if (expon != 1)
        {
          if (expon > 0) sprintf (appstring + 1, "^%d", expon);
            else sprintf (appstring + 1, "^(%d)", expon);
        }
      }
      strcat (appstring, mapappend);
      print_laurentpoly_map (l->stem[i].lx, indetnames, appstring);
    }
  }
}

/*
 * read a laurent polynomial in one indet
 */

struct laurentpoly *
read_laurentpoly1 (FILE *file, char indet_names[2])
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
    l1 = laurentpoly1_addmonom (l1, exp1, coef);
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

struct laurentpoly *
read_laurentpoly2 (FILE *file, char indet_names[2])
{
  char ch;
  int sign, coef, exp1, exp2;
  struct laurentpoly *l2 = 0;

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
 * add two laurent polynomials (k indeterminates)
 */

struct laurentpoly *
laurent_add (struct laurentpoly *a1, struct laurentpoly *a2)
{
  return (laurent_add_scal (a1, a2, 1));
}

/*
 * add two laurent polynomials (k indeterminates) with a2 multiplied by scal
 */

struct laurentpoly *
laurent_add_scal (struct laurentpoly *a1, struct laurentpoly *a2, int scal)
{
  int minexp, maxexp, k;
  struct laurentpoly *res;
  struct laurentpoly *addres;

  if (a2 == 0 || scal == 0) return (laurent_dup(a1));
  if (a1 == 0)
  {
    res = laurent_dup (a2);
    assert (abs (scal) == 1); /* for now */
    if (scal == -1) laurent_negate (res);
    return (res);
  }

  assert (a1->indets == a2->indets);
  minexp = a1->minexpon;
  if (a2->minexpon < minexp) minexp = a2->minexpon;

  maxexp = a1->minexpon + a1->stemdegree;
  if (a2->minexpon + a2->stemdegree > maxexp) maxexp = a2->minexpon + a2->stemdegree;

  res = (struct laurentpoly *) malloc (POLYSIZE(maxexp - minexp + 1));

  res->indets = a1->indets;
  res->minexpon = minexp;
  res->stemdegree = maxexp - minexp;

  if (res->indets == 1)
  {
    for (k = 0; k <= res->stemdegree; k++) res->stem[k].l0 = 0;
    for (k = 0; k <= a1->stemdegree; k++) res->stem[k + a1->minexpon - minexp].l0 = a1->stem[k].l0;
    for (k = 0; k <= a2->stemdegree; k++) res->stem[k + a2->minexpon - minexp].l0 += scal*a2->stem[k].l0;
  } else {
    for (k = 0; k <= res->stemdegree; k++) res->stem[k].lx = 0;
    for (k = 0; k <= a1->stemdegree; k++) res->stem[k + a1->minexpon - minexp].lx = laurent_dup (a1->stem[k].lx);
    for (k = 0; k <= a2->stemdegree; k++)
    {
      if (a2->stem[k].lx == 0) continue;
      addres = laurent_add_scal (res->stem[k + a2->minexpon - minexp].lx, a2->stem[k].lx, scal);
      if (res->stem[k + a2->minexpon - minexp].lx) free_laurentpoly (res->stem[k + a2->minexpon - minexp].lx);
      res->stem[k + a2->minexpon - minexp].lx = addres;
    }
  }

  /* normalize */

  if (res && laurent_normalize (res) == 0)
  {
    free_laurentpoly (res);
    res = 0;
  }

  return (res);
}

/*
 * multiply two laurent polynomials (in one indet)
 */

struct laurentpoly *
laurent_mul1 (struct laurentpoly *f1, struct laurentpoly *f2)
{
  int i, j;
  int resstemdegree;
  struct laurentpoly *res;

  if (f1 == 0 || f2 == 0) return (0);

  assert (f1->indets == 1 && f2->indets == 1);
  resstemdegree = f1->stemdegree + f2->stemdegree;

  res = (struct laurentpoly *) malloc (POLYSIZE(resstemdegree + 1));
  res->indets = 1;
  res->minexpon = f1->minexpon + f2->minexpon;
  res->stemdegree = resstemdegree;

  for (i = 0; i <= resstemdegree; i++) res->stem[i].l0 = 0;

  for (i = 0; i <= f1->stemdegree; i++)
  {
    for (j = 0; j <= f2->stemdegree; j++)
    {
      res->stem[i + j].l0 += f1->stem[i].l0*f2->stem[j].l0;
    }
  }

  return (res);
}

/*
 * multiply two laurent polynomials (k>=1 indeterminates)
 */

struct laurentpoly *
laurent_mul (struct laurentpoly *f1, struct laurentpoly *f2)
{
  int i, j;
  int resstemdegree;
  struct laurentpoly *res;
  struct laurentpoly *mulres, *addres;  /* these have one indet */

  if (f1 == 0 || f2 == 0) return (0);

  assert (f1->indets == f2->indets);
  if (f1->indets == 1) return (laurent_mul1 (f1, f2));

  assert (f1->indets >= 2);
  resstemdegree = f1->stemdegree + f2->stemdegree;

  res = (struct laurentpoly *) malloc (POLYSIZE(resstemdegree + 1));
  res->indets = f1->indets;
  res->minexpon = f1->minexpon + f2->minexpon;
  res->stemdegree = resstemdegree;

  for (i = 0; i <= resstemdegree; i++) res->stem[i].lx = 0;

  for (i = 0; i <= f1->stemdegree; i++)
  {
    for (j = 0; j <= f2->stemdegree; j++)
    {
      mulres = laurent_mul (f1->stem[i].lx, f2->stem[j].lx);
      addres = laurent_add (res->stem[i + j].lx, mulres);
      if (mulres) free_laurentpoly (mulres);
      if (res->stem[i + j].lx) free_laurentpoly (res->stem[i + j].lx);
      res->stem[i + j].lx = addres;
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

  if (l->indets > 1)
  {
    while (l->stem[0].lx == 0)
    {
      if (l->stemdegree == 0)
      {
        return (0);
      }
      for (k = 0; k < l->stemdegree; k++) l->stem[k].lx = l->stem[k+1].lx;
      l->stemdegree--;
      l->minexpon++;
    }

    while (l->stem[l->stemdegree].lx == 0)
    {
      assert (l->stemdegree > 0);
      l->stemdegree--;
    }
  } else {
    while (l->stem[0].l0 == 0)
    {
      if (l->stemdegree == 0)
      {
        return (0);
      }
      for (k = 0; k < l->stemdegree; k++) l->stem[k].l0 = l->stem[k+1].l0;
      l->stemdegree--;
      l->minexpon++;
    }

    while (l->stem[l->stemdegree].l0 == 0)
    {
      assert (l->stemdegree > 0);
      l->stemdegree--;
    }
  }
  return (l);
}

/*
 * duplicate a laurent polynomial (k indeterminates)
 */

struct laurentpoly *
laurent_dup (struct laurentpoly *l)
{
  int k;
  struct laurentpoly *res;

  if (l == 0) return (0);

  res = (struct laurentpoly *) malloc (POLYSIZE(l->stemdegree + 1));

  res->indets = l->indets;
  res->minexpon = l->minexpon;
  res->stemdegree = l->stemdegree;

  if (l->indets > 1)
  {
    for (k = 0; k <= l->stemdegree; k++) res->stem[k].lx = laurent_dup (l->stem[k].lx);
  } else {
    for (k = 0; k <= l->stemdegree; k++) res->stem[k].l0 = l->stem[k].l0;
  }

  return (res);
}

/*
 * canonify with respect to t->1/t
 *
 * note that the sign is such that evaluation in 1 gives
 * 1 (should be careful in general)
 */

void
laurent_canonify1 (struct laurentpoly *l)
{
  int k, kk, sign;

  if (l == 0) return;
  assert (l->indets == 1);

  assert (l->stem[0].l0);
  assert (l->stem[l->stemdegree].l0);
  //if (l->stem[0] < 0) laurent_negate (l);

  sign = 1;
  if (l->stem[l->stemdegree].l0 < 0) sign = -1;

  /* Making lexicographic comparison */

  for (k = 0, kk = l->stemdegree; k < kk; k++, kk--)
  {
    if (l->stem[k].l0 > sign*l->stem[kk].l0)
    {
      laurent_t_to_oneovert (l);
      if (l->stem[0].l0 < 0) laurent_negate (l);
      l->minexpon = 0;
      return;
    }
  }
  l->minexpon = 0; // we are interested only to the stem
  return;
}

/*
 * laurent_canonifyx
 * for now simply shift exponents to minimal nonnegative
 */

void
laurent_canonifyx (struct laurentpoly *l)
{
  int i, minexpon1;
  struct laurentpoly *l1;

  if (l == 0) return;

  assert (l->indets == 2);
  assert (l->stem[0].lx);

  l->minexpon = 0;

  minexpon1 = (l->stem[0].lx)->minexpon;
  for (i = 0; i <= l->stemdegree; i++)
  {
    l1 = l->stem[i].lx;
    if (l1 == 0) continue;
    if (l1->minexpon < minexpon1) minexpon1 = l1->minexpon;
  }
  for (i = 0; i <= l->stemdegree; i++)
  {
    l1 = l->stem[i].lx;
    if (l1 == 0) continue;
    l1->minexpon -= minexpon1;
  }
}

/*
 * transform t->1/t (one indet laurentpolyx)
 */

void
laurent_t_to_oneovert (struct laurentpoly *l)
{
  int k, kk, saved;

  if (l == 0) return;
  assert (l->indets == 1);

  l->minexpon = -(l->minexpon + l->stemdegree);

  for (k = 0, kk = l->stemdegree; k < kk; k++, kk--)
  {
    saved = l->stem[k].l0;
    l->stem[k].l0 = l->stem[kk].l0;
    l->stem[kk].l0 = saved;
  }
  return;
}

/*
 * free data associated to a laurent polynomial with 2 indeterminates
 */

void
free_laurentpoly (struct laurentpoly *l)
{
  int i;

  if (l == 0) return;

  if (l->indets > 1)
  {
    for (i = 0; i <= l->stemdegree; i++)
    {
      if (l->stem[i].lx) free_laurentpoly (l->stem[i].lx);
    }
  }
  free (l);
}

/*
 * add a monomial to a laurent_poly in one indet
 */

struct laurentpoly *
laurentpoly1_addmonom (struct laurentpoly *l, int expon, int coef)
{
  struct laurentpoly *m = 0;
  struct laurentpoly *res;

  if (coef == 0) return (l);

  if (l) assert (l->indets == 1);
  if (l && expon >= l->minexpon && expon <= l->minexpon + l->stemdegree)
  {
    l->stem[expon - l->minexpon].l0 += coef;
    if (laurent_normalize(l) == 0)
    {
      free_laurentpoly (l);
      return (0);
    }
    return (l);
  }

  m = (struct laurentpoly *) malloc (POLYSIZE(1));
  m->indets = 1;
  m->minexpon = expon;
  m->stemdegree = 0;
  m->stem[0].l0 = coef;

  if (l == 0) return (m);

  res = laurent_add (l, m);
  free_laurentpoly (m);
  free_laurentpoly (l);

  return (res);
}

/*
 * add a monomial to a laurentpolyx (two indeterminates)
 */

struct laurentpoly *
laurentpoly2_addmonom (struct laurentpoly *l, int degu, int degv, int coef)
{
  int exponvec[2];

  if (coef == 0) return (l);

  if (l) assert (l->indets == 2);

  exponvec[0] = degu;
  exponvec[1] = degv;
  return (laurentpoly_addmonom (l, 2, exponvec, coef));
}

/*
 * add monomial to a generic laurent polinomial
 */

struct laurentpoly *
laurentpoly_addmonom (struct laurentpoly *l, int indets, int *exponvec, int scal)
{
  struct laurentpoly *m, *res;
  int iv, expon;

  if (scal == 0) return (l);
  expon = exponvec[indets-1];

  if (l && expon >= l->minexpon && expon <= l->minexpon + l->stemdegree)
  {
    if (indets <= 1)
    {
      l->stem[expon - l->minexpon].l0 += scal;
    } else {
      iv = expon - l->minexpon;
      l->stem[iv].lx = laurentpoly_addmonom (l->stem[iv].lx, indets - 1, exponvec, scal);
    }
    if (laurent_normalize(l) == 0)
    {
      free_laurentpoly (l);
      return (0);
    }
    return (l);
  }

  m = (struct laurentpoly *) malloc (POLYSIZE(1));
  m->indets = indets;
  m->minexpon = expon;
  m->stemdegree = 0;
  if (indets <= 1)
  {
    m->stem[0].l0 = scal;
  } else {
    m->stem[0].lx = laurentpoly_addmonom (0, indets-1, exponvec, scal);
  }

  if (l == 0) return (m);

  res = laurent_add (l, m);
  free_laurentpoly (m);
  free_laurentpoly (l);

  return (res);
}

/*
 * add up all coefficients of a laurentpoly in one indet
 * (equivalent to evaluation in t=1)
 */

int
laurent_sum_coefficients1 (struct laurentpoly *l)
{
  int k;
  int res = 0;

  if (l == 0) return (0);
  assert (l->indets == 1);

  for (k = 0; k <= l->stemdegree; k++) res += l->stem[k].l0;
  return (res);
}

/*
 * add up all coefficients of a laurentpolyx in two indets
 * (equivalent to evaluation in v=1)
 */

struct laurentpoly *
laurent_sum_coefficients2 (struct laurentpoly *l)
{
  int k;
  struct laurentpoly *res, *addres;

  if (l == 0) return (0);
  assert (l->indets == 2);
  res = laurent_dup(l->stem[0].lx);

  for (k = 1; k <= l->stemdegree; k++)
  {
    addres = laurent_add (res, l->stem[k].lx);
    if (res) free_laurentpoly (res);
    res = addres;
  }
  return (res);
}

/*
 * add up all coefficients of each coefficient or a laurentpolyx in two indets
 * (equivalent to evaluation in u=1)
 */

struct laurentpoly *
laurent_sum_each_coefficient2 (struct laurentpoly *l)
{
  int k;
  struct laurentpoly *res;

  if (l == 0) return (0);
  assert (l->indets == 2);

  res = (struct laurentpoly *) malloc (POLYSIZE(l->stemdegree + 1));
  res->indets = 1;
  res->minexpon = l->minexpon;
  res->stemdegree = l->stemdegree;

  for (k = 0; k <= l->stemdegree; k++)
  {
    res->stem[k].l0 = laurent_sum_coefficients1 (l->stem[k].lx);
  }
  if (laurent_normalize (res) == 0)
  {
    free_laurentpoly (res);
    return (0);
  }

  return (res);
}

/*
 * get total degre (after shift) of a laurentpolyx (two indets)
 */

int
laurentx_totdegree (struct laurentpoly *l)
{
  struct laurentpoly *lu;  /* one indet poly */
  int minexpv, minexpu, maxtotdegree, totdegree;
  int iv;

  if (l == 0) return (-1);
  assert (l->indets == 2);
  assert (l->stem[0].lx);
  minexpv = l->minexpon;
  lu = l->stem[0].lx;
  minexpu = lu->minexpon;
  assert (lu->stem[lu->stemdegree].l0);
  maxtotdegree = minexpv + minexpu + lu->stemdegree;
  for (iv = 1; iv <= l->stemdegree; iv++)
  {
    lu = l->stem[iv].lx;
    if (lu == 0) continue;
    if (lu->minexpon < minexpu) minexpu = lu->minexpon;
    assert (lu->stem[lu->stemdegree].l0);
    totdegree = minexpv + iv + lu->minexpon + lu->stemdegree;
    if (totdegree > maxtotdegree) maxtotdegree = totdegree;
  }
  return (maxtotdegree - minexpu - minexpv);
}

/*
 * lexicographic comparison between two laurent polynomials in two indets
 * assumptions:
 * 1. They are canonified by monomial multiplication
 * 2. They have the same total degree
 * comparison changes sign of poly if first monomial is negative
 */

int
laurentx_lexicocompare (struct laurentpoly *p1, struct laurentpoly *p2)
{
  int sign1 = 1, sign2 = 1;
  int k, ku;
  struct laurentpoly *pu1, *pu2;  /* one indet */

  if (p1 == 0 && p2 == 0) return (0);
  if (p1 == 0) return (-1);
  if (p2 == 0) return (1);

  assert (p1->indets == 2 && p2->indets == 2);
  if ((p1->stem[0].lx)->stem[0].l0 < 0) sign1 = -1;
  if ((p2->stem[0].lx)->stem[0].l0 < 0) sign2 = -1;

  /* first test: degree with respect to v */

  if (p1->stemdegree < p2->stemdegree) return (-1);
  if (p1->stemdegree > p2->stemdegree) return (1);

  /* second test: lexicographic degree in u */

  for (k = p1->stemdegree; k >= 0; k--)
  {
    pu1 = p1->stem[k].lx;
    pu2 = p2->stem[k].lx;
    if (pu1 == 0 && pu2 == 0) continue;
    if (pu1 == 0) return (-1);
    if (pu2 == 0) return (1);
    if (pu1->minexpon + pu1->stemdegree < pu2->minexpon + pu2->stemdegree) return (-1);
    if (pu1->minexpon + pu1->stemdegree > pu2->minexpon + pu2->stemdegree) return (1);
  }

  /* third test: lexicographic minexpon in u */

  for (k = 0; k <= p1->stemdegree; k++)
  {
    pu1 = p1->stem[k].lx;
    pu2 = p2->stem[k].lx;
    if (pu1 == 0) continue;
    if (pu1->minexpon > pu2->minexpon) return (-1);
    if (pu1->minexpon < pu2->minexpon) return (1);
  }

  /* finally full lexicographic comparison */

  for (k = 0; k <= p1->stemdegree; k++)
  {
    pu1 = p1->stem[k].lx;
    pu2 = p2->stem[k].lx;

    if (pu1 == 0) continue;
    for (ku = 0; ku <= pu1->stemdegree; ku++)
    {
      if (abs(pu1->stem[ku].l0) < abs(pu2->stem[ku].l0)) return (-1);
      if (abs(pu1->stem[ku].l0) > abs(pu2->stem[ku].l0)) return (1);
      if (pu1->stem[ku].l0 == 0) continue;
      if (sign1*pu1->stem[ku].l0 > sign2*pu2->stem[ku].l0) return (-1);
      if (sign1*pu1->stem[ku].l0 < sign2*pu2->stem[ku].l0) return (1);
    }
  }
  return (0);
}

/*
 * canonify sign such that either p(1) > 0
 * or the first nonzero coefficient is positive (one indet)
 */

void
laurent_canonifysign1 (struct laurentpoly *p)
{
  int val1 = 0;
  int k;

  if (p == 0) return;
  assert (p->indets == 1);

  /* evaluate in (1) */
  for (k = 0; k <= p->stemdegree; k++)
    val1 += p->stem[k].l0;

  if (val1 > 0) return;

  if (val1 == 0)
  {
    assert (p->stem[0].l0 != 0);
    if (p->stem[0].l0 > 0) return;
  }

  /* if here change sign */
  for (k = 0; k <= p->stemdegree; k++)
    p->stem[k].l0 = - p->stem[k].l0;
}

/*
 * canonify sign such that either p(1,1) > 0
 * or the first nonzero coefficient is positive
 */

void
laurent_canonifysign2 (struct laurentpoly *p)
{
  int val11 = 0;
  int kv, ku;
  struct laurentpoly *pu; /* one indet */

  if (p == 0) return;
  assert (p->indets == 2);

  /* evaluate in (1,1) */
  for (kv = 0; kv <= p->stemdegree; kv++)
  {
    pu = p->stem[kv].lx;
    if (pu == 0) continue;
    for (ku = 0; ku <= pu->stemdegree; ku++) val11 += pu->stem[ku].l0;
  }

  if (val11 > 0) return;

  if (val11 == 0)
  {
    assert ((p->stem[0].lx)->stem[0].l0 != 0);
    if ((p->stem[0].lx)->stem[0].l0 > 0) return;
  }

  /* if here change sign */
  for (kv = 0; kv <= p->stemdegree; kv++)
  {
    pu = p->stem[kv].lx;
    if (pu == 0) continue;
    for (ku = 0; ku <= pu->stemdegree; ku++) pu->stem[ku].l0 = -pu->stem[ku].l0;
  }
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
  if (p1) {for (i = 0; i <= p1->stemdegree; i++) p1->stem[i].l0 *= inspread*cont_p1;}  /* restore p1 */
  resgcd = laurent_extended_euclid (p1, p2);
  cont_p2 = laurent_factor_content (p2);

  alpha = laurent_factor_content (resgcd);

  if (internalcheck)
  {
    resgcd2 = laurent_euclid (p1, p2);
    laurent_factor_content (resgcd2);
    assert (resgcd2->stemdegree == resgcd->stemdegree);
    for (i = 0; i <= resgcd2->stemdegree; i++) assert (resgcd2->stem[i].l0 == resgcd->stem[i].l0);
    if (resgcd2) free (resgcd2);
  }
  cont_gcd = gcd (cont_p1, cont_p2);
  if (cont_gcd < 0) cont_gcd = -cont_gcd;
  if (alpha < 0) alpha = -alpha;
  *spreadpt = alpha/cont_gcd;
  assert (alpha == (*spreadpt)*cont_gcd);
  for (i = 0; i <= p1->stemdegree; i++) p1->stem[i].l0 /= inspread;
  for (i = 0; i <= p2->stemdegree; i++) p2->stem[i].l0 *= cont_p2;
  for (i = 0; i <= resgcd->stemdegree; i++) resgcd->stem[i].l0 *= cont_gcd;  /* obtain real gcd */

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
  for (i = 0; i <= dega; i++) a[i] = p1->stem[dega - i].l0;
  for (i = 0; i <= degb; i++) b[i] = p2->stem[degb - i].l0;

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

  resgcd = (struct laurentpoly *) malloc (POLYSIZE(degb+1));
  resgcd->indets = 1;
  for (i = 0; i <= degb; i++)
  {
    resgcd->stem[i].l0 = b[degb - i];
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
  assert (p->stem[0].l0 != 0);
  sign = 1;
  if (p->stem[0].l0 < 0) sign = -1;
  cont = p->stem[0].l0;
  for (i = 1; i <= p->stemdegree; i++) cont = gcd (cont, p->stem[i].l0);
  if (cont*sign < 0) cont = -cont;
  for (i = 0; i <= p->stemdegree; i++) p->stem[i].l0 /= cont;

  return (cont);
}

/*
 * divide a polynomial in k (2 or more) indeterminates by (1-w) where w indicates
 * the last unknown.
 * return the quotient (in k indets) and the rest (in k-1 indets)
 */

struct laurentpoly *
laurent_divide_by_1minusw (struct laurentpoly *l, struct laurentpoly **rpt)
{
  struct laurentpoly *q;
  struct laurentpoly *ltemp;
  struct laurentpoly *rest;  /* one less indeterminate */
  int i;

  *rpt = 0;
  if (l == 0) return (0);

  assert (l->indets > 1);  /* we shall call this with 3 indets, but it works anyway */

  q = (struct laurentpoly *) malloc (POLYSIZE (l->stemdegree));
  q->indets = l->indets;
  q->minexpon = l->minexpon;
  q->stemdegree = l->stemdegree - 1;

  for (i = q->stemdegree; i >= 0; i--)
  {
    q->stem[i].lx = laurent_dup (l->stem[i + 1].lx);
    laurent_negate (q->stem[i].lx);
    if (i < q->stemdegree)
    {
      ltemp = laurent_add (q->stem[i].lx, q->stem[i+1].lx);
      if (q->stem[i].lx) free_laurentpoly (q->stem[i].lx);
      q->stem[i].lx = ltemp;
    }
  }
  if (q->stem[0].lx == 0)
  {
    rest = laurent_dup (l->stem[0].lx);
  } else {
    rest = laurent_add_scal (l->stem[0].lx, q->stem[0].lx, -1);
  }
  if (q)
  {
    ltemp = laurent_normalize (q);
    if (ltemp == 0) free_laurentpoly (q);
    q = ltemp;
  }
  *rpt = rest;
  return (q);
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
