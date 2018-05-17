#ifndef FUNDAMENTAL_H_INCLUDED
  #include "fundamental.h"   // needed for structure definition
#endif

typedef int Stemint;

struct stem {
  int dim;
  int degree;
  Stemint coef[];
};

struct stemideal {
  int dim;
  int num;
  struct stem *stem[];
};

/*
 * prototypes
 */

/* laurent polynomials in one indeterminate (integral coefficient */

struct alexanderideal * groebner1 (struct alexanderideal *ai);
struct stemideal * ai2si (struct alexanderideal *ai);
struct stem * lp2stem (struct laurentpoly *lp);
struct alexanderideal * si2ai (struct stemideal *si, struct alexanderideal *ai);
struct laurentpoly * stem2lp (struct stem *stem);
int groebner1_tryreduce (struct stemideal *si);
int groebner1_tryreduce_once (struct stemideal *si);
struct stem * groebner1_reduce_using_rule (struct stem *p,
                                                  struct stem *rule, int *retcode);
int groebner1_dropzeros (struct stemideal *si);
int groebner1_add_spolynomials (struct stemideal *si);
void printout_si (struct stemideal *si);
void printout_stem (struct stem *si);
void free_stemideal (struct stemideal *si);
int gb_int_div (int dividend, int divisor);

/*
 * struct for working with stems in one indeterminate
 * TODO we should eventually be able to deal with integers of any size...
 */

#define STEMSIZE(n) (sizeof(struct stem) + (n)*sizeof(Stemint))
#define STEMIDEALSIZE(n) (sizeof(struct stemideal) + (n)*sizeof(struct stem *))
/* provide for extra space to avoid reallocation */
#define GB_EXTRAROOM(n) (0)
