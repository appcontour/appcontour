#ifndef FUNDAMENTAL_H_INCLUDED
  #include "fundamental.h"   // needed for structure definition
#endif

/*
 * prototypes
 */

/* laurent polynomials in one indeterminate (integral coefficient */

struct alexanderideal * groebner1 (struct alexanderideal *ai);
struct stemideal * ai2si (struct alexanderideal ai);
struct alexanderideal * si2ai (struct stemideal *si, struct alexanderideal *ai);
int groebner1_tryreduce (struct laurentpoly *l[], int gennum);
int groebner1_tryreduce_once (struct laurentpoly *l[], int gennum);
struct laurentpoly * groebner1_reduce_using_rule (struct laurentpoly *p,
                                                  struct laurentpoly *rule, int *retcode);
int groebner1_dropzeros (struct laurentpoly *l[], int gennum);
int groebner1_add_spolynomials (struct laurentpoly *l[], int gennum, int gendim);

/*
 * struct for working with stems in one indeterminate
 * TODO we should eventually be able to deal with integers of any size...
 */

typedef int Stemint;
#define STEMSIZE(n) (sizeof(struct stem) + (n)*sizeof(Stemint))

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
