#ifndef FUNDAMENTAL_H_INCLUDED
  #include "fundamental.h"   // needed for structure definition
#endif

/*
 * prototypes
 */

/* laurent polynomials in one indeterminate (integral coefficient */

struct alexanderideal * groebner1 (struct alexanderideal *ai);
int groebner1_tryreduce (struct laurentpoly *l[], int gennum);
int groebner1_tryreduce_once (struct laurentpoly *l[], int gennum);
struct laurentpoly * groebner1_reduce_using_rule (struct laurentpoly *p,
                                                  struct laurentpoly *rule, int *retcode);
int groebner1_dropzeros (struct laurentpoly *l[], int gennum);
int groebner1_add_spolynomials (struct laurentpoly *l[], int gennum, int gendim);
