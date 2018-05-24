#ifndef FUNDAMENTAL_H_INCLUDED
  #include "fundamental.h"   // needed for structure definition
#endif

typedef long long int Stemint;

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
int build_S_pols (struct stem *p1, struct stem *p2, struct stem **spoltoppt, struct stem **spolbottompt);
struct stem * reduce_pol_si (struct stem *p1, struct stemideal *si);
struct stem * reduce_pol_si_cycle (struct stem *p1, struct stemideal *si, int *statuspt);
struct stem * stem_normalize (struct stem *stem);
struct stem * stemideal_gcd (struct stemideal *si);
struct stem * stem_euclid (struct stem *p1, struct stem *p2);
Stemint stem_division (struct stem *dividend, struct stem *divisor, struct stem **quotientpt, struct stem **remainderpt);
void stem_canonify_sign (struct stem *p);
void printout_si (struct stemideal *si);
void printout_stem (struct stem *si);
void free_stemideal (struct stemideal *si);
Stemint gb_int_div (Stemint dividend, Stemint divisor);
Stemint stem_linf (struct stem *stem);
struct stem * stem_dup (struct stem *stem);
Stemint llgcd (Stemint a, Stemint b);
Stemint stem_factor_content (struct stem *stem);

/*
 * struct for working with stems in one indeterminate
 * TODO we should eventually be able to deal with integers of any size...
 */

#define STEMSIZE(n) (sizeof(struct stem) + (n)*sizeof(Stemint))
#define STEMIDEALSIZE(n) (sizeof(struct stemideal) + (n)*sizeof(struct stem *))
/* provide for extra space to avoid reallocation */
#define GB_EXTRAROOM(n) (1)
