/*
 * definitions for alexander polynomial
 */

union intorpointer {
  struct laurentpoly *lx;
  int                  l0;
};

struct laurentpoly {
  int indets;
  int minexpon;
  int stemdegree;
  union intorpointer stem[];
};

#define POLYSIZE(n) (sizeof(struct laurentpoly) + (n)*sizeof(union intorpointer))

struct laurentpoly *laurent_buildconstant (int indets, int cnst);
void laurent_negate (struct laurentpoly *term);
void laurent_mulscal (struct laurentpoly *term, int scal);
void laurent_mulu (struct laurentpoly *l);
void laurent_mulv (struct laurentpoly *l);
void print_laurentpoly (struct laurentpoly *l, char *indetlist);
struct laurentpoly *read_laurentpoly1 (FILE *file, char indet_names[2]);
struct laurentpoly *read_laurentpoly2 (FILE *file, char indet_names[2]);
struct laurentpoly *laurent_add (struct laurentpoly *add1, struct laurentpoly *add2);
struct laurentpoly *laurent_add_scal (struct laurentpoly *add1, struct laurentpoly *add2, int scal);
struct laurentpoly *laurent_addto_scal (struct laurentpoly *add1, struct laurentpoly *add2, int scal);
struct laurentpoly *laurent_mul1 (struct laurentpoly *fact1, struct laurentpoly *fact2);
struct laurentpoly *laurent_mul (struct laurentpoly *fact1, struct laurentpoly *fact2);
struct laurentpoly *laurent_normalize (struct laurentpoly *l);
struct laurentpoly *laurent_dup (struct laurentpoly *l);
void laurent_canonify1 (struct laurentpoly *l);
void laurent_canonify_exponents (struct laurentpoly *l);
void laurent_canonify_exponent (struct laurentpoly *l, int indet);
int laurent_compute_minexpon (struct laurentpoly *l, int indet);
void laurent_mulpowindet (struct laurentpoly *l, int indet, int expon);
void laurent_t_to_oneovert (struct laurentpoly *l);
void free_laurentpoly (struct laurentpoly *l);
struct laurentpoly *laurentpoly_addmonom (struct laurentpoly *l, int rank, int *expuv, int coef);
struct laurentpoly *laurentpoly1_addmonom (struct laurentpoly *l, int expon, int coef);
struct laurentpoly *laurentpoly2_addmonom (struct laurentpoly *l, int expu, int expv, int coef);
int laurent_sum_coefficients1 (struct laurentpoly *l);
struct laurentpoly *laurent_sum_coefficients2 (struct laurentpoly *l);
struct laurentpoly *laurent_sum_each_coefficient2 (struct laurentpoly *l);

struct laurentpoly *laurent_gcd (int inspread, struct laurentpoly *p1, struct laurentpoly *p2, int *spreadpt);

int laurent_factor_content (struct laurentpoly *p);
struct laurentpoly *laurent_euclid (struct laurentpoly *p1, struct laurentpoly *p2);
struct laurentpoly *laurent_extended_euclid (struct laurentpoly *p1, struct laurentpoly *p2);
int laurentx_totdegree (struct laurentpoly *l);
int laurentx_lexicocompare (struct laurentpoly *p1, struct laurentpoly *p2);
void laurent_canonifysign1 (struct laurentpoly *p);
void laurent_canonifysign2 (struct laurentpoly *p);

struct laurentpoly *laurent_divide_by_1minusw (struct laurentpoly *l, struct laurentpoly **rpt);

int gcd (int a, int b);
int binomial (int n, int k);
