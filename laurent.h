/*
 * definitions for alexander polynomial
 */

struct laurentpoly {
  int minexpon;
  int stemdegree;
  int stem[];
};

struct laurentpoly2 {
  int minexpon;
  int stemdegree;
  struct laurentpoly *stem[];
};

void laurent_negate (struct laurentpoly *term);
void laurent_negate2 (struct laurentpoly2 *term);
void print_laurentpoly (struct laurentpoly *l, char indet);
void print_laurentpoly2 (struct laurentpoly2 *l, char indet1, char indet2);
struct laurentpoly *read_laurentpoly (FILE *file, char indet_names[2]);
struct laurentpoly2 *read_laurentpoly2 (FILE *file, char indet_names[2]);
struct laurentpoly *laurent_add (struct laurentpoly *add1, struct laurentpoly *add2);
struct laurentpoly2 *laurent_add2 (struct laurentpoly2 *add1, struct laurentpoly2 *add2);
struct laurentpoly *laurent_mul (struct laurentpoly *fact1, struct laurentpoly *fact2);
struct laurentpoly2 *laurent_mul2 (struct laurentpoly2 *fact1, struct laurentpoly2 *fact2);
struct laurentpoly *laurent_normalize (struct laurentpoly *l);
struct laurentpoly2 *laurent_normalize2 (struct laurentpoly2 *l);
struct laurentpoly *laurent_dup (struct laurentpoly *l);
struct laurentpoly2 *laurent_dup2 (struct laurentpoly2 *l);
void laurent_canonify (struct laurentpoly *l);
void laurent_canonify2 (struct laurentpoly2 *l);
void laurent_t_to_oneovert (struct laurentpoly *l);
void free_laurentpoly2 (struct laurentpoly2 *l);
struct laurentpoly *laurentpoly_addmonom (struct laurentpoly *l, int expon, int coef);
struct laurentpoly2 *laurentpoly2_addmonom (struct laurentpoly2 *l, int expu, int expv, int coef);
int laurent_sum_coefficients (struct laurentpoly *l);
struct laurentpoly *laurent_sum_coefficients2 (struct laurentpoly2 *l);
struct laurentpoly *laurent_sum_each_coefficient2 (struct laurentpoly2 *l);

struct laurentpoly *laurent_gcd (int inspread, struct laurentpoly *p1, struct laurentpoly *p2, int *spreadpt);
int laurent_factor_content (struct laurentpoly *p);
struct laurentpoly *laurent_euclid (struct laurentpoly *p1, struct laurentpoly *p2);
struct laurentpoly *laurent_extended_euclid (struct laurentpoly *p1, struct laurentpoly *p2);

int gcd (int a, int b);
