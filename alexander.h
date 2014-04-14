/*
 * definitions for alexander polynomial
 */

struct laurentpoly {
  int minexpon;
  int stemdegree;
  int denom;
  int stem[];
};

/* prototypes */

int alexander (struct presentation *p);
int corank_one_alexander (struct presentation *p);
struct laurentpoly *laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate);
struct laurentpoly *laurent_get_exp_sum (struct presentationrule *r, int g, int gconj);
struct laurentpoly *laurent_compute_determinant (struct laurentpoly ***matrix, int n);
void print_laurentpoly (struct laurentpoly *l);
struct laurentpoly *laurent_add (struct laurentpoly *add1, struct laurentpoly *add2);
struct laurentpoly *laurent_mul (struct laurentpoly *fact1, struct laurentpoly *fact2);
void laurent_negate (struct laurentpoly *term);
struct laurentpoly *laurent_normalize (struct laurentpoly *l);
struct laurentpoly *laurent_dup (struct laurentpoly *l);
void laurent_canonify (struct laurentpoly *l);
void laurent_t_to_oneovert (struct laurentpoly *l);
