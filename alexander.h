/*
 * definitions for alexander polynomial
 */

struct laurentpoly {
  int minexpon;
  int stemdegree;
  int denom;
  int stem[];
};

struct laurentpoly2 {
  int minexpon;
  int stemdegree;
  int denom;
  struct laurentpoly *stem[];
};

/* prototypes */

int alexander (struct presentation *p);
int linkingnumber (struct presentation *p);
int corank_one_alexander (struct presentation *p);
struct laurentpoly *laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate);
struct laurentpoly2 *laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2,
                      struct laurentpoly2 ***extrapt);
struct laurentpoly *laurent_get_exp_sum (struct presentationrule *r, int g, int gconj);
struct laurentpoly2 *laurent_get_exp_sum2 (struct presentationrule *r, int g, int e1, int e2);
struct laurentpoly2 *laurent_mixed_derivative2 (struct presentationrule *r, int x1, int x2);
struct laurentpoly2 *laurent_common_factor2 (struct presentationrule *r, int x1, int x2);
struct laurentpoly *laurent_compute_determinant (struct laurentpoly ***matrix, int n);
struct laurentpoly2 *laurent_compute_determinant2 (struct laurentpoly2 ***matrix, int n);
struct laurentpoly2 *laurent_minor_determinant2 (struct laurentpoly2 ***matrix, int n,
                          int row_to_substitute);
void print_laurentpoly (struct laurentpoly *l, char indet);
void print_laurentpoly2 (struct laurentpoly2 *l, char indet1, char indet2);
struct laurentpoly *laurent_add (struct laurentpoly *add1, struct laurentpoly *add2);
struct laurentpoly2 *laurent_add2 (struct laurentpoly2 *add1, struct laurentpoly2 *add2);
struct laurentpoly *laurent_mul (struct laurentpoly *fact1, struct laurentpoly *fact2);
struct laurentpoly2 *laurent_mul2 (struct laurentpoly2 *fact1, struct laurentpoly2 *fact2);
void laurent_negate (struct laurentpoly *term);
void laurent_negate2 (struct laurentpoly2 *term);
struct laurentpoly *laurent_normalize (struct laurentpoly *l);
struct laurentpoly2 *laurent_normalize2 (struct laurentpoly2 *l);
struct laurentpoly *laurent_dup (struct laurentpoly *l);
struct laurentpoly2 *laurent_dup2 (struct laurentpoly2 *l);
void laurent_canonify (struct laurentpoly *l);
void laurent_canonify2 (struct laurentpoly2 *l);
int canonify_ideal2 (struct laurentpoly2 **lpt, struct laurentpoly2 **extradets, int extranum);
void laurent_t_to_oneovert (struct laurentpoly *l);
void free_laurentpoly2 (struct laurentpoly2 *l);
struct laurentpoly *laurentpoly_addmonom (struct laurentpoly *l, int expon, int coef);

int laurent_sum_coefficients (struct laurentpoly *l);
struct laurentpoly *laurent_sum_coefficients2 (struct laurentpoly2 *l);
struct laurentpoly *laurent_sum_each_coefficient2 (struct laurentpoly2 *l);
int laurent_suppdim2 (struct laurentpoly2 *l);
int base_canonify2 (struct laurentpoly2 **lpt);
int base_canonify2_onedim (struct laurentpoly2 **lpt);
int base_canonify2_twodim (struct laurentpoly2 **lpt);

int mcd (int a, int b);
