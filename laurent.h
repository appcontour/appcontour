/*
 * definitions for alexander polynomial
 */

union intorpointer {
  struct laurentpolyx *lx;
  int                  l0;
};

struct laurentpolyx {
  int indets;
  int minexpon;
  int stemdegree;
  union intorpointer stem[];
};

#define POLYSIZE(n) (sizeof(struct laurentpolyx) + (n)*sizeof(union intorpointer))

void laurent_negate (struct laurentpolyx *term);
void print_laurentpoly (struct laurentpolyx *l, char *indetlist);
struct laurentpolyx *read_laurentpoly1 (FILE *file, char indet_names[2]);
struct laurentpolyx *read_laurentpoly2 (FILE *file, char indet_names[2]);
struct laurentpolyx *laurent_add (struct laurentpolyx *add1, struct laurentpolyx *add2);
struct laurentpolyx *laurent_mul1 (struct laurentpolyx *fact1, struct laurentpolyx *fact2);
struct laurentpolyx *laurent_mul2 (struct laurentpolyx *fact1, struct laurentpolyx *fact2);
struct laurentpolyx *laurent_normalize (struct laurentpolyx *l);
struct laurentpolyx *laurent_dup (struct laurentpolyx *l);

/* TODO from here clean up mess with old laurentpoly struct */

void laurent_canonify (struct laurentpolyx *l);
void laurent_canonifyx1 (struct laurentpolyx *l);
void laurent_canonifyx (struct laurentpolyx *l);
void laurent_t_to_oneovert (struct laurentpolyx *l);
void laurentx1_t_to_oneovert (struct laurentpolyx *l);
void free_laurentpolyx (struct laurentpolyx *l);
struct laurentpolyx *laurentpoly_addmonom (struct laurentpolyx *l, int expon, int coef);
struct laurentpolyx *laurentpolyx1_addmonom (struct laurentpolyx *l, int expon, int coef);
struct laurentpolyx *laurentpolyx2_addmonom (struct laurentpolyx *l, int expu, int expv, int coef);
struct laurentpolyx *laurentpolyx_addmonom (struct laurentpolyx *l, int rank, int *expuv, int coef);
int laurent_sum_coefficients (struct laurentpolyx *l);
int laurent_sum_coefficientsx1 (struct laurentpolyx *l);
struct laurentpolyx *laurent_sum_coefficientsx2 (struct laurentpolyx *l);
struct laurentpolyx *laurent_sum_each_coefficientx2 (struct laurentpolyx *l);

struct laurentpolyx *laurent_gcd (int inspread, struct laurentpolyx *p1, struct laurentpolyx *p2, int *spreadpt);
int laurent_factor_content (struct laurentpolyx *p);
struct laurentpolyx *laurent_euclid (struct laurentpolyx *p1, struct laurentpolyx *p2);
struct laurentpolyx *laurent_extended_euclid (struct laurentpolyx *p1, struct laurentpolyx *p2);

int laurentx_totdegree (struct laurentpolyx *l);
int laurentx_lexicocompare (struct laurentpolyx *p1, struct laurentpolyx *p2);
void laurent_canonifysign (struct laurentpolyx *p);
void laurentx1_canonifysign (struct laurentpolyx *p);
void laurentx_canonifysign (struct laurentpolyx *p);

int gcd (int a, int b);
