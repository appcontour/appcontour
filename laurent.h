/*
 * definitions for alexander polynomial
 */

struct laurentpoly {
  int minexpon;
  int stemdegree;
  int stem[];
};

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

#define POLYXSIZE(n) (sizeof(struct laurentpolyx) + (n)*sizeof(union intorpointer))

void laurent_negate (struct laurentpoly *term);
void laurent_negatex (struct laurentpolyx *term);
void print_laurentpoly (struct laurentpoly *l, char indet);
void print_laurentpolyx (struct laurentpolyx *l, char *indetlist);
struct laurentpoly *read_laurentpoly (FILE *file, char indet_names[2]);
struct laurentpolyx *read_laurentpolyx (FILE *file, char indet_names[2]);
struct laurentpoly *laurent_add (struct laurentpoly *add1, struct laurentpoly *add2);
struct laurentpolyx *laurent_addx (struct laurentpolyx *add1, struct laurentpolyx *add2);
struct laurentpoly *laurent_mul (struct laurentpoly *fact1, struct laurentpoly *fact2);
struct laurentpolyx *laurent_mulx1 (struct laurentpolyx *fact1, struct laurentpolyx *fact2);
struct laurentpolyx *laurent_mulx2 (struct laurentpolyx *fact1, struct laurentpolyx *fact2);
struct laurentpoly *laurent_normalize (struct laurentpoly *l);
struct laurentpolyx *laurent_normalizex (struct laurentpolyx *l);
struct laurentpoly *laurent_dup (struct laurentpoly *l);
struct laurentpolyx *laurent_dupx (struct laurentpolyx *l);
void laurent_canonify (struct laurentpoly *l);
void laurent_canonifyx1 (struct laurentpolyx *l);
void laurent_canonifyx (struct laurentpolyx *l);
void laurent_t_to_oneovert (struct laurentpoly *l);
void laurentx1_t_to_oneovert (struct laurentpolyx *l);
void free_laurentpolyx (struct laurentpolyx *l);
struct laurentpoly *laurentpoly_addmonom (struct laurentpoly *l, int expon, int coef);
struct laurentpolyx *laurentpolyx1_addmonom (struct laurentpolyx *l, int expon, int coef);
struct laurentpolyx *laurentpolyx2_addmonom (struct laurentpolyx *l, int expu, int expv, int coef);
struct laurentpolyx *laurentpolyx_addmonom (struct laurentpolyx *l, int rank, int *expuv, int coef);
int laurent_sum_coefficients (struct laurentpoly *l);
int laurent_sum_coefficientsx1 (struct laurentpolyx *l);
struct laurentpolyx *laurent_sum_coefficientsx2 (struct laurentpolyx *l);
struct laurentpolyx *laurent_sum_each_coefficientx2 (struct laurentpolyx *l);

struct laurentpoly *laurent_gcd (int inspread, struct laurentpoly *p1, struct laurentpoly *p2, int *spreadpt);
int laurent_factor_content (struct laurentpoly *p);
struct laurentpoly *laurent_euclid (struct laurentpoly *p1, struct laurentpoly *p2);
struct laurentpoly *laurent_extended_euclid (struct laurentpoly *p1, struct laurentpoly *p2);

int laurentx_totdegree (struct laurentpolyx *l);
int laurentx_lexicocompare (struct laurentpolyx *p1, struct laurentpolyx *p2);
void laurent_canonifysign (struct laurentpoly *p);
void laurentx1_canonifysign (struct laurentpolyx *p);
void laurentx_canonifysign (struct laurentpolyx *p);

int gcd (int a, int b);
