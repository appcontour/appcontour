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
struct laurentpoly laurent_get_exp_sum (struct presentationrule *r, int g, int gconj);
void print_laurentpoly (struct laurentpoly *l);
