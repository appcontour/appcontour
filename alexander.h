/*
 * definitions for alexander polynomial
 */

#define IDEAL_DEF_GENERATORS_NUM 10

/*
 * if indets == 0, then val holds the constant value
 * if indets == 1, then l1[l1num] is a vector with
 * generators of the ideal (for now l1num must be 1)
 * if indets == 2, then l2[l2num] is a vector with
 * generic generators, whereas fl2[fl2num] is a vector
 * with generators that must be multiplied by the
 * fundamental ideal
 */

union lxunion {
  struct laurentpoly *l1;
  struct laurentpoly2 *l2;
};

struct alexanderideal {
  int indets;
  int l1num;
  int spread;
  int l2num;
  int fl2offset;
  int fl2num;
  int val;
  int max_generators_num;
  union lxunion lx[]; /* followed by ex fl2 at fl2offset */
  //struct laurentpoly *l1[IDEAL_MAX_GENERATORS_NUM];
  //struct laurentpoly2 *l2[IDEAL_MAX_GENERATORS_NUM];  /* followed by ex fl2 at fl2offset */
  //struct laurentpoly2 *fl2[IDEAL_MAX_GENERATORS_NUM];
};

#define AI_DIM0 (sizeof(struct alexanderideal))
#define AI_DIMX (sizeof(union lxunion))
#define AI_DIM(size) ((size)*AI_DIMX + AI_DIM0)

/*
 * laurentmatrix
 */

struct laurentmatrix {
  int numrows;
  int numcols;
  struct laurentpoly ***columns;
};

/*
 * laurentmatrix in two indeterminates
 */

struct laurentmatrix2 {
  int numrows;
  int numcols;
  struct laurentpoly2 ***columns;
};

/* prototypes */

int alexander (struct presentation *p);
int alexander_fromideal (struct alexanderideal *ai);
void printout_ideal1 (struct alexanderideal *ai, struct laurentpoly *principal);
void printout_ideal2 (struct alexanderideal *ai, struct laurentpoly2 *principal,
                     struct laurentpoly2 **fundamentals, int fnum, int printprincipal);
int printout_constant_ideal (char *msg, int val);
int linkingnumber (struct presentation *p);
int linkingnumber_fromideal (struct alexanderideal *ai);
int corank_one_alexander (struct presentation *p);
struct laurentpoly *laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate);
struct alexanderideal *laurent_notfirst_elementary_ideal (struct presentation *p, int eliminate, int corank);
struct laurentpoly2 *laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2,
                      struct laurentpoly2 ***extrapt);
struct alexanderideal *laurent_notfirst_elementary_ideal2 (struct presentation *p, int e1, int e2, int corank);
struct laurentpoly *laurent_get_exp_sum (struct presentationrule *r, int g, int gconj);
struct laurentpoly2 *laurent_get_exp_sum2 (struct presentationrule *r, int g, int e1, int e2);
struct laurentpoly2 *laurent_mixed_derivative2 (struct presentationrule *r, int x1, int x2);
struct laurentpoly2 *laurent_common_factor2 (struct presentationrule *r, int x1, int x2);
struct laurentpoly *laurent_compute_determinant (struct laurentpoly ***matrix, int n);
struct laurentpoly2 *laurent_compute_determinant2 (struct laurentpoly2 ***matrix, int n);
struct laurentpoly2 *laurent_minor_determinant2 (struct laurentpoly2 ***matrix, int n,
                          int row_to_substitute);
int canonify_ideal2 (struct laurentpoly2 **lpt, struct laurentpoly2 **extradets, int extranum);
struct laurentmatrix *laurent_build_matrix (struct presentation *p, int eliminate);
struct laurentmatrix2 *laurent_build_matrix2 (struct presentation *p, int e1, int e2);
struct laurentmatrix *minor_matrix_corank1 (struct laurentmatrix *matrix, int i, int j);
struct laurentmatrix2 *minor_matrix2_corank1 (struct laurentmatrix2 *matrix, int i, int j);
void laurent_free_matrix (struct laurentmatrix *matrix);
void laurent_free_matrix2 (struct laurentmatrix2 *matrix);
struct alexanderideal *laurent_simplify_ideal (struct alexanderideal *ai);
int laurent_try_simplify_ideal (struct alexanderideal *ai);

int laurent_suppdim2 (struct laurentpoly2 *l);
void laurent_getthree2 (struct laurentpoly2 *l, int *x1, int *x2, int *x3);
int base_canonify2 (struct laurentpoly2 **lpt);
int base_canonify2_onedim (struct laurentpoly2 **lpt);
int base_canonify2_twodim (struct laurentpoly2 **lpt);
void shuffle_poly2 (struct laurentpoly2 **lpt, struct laurentpoly2 **extradets, int extranum);
struct laurentpoly2 *base_change2 (struct laurentpoly2 *l, int matrixb[2][2]);

void base_random (int matrixb[2][2]);
int isinvertible_base (int b[2][2]);

struct alexanderideal *read_alexander_ideal (FILE *file);
struct alexanderideal *ai_increase_size (struct alexanderideal *ai);
