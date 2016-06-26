/*
 * definitions for alexander polynomial
 */

#define IDEAL_DEF_GENERATORS_NUM 10

/*
 * if indets == 0, then val holds the constant value
 * if indets == 1, then l[l1num] is a vector with
 * generators of the ideal (for now l1num must be 1)
 * if indets == 2, then l[l2num] is a vector with
 * generic generators, whereas fl[fl2num] is a vector
 * with generators that must be multiplied by the
 * fundamental ideal (it is contained in l[] with
 * an offset (fl2offset)
 */

struct alexanderideal {
  int indets;
  int l1num;
  int spread;
  int l2num;
  int fl2offset;
  int fl2num;
  int val;
  int max_generators_num;
  struct laurentpoly *l[]; /* followed by ex fl2 at fl2offset */
};

#define AI_DIM0 (sizeof(struct alexanderideal))
#define AI_DIMX (sizeof(struct laurentpollyx *))
#define AI_DIM(size) ((size)*AI_DIMX + AI_DIM0)

/*
 * laurentmatrix (generic number of indets)
 */

struct laurentmatrix {
  int numrows;
  int numcols;
  struct laurentpoly ***columns;
};

/* prototypes */

int alexander (struct presentation *p);
int alexander_fromideal (struct alexanderideal *ai);
void printout_ideal1 (struct alexanderideal *ai, struct laurentpoly *principal);
void printout_ideal (struct alexanderideal *ai, struct laurentpoly *principal,
                     struct laurentpoly **fundamentals, int fnum, int printprincipal);
int printout_constant_ideal (char *msg, int val);
int linkingnumber (struct presentation *p);
int linkingnumber_fromideal (struct alexanderideal *ai);
int corank_one_alexander (struct presentation *p);
struct laurentpoly *laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate);
struct alexanderideal *laurent_notfirst_elementary_ideal (struct presentation *p, int eliminate, int corank);
struct laurentpoly *laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2,
                      struct laurentpoly ***extrapt);
struct alexanderideal *laurent_notfirst_elementary_ideal2 (struct presentation *p, int e1, int e2, int corank);
struct laurentpoly *laurent_get_exp_sum1 (struct presentationrule *r, int g, int gconj);
struct laurentpoly *laurent_get_exp_sum (struct presentationrule *r, int g, int e1, int e2);
struct laurentpoly *laurent_mixed_derivative2 (struct presentationrule *r, int x1, int x2);
struct laurentpoly *laurent_common_factor (struct presentationrule *r, int x1, int x2);
struct laurentpoly *laurent_compute_determinant1 (struct laurentpoly ***matrix, int n);
struct laurentpoly *laurent_compute_determinant (struct laurentpoly ***matrix, int n, int indets);
struct laurentpoly *laurent_minor_determinant (struct laurentpoly ***matrix, int n,
                          int row_to_substitute, int indets);
int canonify_ideal2 (struct laurentpoly **lpt, struct laurentpoly **extradets, int extranum);
struct laurentmatrix *laurent_build_matrix1 (struct presentation *p, int eliminate);
struct laurentmatrix *laurent_build_matrix2 (struct presentation *p, int e1, int e2);
struct laurentmatrix *minor_matrix1_corank1 (struct laurentmatrix *matrix, int i, int j);
struct laurentmatrix *minor_matrix2_corank1 (struct laurentmatrix *matrix, int i, int j);
void laurent_free_matrix (struct laurentmatrix *matrix);
struct alexanderideal *laurent_simplify_ideal (struct alexanderideal *ai);
int laurent_try_simplify_ideal (struct alexanderideal *ai);

int laurent_suppdimx2 (struct laurentpoly *l);
void laurent_getthree (struct laurentpoly *l, int *x1, int *x2, int *x3);
int base_canonify2 (struct laurentpoly **lpt);
int base_canonify2_onedim (struct laurentpoly **lpt);
int base_canonify2_twodim (struct laurentpoly **lpt);
void shuffle_poly2 (struct laurentpoly **lpt, struct laurentpoly **extradets, int extranum);
struct laurentpoly *base_change2 (struct laurentpoly *l, int matrixb[2][2]);

void base_random (int matrixb[2][2]);
int isinvertible_base (int b[2][2]);

struct alexanderideal *read_alexander_ideal (FILE *file);
struct alexanderideal *ai_increase_size (struct alexanderideal *ai);
struct alexanderideal *compute_invariant_factor (struct laurentpoly ***columns, int numrows, int numcols, int minordim, int indets);

void print_matrix1 (struct laurentmatrix *matrix);
void print_matrix (struct laurentmatrix *matrix, int indets);
