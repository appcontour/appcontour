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
  struct laurentpolyx *l2;
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

struct laurentmatrixx {
  int numrows;
  int numcols;
  struct laurentpolyx ***columns;
};

/* prototypes */

int alexander (struct presentation *p);
int alexander_fromideal (struct alexanderideal *ai);
void printout_ideal1 (struct alexanderideal *ai, struct laurentpoly *principal);
void printout_idealx (struct alexanderideal *ai, struct laurentpolyx *principal,
                     struct laurentpolyx **fundamentals, int fnum, int printprincipal);
int printout_constant_ideal (char *msg, int val);
int linkingnumber (struct presentation *p);
int linkingnumber_fromideal (struct alexanderideal *ai);
int corank_one_alexander (struct presentation *p);
struct laurentpoly *laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate);
struct alexanderideal *laurent_notfirst_elementary_ideal (struct presentation *p, int eliminate, int corank);
struct laurentpolyx *laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2,
                      struct laurentpolyx ***extrapt);
struct alexanderideal *laurent_notfirst_elementary_ideal2 (struct presentation *p, int e1, int e2, int corank);
struct laurentpoly *laurent_get_exp_sum (struct presentationrule *r, int g, int gconj);
struct laurentpolyx *laurent_get_exp_sumx (struct presentationrule *r, int g, int e1, int e2);
struct laurentpolyx *laurent_mixed_derivativex2 (struct presentationrule *r, int x1, int x2);
struct laurentpolyx *laurent_common_factorx (struct presentationrule *r, int x1, int x2);
struct laurentpoly *laurent_compute_determinant (struct laurentpoly ***matrix, int n);
struct laurentpolyx *laurent_compute_determinantx (struct laurentpolyx ***matrix, int n);
struct laurentpolyx *laurent_minor_determinantx (struct laurentpolyx ***matrix, int n,
                          int row_to_substitute);
int canonify_idealx (struct laurentpolyx **lpt, struct laurentpolyx **extradets, int extranum);
struct laurentmatrix *laurent_build_matrix (struct presentation *p, int eliminate);
struct laurentmatrixx *laurent_build_matrixx (struct presentation *p, int e1, int e2);
struct laurentmatrix *minor_matrix_corank1 (struct laurentmatrix *matrix, int i, int j);
struct laurentmatrixx *minor_matrixx_corank1 (struct laurentmatrixx *matrix, int i, int j);
void laurent_free_matrix (struct laurentmatrix *matrix);
void laurent_free_matrixx (struct laurentmatrixx *matrix);
struct alexanderideal *laurent_simplify_ideal (struct alexanderideal *ai);
int laurent_try_simplify_ideal (struct alexanderideal *ai);

int laurent_suppdimx2 (struct laurentpolyx *l);
void laurent_getthreex (struct laurentpolyx *l, int *x1, int *x2, int *x3);
int base_canonifyx2 (struct laurentpolyx **lpt);
int base_canonifyx2_onedim (struct laurentpolyx **lpt);
int base_canonifyx2_twodim (struct laurentpolyx **lpt);
void shuffle_polyx (struct laurentpolyx **lpt, struct laurentpolyx **extradets, int extranum);
struct laurentpolyx *base_changex (struct laurentpolyx *l, int matrixb[2][2]);

void base_random (int matrixb[2][2]);
int isinvertible_base (int b[2][2]);

struct alexanderideal *read_alexander_ideal (FILE *file);
struct alexanderideal *ai_increase_size (struct alexanderideal *ai);
