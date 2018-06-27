
#ifndef FUNDAMENTAL_H_INCLUDED
  #include "fundamental.h"   // needed for structure definition
#endif

#define MAXGENNUM 20

struct sl2elem {
  int a[2][2];
  };

int cccountsl2zp (struct presentation *p);
int count_sl2zp_cclasses (struct presentation *pst, int p);
void sl2_clear (int m[2][2]);
int sl2_next_det1 (int m[2][2], int p);
int sl2_next (int m[2][2], int p);
int sl2_det (int m[2][2], int p);
void sl2_invert (int m[2][2], int minv [2][2], int p);
int inv_modp (int n, int p);
void sl2_print (int m[2][2]);
int sl2_nextmap (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, int gennum, int p);
int sl2_checkrelators(struct sl2elem *sl2vec, struct sl2elem *sl2vecinv,
                      struct presentation *pst, int p);
int sl2_checkrelator (struct sl2elem *sl2vec, struct sl2elem *sl2vecinv, int gennum,
                  struct presentationrule *rule, int p);
void sl2_matmul (int acc[2][2], int mat2[2][2], int p);
void sl2_matcopy (int acc[2][2], int mat[2][2]);
int sl2_matcmp (int m1[2][2], int m2[2][2]);
int sl2_iscanon(struct sl2elem *sl2vec, int gennum, int p);
