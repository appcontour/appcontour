
#ifndef FUNDAMENTAL_H_INCLUDED
  #include "fundamental.h"   // needed for structure definition
#endif

#define MAXGENNUM 20
#define REPR_MAXN 20

struct sl2elem {
  int a[2][2];
  };

struct snelem {
  int n;
  int perm[REPR_MAXN];
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
void sl2_set (int m[2][2], int p);
int sl2_iscanon(struct sl2elem *sl2vec, int gennum, int p);
int sl2_isnotcanon(struct sl2elem *sl2vec, int gennum, int p);

int cccountsn (struct presentation *p);
int cccountsn_list (struct presentationlist *pstlist);
void count_sn_cclasses (struct presentationlist *pstlist, int n, int *results);
void sn_init (struct snelem *perm, int n);
int sn_next_cond (struct snelem *perm);
void sn_invert (int *perm, int *perminv, int n);
int sn_isnotcanon(struct snelem *perms, int gennum, int n);
void sn_setlast (int *perm, int n);
int sn_checkrelators(struct snelem *perms, struct snelem *permsinv,
                      struct presentation *pst, int n);
int sn_checkrelator (struct snelem *perms, struct snelem *permsinv, int gennum,
                     struct presentationrule *rule, int n);
void sn_print (int *perm, int n);
int sn_nextmap (struct snelem *perms, struct snelem *permsinv, int gennum);
int sn_next (int *perm, int n);
int sn_permcmp (int *p1, int *p2, int n);
int sn_iseven (int *perm, int n);
void sn_permmul (int *acc, int *perm, int n);

/* various utilities */

int isprime (int p);
int isprimepower (int q);
