#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define RULENAMES_OLD 1
#define RULENAMES_NEW 2

#define ACTION_NONE 0
#define ACTION_PRINTSKETCH 1
#define ACTION_COMPARE 2
#define ACTION_ISCONTOUR 3
#define ACTION_CANONIFY 4
#define ACTION_APPLYRULE 5
#define ACTION_COUNTCC 6
#define ACTION_EXTRACTCC 7
#define ACTION_TESTALLRULES 8
#define ACTION_KNOT2MORSE 9
#define ACTION_PRINTMORSE 10
#define ACTION_CHARACTERISTIC 11
#define ACTION_REMOVECC 12
#define ACTION_INFO 13
#define ACTION_FRONTBACK 14
#define ACTION_LEFTRIGHT 15
#define ACTION_EVERT 16
#define ACTION_MENDES 17
#define ACTION_ISHUFFMAN 18
#define ACTION_MERGEARCS 19
#define ACTION_LISTMA 20
#define ACTION_WRINKLE 21
#define ACTION_LISTWR 22
#define ACTION_SWALLOWTAIL 23
#define ACTION_LISTST 24
#define ACTION_PUNCTURE 25
#define ACTION_LISTPU 26
#define ACTION_ADDSPHERE 27
#define ACTION_PUNCHHOLE 28
#define ACTION_LISTSTRATA 29
#define ACTION_LISTHOLES 30
#define ACTION_REMOVEHOLE 31
#define ACTION_LISTSPHERES 32
#define ACTION_REMOVESPHERE 33
#define ACTION_GLUEARCS 34
#define ACTION_PINCHNECK 35
#define ACTION_CCORIENTATION 36
#define ACTION_FINDCCPARENT 37
#define ACTION_CCORDERING 38
#define ACTION_CCCHILDS 39
#define ACTION_WRAP 40
#define ACTION_3DEVERT 41
#define ACTION_CELLCOMPLEX 42
#define ACTION_FUNDAMENTAL 43
#define ACTION_UNION 44
#define ACTION_CONNECTEDSUM 45
#define ACTION_ANY2MORSE 46
#define ACTION_GIOVECANONIFY 47
#define ACTION_FILEPATH 48
#define ACTION_AFUNDAMENTAL 49
#define ACTION_ALEXANDER 50
#define ACTION_LINKINGNUMBER 51
#define ACTION_SUGGEST_P_SURGERY 52
#define ACTION_FOXJACOBIAN 53
#define ACTION_KNOTSUM 54
#define ACTION_KNOTNAME2DTCODE 55
#define ACTION_KNOTNAME2GAUSSCODE 56

#define MAXKNOTNAMELENGTH 1000

struct global_data {
  int rulenames;
  int focus_on_fundamental;
  int principal;
  int factorideal;
  int internalcheck;
  int abelianize;
  int experimental;
  int userwantscode;
  int knotname_fallback;
  char knotname1[MAXKNOTNAMELENGTH];
  char knotname2[MAXKNOTNAMELENGTH];
};

#define MRPTMAX 4
#define MAPTMAX 4
#define MSPTMAX 4

#define FOXD_UNDEF (-1000)
#define FOXD_MAXINTERESTING (-1001)
#define FOXD_MINVALID (-999)
#define FOXD_KNOT 1
#define FOXD_LINK 1
#define FOXD_GENUS2 2
#define FOXD_GENUS3 3

#define OUTFORMAT_DEFAULT 0
#define OUTFORMAT_APPCONTOUR 1
#define OUTFORMAT_MACAULAY2 2

struct tagged_data {
  int mrnum;
  int region[MRPTMAX];
  int manum;
  int arc[MAPTMAX];
  int arcl[MAPTMAX];
  int stnum;
  int stratum[MSPTMAX];
};
  
/*
 * definizioni per la descrizione tramite regioni
 */

struct sketch {
  struct arc *arcs;
  struct region *regions;
  struct region *extregion;
  int huffman_labelling;
  int isempty;
  int arccount;
  int regioncount;
  int cc_tagged;
  int ccnum;
  int *cc_characteristics;
};

/*
 * please note: the value of f is also twice the
 * link number associated to this region, thus
 * we do not compute it
 */

struct region {
  int f;
  int tag;
  int *strati;
  struct borderlist *border;
  struct region *next;
};

struct borderlist {
  int isexternal;
  struct region *region;
  struct border *sponda;
  struct borderlist *next;
};

struct border {
  int orientation;
  struct arc *info;
  struct borderlist *border;
  struct border *next;
};

struct arc {
  int tag;
  int endpoints;
  int transparent;
  int link_number;
  struct border *regionleft;
  struct border *regionright;
  int cusps;
  int dvalues;
  int depthsdim;
  int *depths;
  struct arc *next;
};

#define F_UNDEF (-9999)

#define TC_UNKNOWN 0
#define TC_EROSION 1
#define TC_DILATION 2

/* definizioni per la lettura di una descrizione 'morse' */

#define ORIENT_LD (-1)
#define ORIENT_RU (1)
#define ORIENT_EMPTY 0

#define INFINITY_ARC (-9999)

#define BIG_INT   10000

/* prototipi */

int checkconsistency (struct sketch *sketch);
int testallrules (struct sketch *sketch);
int testsinglerule (char *rname, int (*rulefunc)(struct sketch *s, int rc), 
                    int exitcode, struct sketch *sketch);
int apply_rule (char *rule, struct sketch *s);
int rule_n14 (struct sketch *s, int rule, int count);
int rule_n5 (struct sketch *s, int count);
int rule_cr2 (struct sketch *s, int count);
int rule_c1 (struct sketch *s, int count);
int rule_c2 (struct sketch *s, int count);
int rule_t1 (struct sketch *s, int count);
int rule_a1 (struct sketch *s, int count);
int rule_a2 (struct sketch *s, int count);
int rule_a12 (struct sketch *s, int count, int ddiff);
int rule_cn1 (struct sketch *s, int count);
int rule_cn1_cr2 (struct sketch *s, int count, int iscr2);
int rule_cn2l (struct sketch *s, int count);
int rule_cn2r (struct sketch *s, int count);
int rule_cn2lb (struct sketch *s, int count);
int rule_cn2rb (struct sketch *s, int count);
int rule_cn2lr (struct sketch *s, int count, int isleft, int isback);
int rule_cn3 (struct sketch *s, int count);
int rule_cr3l (struct sketch *s, int count);
int rule_cr3r (struct sketch *s, int count);
int rule_cr3lr (struct sketch *s, int count, int ori);
int rule_cr1 (struct sketch *s, int count);
int rule_cr1b (struct sketch *s, int count);
int rule_cr11b (struct sketch *s, int count, int isback);
int rule_cr4l (struct sketch *s, int count);
int rule_cr4r (struct sketch *s, int count);
int rule_cr4lb (struct sketch *s, int count);
int rule_cr4rb (struct sketch *s, int count);
int rule_cr4lr (struct sketch *s, int count, int isleft, int isback);
struct border *get_ith_cusp (struct region *r, int i, int *cp);
struct border *rimuovi_losanga (struct region *r, struct sketch *sketch);
void taglia_nodo (struct border *b1n, struct sketch *sketch,
                  struct border **bleft, struct border **bright);
void triple_switch (struct border *b1);
void remove_s1 (struct border *b, struct sketch *sketch);
void spezza_bordo (struct border *cusp, int cpos, struct sketch *sketch,
		 int removecusp);
void pinch_arcs (struct border **bp1pt, int cusp1pos,
                 struct border **bp2pt, int cusp2pos,
                 struct sketch *sketch, int removecusp1, int removecusp2);
/*
 * ora sostituita con la precedente con "removecusps = 1"
 *
void join_cusps (struct border *cusp1, int cusp1pos,
                 struct border *cusp2, int cusp2pos,
                 struct sketch *sketch);
 */
void remove_annulus (struct region *r, struct sketch *sketch);
void remove_ear (struct region *r, struct sketch *sketch);
void applyrulecn2 (struct border *b1n, struct arc *arc,
                   int ori, int orib, struct sketch *sketch);
void remove_cusp (struct region *r, struct sketch *sketch);
void applyrulecr3 (struct border *b1, struct border *b2, int dindex,
                                              struct sketch *sketch);

/* invrules.c */

int gluearcs_or_pinchneck (struct sketch *s, struct arc *a1,
      struct arc *a2, int a1l, int a2l, int ori);
int rule_puncture (struct sketch *s, int rcount);
int list_punctures (struct sketch *s);
int apply_puncture (struct sketch *s, 
		     struct arc *a1, struct arc *a2,
		     int a1l, int a2l, int test);
int rule_createswallowtail (struct sketch *s, int rcount);
int rule_createswallowtailb (struct sketch *s, int rcount);
int list_swallowtails (struct sketch *s);
int apply_createswallowtail (struct sketch *s, struct arc *a,
		     int arcl, int sign, int test);
int rule_createwrinkle (struct sketch *s, int rcount);
int list_strata (struct sketch *s);
int apply_createwrinkle (struct sketch *s, struct region *r,
		     int stratum, int test);
int lookup_mergearcs (char *rule);
int rule_mergearcs (struct sketch *s, int code, int rcount);
int list_mergearcs (struct sketch *s);
int apply_mergearcs (struct sketch *s, struct region *r,
		     struct arc *a1, struct arc *a2,
		     int a1l, int a2l, int test);

int sketchcmp (struct sketch *s1, struct sketch *s2);
int regioncmp (struct region *r1, struct region *r2);
int regiontotcmp (struct region *r1, struct region *r2);
int bordercmp (struct border *b1, struct border *b2);
int bordertotcmp (struct border *b1, struct border *b2);
int singlebordercmp (struct border *b1, struct border *b2);
int singlebordertotcmp (struct border *b1, struct border *b2);
int arccmp (struct arc *a1, struct arc *a2, int depth);
int arctotcmp (struct arc *a1, struct arc *a2, int depth);
int count_hidden_arcs (struct sketch *s1);

/* prototipi appcontour */

int list_hor_sur (struct sketch *s, int ori);
struct region *findregion (struct sketch *s, int tag);
struct arc *findarc (struct sketch *s, int tag);
void showinfo (struct sketch *sketch);
int euler_characteristic (struct sketch *sketch);
int cc_euler_characteristic (struct sketch *sketch);
int appcontourcheck (struct sketch *s, int huffman, int verbose);
int checkorientationborder (struct border *b);
int count_connected_components (struct sketch *sketch);
int connected_component_orientation (int ccid, struct sketch *sketch);
int find_connected_component_parent (int ccid, struct sketch *sketch);
void print_connected_components_ordering (struct sketch *sketch);
void print_connected_component_childs (int ccid, struct sketch *sketch);
int tag_connected_components (struct sketch *sketch);
int free_connected_components (struct sketch *sketch);
int suggest_p_surgery (struct sketch *s, struct region **r, int *stratum);

int extract_connected_components (int ncc, int *ccids, struct sketch *sketch);
//int extract_connected_component (int ccid, struct sketch *sketch);
int remove_connected_component (int ccid, struct sketch *sketch);

int add_s1 (struct sketch *s, struct region *r,
            int stratum, int ori);
int put_in_s1 (struct sketch *s);
int frontback (struct sketch *s);
int leftright (struct sketch *s);
int changeextregion (struct sketch *s, int tag);
int evert3d (struct sketch *s, int tag);
int compute_ohmoto (struct sketch *s);
void compute_link_num_arcs (struct sketch *s);
int count_link_components (struct sketch *s);
int sketch_union (struct sketch *s1, struct sketch *s2);
int sketch_sum (struct sketch *s1, struct sketch *s2);
int sketch_knotsum (struct sketch *s1, struct sketch *s2);

struct region *dupregions (struct region *r);
struct borderlist *dupborderlist (struct borderlist *bl, struct region *r);
struct border *dupborder (struct border *b, struct borderlist *bl);
struct border *dupborder_opened (struct border *b, struct borderlist *bl);
void remregions (struct region *r);
struct arc *reorderbytags (struct arc *r);

void canonify (struct sketch *s);
void oldcanonify (struct sketch *s);
void compat_canonify (struct sketch *s);
void canonifyarc (struct arc *arc);
struct border *canonifyborder (struct border *b);

void sortarcs (struct sketch *s);
struct arc *sortarclist (struct arc *arc);
struct arc *sortequivarcs (struct arc *arc);
struct arc *mergeequivarcs (struct arc *arc, struct arc *rest);
struct borderlist *sortholelist (struct borderlist *hl);
struct region *sortregionlist (struct region *region);

struct sketch *readcontour (FILE *file);
struct sketch *readsketch (FILE *file);
struct sketch *readmorse (FILE *file);
struct sketch *readdtcode (FILE *file);
struct sketch *readgausscode (FILE *file);
struct sketch *readknotscape (FILE *file);
void insert_arc_in_list (struct arc *a, struct arc *list);
void insert_region_in_list (struct region *region, struct region *list);
int readrow (FILE *file, struct sketch *sketch, 
             struct border *actregions[], int num, int vecdim);
int getcrossinfo (FILE *file);
int getarcinfo (int key, FILE *file, 
                struct border *bleft,
                struct border *bright);
int adjustarcinfo (struct arc *narc, struct arc *upperarc, int incr, int oriup);
int *concatenate_depths (struct arc *arcleft, struct arc *arcright);

int arcmult (struct arc *arc);
int cc_arcmult (struct arc *arc, int ccid);
void postprocesssketch (struct sketch *s);
void make_extregion_first (struct sketch *s);
void computefvalue (struct sketch *s, struct region *extregion, int finfinity);
int adjust_isexternalinfo (struct sketch *s);
int iei_process_region (struct region *r);
void defineregionleftright (struct border *border);
void printsketch (struct sketch *s);
void printborder (struct border *b, struct region *r);

int get_d_increase_across_node (struct arc *arc, int ori);
struct border *gettransborder (struct border *b);
struct border *prevborder (struct border *b);
int findinborder (struct border *, struct border *);

struct arc *mergearcs (struct arc *arc1, struct arc *arc2, struct sketch *sketch);
struct arc *mergearcsc (struct arc *arc1, struct arc *arc2, int dincr,
                        struct sketch *sketch);
int topo_change_g (struct border *b1, struct border *b2, 
                   int type, struct sketch *sketch);
int topo_change (struct border *b1, struct border *b2);
void topo_change_l (struct border *b1, struct border *b2);
struct region *regionunion (struct region *r1, struct region *r2, 
                            struct sketch *sketch);
struct borderlist *extractborderlist (struct borderlist *bl);
void redefineregion (struct region *, struct region *);
void redefineborder (struct border *, struct borderlist *);
void transfer_isles (struct borderlist *bl1, struct borderlist *bl2, int hei);

struct sketch *newsketch (void);
struct region *newregion (struct sketch *s);
struct borderlist *newborderlist (struct region *r);
struct borderlist *newborderlist_tail (struct region *r);
struct border *newborder (struct borderlist *bl);
struct arc *newarc (struct sketch *s);

void removearc (struct arc *, struct sketch *);
void removeregion (struct region *, struct sketch *);
struct border *removeborder (struct border *b);
void ensurecanremoveborder (struct border *b);

void freesketch (struct sketch *);
void freeregion (struct region *);
void freeborder (struct border *);
void freeborderdl (struct border *);
void freeborderlist (struct borderlist *);
void freearc (struct arc *);

int readsketch_arc (int arcid, struct sketch *sketch, FILE *file);
int readsketch_region (int regionid, struct sketch *sketch, FILE *file);
struct borderlist *readsketch_bl (struct region *r, struct sketch *sketch, FILE *file);

char check_basic_knotname (char *param);

void knotname2code (FILE *file, int action);

/* knot2morse */
int knot2morse (FILE *file);
int any2morse (FILE *file);

/* printmorse */
void printmorse (struct sketch *sketch);
int morse_ohmoto (struct sketch *sketch);

/* contour */
void ccid_isvalidp (int ncc, int *ccids, int count);

