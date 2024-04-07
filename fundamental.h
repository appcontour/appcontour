/*
 * definitions for the computation of the fundamental group
 */

#define FUNDAMENTAL_H_INCLUDED 1
#define FG_UNDEF (-2)
#define FG_SURFACE 0
#define FG_INTERNAL 1
#define FG_EXTERNAL (-1)

#define CC_REMOVED 9999

/* from here... */
#define CC_NODETYPE_FOLD 1
#define CC_NODETYPE_CUT 2
#define CC_NODETYPE_VIRTUALFOLD 3
#define CC_NODETYPE_VIRTUALCUT 4
/* ...to here must be kept contiguous */
/* see function fund_findnode */
#define CC_NODETYPE_CUSP 5

#define CC_ARCTYPE_CUT 1
#define CC_ARCTYPE_FOLD 2
#define CC_ARCTYPE_VIRTUAL 3
#define CC_ARCTYPE_COLUMN 4
#define CC_ARCTYPE_VCOLUMN 5

#define CC_FACETYPE_HORIZONTAL 1
#define CC_FACETYPE_WALL 2

#define BETTI_UNDEF (-1)

struct ccomplex {
    struct sketch *sketch;
    int type;
    int nodenum;
    int arcnum;
    int facenum;
    int nodedim;
    struct ccomplexnode *nodes;
    int arcdim;
    struct ccomplexarc *arcs;
    int facedim;
    struct ccomplexface *faces;
    int ccnum;
    int surfccnum;
    int *cc_characteristics;
    struct ccomplexcc *cc;    // connected components list
  };

struct ccomplexnode {
    int type;
    int stratum;
    int cusp;
    int surfcc;
    struct ccomplexcc *cc;
    struct arc *ne;   // with both arcs oriented from left to right
    struct arc *se;
    int refcount;     // number of arcs adjacent to this node
  };

struct ccomplexarc {
    int type;
    int enda;
    int endb;
    int cusp1;
    int cusp2;
    int isinspanningtree;
    int stratum;
    struct arc *arc;
    struct borderlist *bl;  // for virtual arcs, this points to the connecting island
    int refcount;           // number of faces adjacent to this arc
  };

/*
 * integer vector faceborder is dimensioned facebordernum
 * and contains positive values "na" for positively oriented arcs
 * negative values for negatively oriented arcs
 * the integer value |na|-1 points to the actual arc in the
 * cell structure
 */

struct ccomplexface {
    int type;
    int facebordernum;
    int *faceborder;
    int stratum;
  };

struct ccomplexcc {
    int tag;
    int basenode;
    int betti_2;
    int spherical_voids;
    struct presentation *p;
    struct ccomplexcc *next;
  };

/*
 * data structure for group presentation
 * it can contain a list of selected elements (meridians, longitudes,...)
 */

struct presentation {
    int gennum;
    int characteristic;
    int espected_deficiency;
    struct presentationrule *rules;
    struct presentationrule *elements;
  };

struct presentationlist {
    struct presentation *p;
    struct presentationlist *next;
  };

struct presentationrule {
    int length;
    struct presentationrule *next;
    int var[];   // negative for inverses
  };

/*
 * prototypes
 */

void compute_fundamental (struct ccomplex *cc, int action);
void fundamental_group (struct presentation *p);
void print_deficiency (struct presentation *p);
void abelianized_fundamental_group (struct presentation *p);
int complex_characteristic (struct ccomplex *cc);
struct presentation *compute_fundamental_single (struct ccomplex *cc, struct ccomplexcc *cccc);
void print_presentation (struct presentation *p);
void print_rule_list (struct presentationrule *r, int gennum);
void print_single_rule (struct presentationrule *r, int gennum);
void print_invariant_factors (struct presentation *p);
void print_exponent_matrix (struct presentation *p);
void fg_interactive (struct presentation *p);
int simplify_presentation (struct presentation *p);
int sp_eliminatevar (struct presentation *p);
int sp_removeemptyrules (struct presentation *p);
void topreabelian (struct presentation *p);
int compute_fg_rank (struct presentation *p);
void remove_all_relators (struct presentation *p);
void read_group_presentation (FILE *file, struct presentation *p);
struct presentationlist *read_group_presentation_list (FILE *file);
int get_exp_sum (struct presentationrule *r, int n);
int complex_melt (struct ccomplex *cc);
int complex_facemelt (struct ccomplex *cc);
int complex_faceremovekink (struct ccomplex *cc);
void complex_do_melt_faces (struct ccomplex *cc, int, int, int);
void complex_do_removekink (struct ccomplex *cc, int, int, int, int);
int complex_collapse (struct ccomplex *cc);
int complex_collapse_faces (struct ccomplex *cc);
int complex_collapse_arcs (struct ccomplex *cc);
void complex_remove_face (struct ccomplex *cc, int n);
void complex_remove_face_nd (struct ccomplex *cc, int n);
void complex_remove_arc (struct ccomplex *cc, int n);
void complex_remove_node (struct ccomplex *cc, int n);
void complex_countreferences (struct ccomplex *cc);
struct ccomplex *compute_cellcomplex (struct sketch *s, int which);
int fundamental_countnodes (struct sketch *s);
int fundamental_countarcs (struct sketch *s, int which);
int fundamental_countfaces (struct sketch *s, int which);
void fundamental_fillnodes (struct ccomplex *cc);
void fundamental_fillarcs (struct ccomplex *cc);
void fundamental_fillfaces (struct ccomplex *cc);
int fund_findnode (struct ccomplex *cc, struct arc *a, int stratum);
int find_spanning_tree (struct ccomplex *cc);
void cc_revert_face (struct ccomplex *cc, int nface);
void cellcomplex_print (struct ccomplex *cc, int verbose);
void cellcomplex_printarcs (struct ccomplex *cc, int verbose, int *noderemap, int *arcremap);
void cellcomplex_printnodes (struct ccomplex *cc, int verbose, int *noderemap);
void cellcomplex_printfaces (struct ccomplex *cc, int verbose, int *arcremap);
int onarc2narc (int);
int cellcomplex_checkconsistency (struct ccomplex *cc);

int read_generators_list (FILE *file, char *gennames, int maxgennum);
