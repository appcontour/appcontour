/*
 * definitions for the computation of the fundamental group
 */

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
    struct ccomplexcc *cc;    // connected components list
  };

struct ccomplexnode {
    int type;
    int stratum;
    int cusp;
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
    struct ccomplexcc *next;
  };

/*
 * prototypes
 */

void compute_fundamental (struct ccomplex *cc);
int compute_fundamental_single (struct ccomplex *cc, struct ccomplexcc *cccc);
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
void cellcomplex_printarcs (struct ccomplex *cc, int verbose);
void cellcomplex_printnodes (struct ccomplex *cc, int verbose);
void cellcomplex_printfaces (struct ccomplex *cc, int verbose);
int onarc2narc (int);
int cellcomplex_checkconsistency (struct ccomplex *cc);
