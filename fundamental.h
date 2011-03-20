/*
 * definitions for the computation of the fundamental group
 */

#define FG_SURFACE 0
#define FG_INTERNAL 1
#define FG_EXTERNAL (-1)
#define CC_NODETYPE_FOLD 1
#define CC_NODETYPE_CUT 2
#define CC_NODETYPE_VIRTUALFOLD 3
#define CC_NODETYPE_VIRTUALCUT 4

#define CC_ARCTYPE_CUT 1
#define CC_ARCTYPE_FOLD 2

struct ccomplex {
    struct sketch *sketch;
    int type;
    int nodenum;
    struct ccomplexnode *nodes;
    int arcnum;
    struct ccomplexarc *arcs;
    int facenum;
    int ccnum;
    struct ccomplexcc *cc;    // connected components list
  };

struct ccomplexnode {
    int type;
    int stratum;
    struct ccomplexcc *cc;
    struct arc *under;
    struct arc *over;
  };

struct ccomplexarc {
    int type;
    int enda;
    int endb;
    int isinspanningtree;
    int stratum;
    struct arc *arc;
  };

struct ccomplexcc {
    int tag;
    int basenode;
    struct ccomplexcc *next;
  };

/*
 * prototypes
 */

struct ccomplex *compute_fundamental (struct sketch *s, int which);
int fundamental_countnodes (struct sketch *s);
int fundamental_countarcs (struct sketch *s, int which);
int fundamental_countfaces (struct sketch *s, int which);
void fundamental_fillnodes (struct ccomplex *cc);
void fundamental_fillarcs (struct ccomplex *cc);
int find_spanning_tree (struct ccomplex *cc);
