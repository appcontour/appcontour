struct emb_node {
  int id;
  int valency;
  struct emb_node *next;
  int ping[4]; /* these are the adjacent nodes, ordered counterclockwise */
  int direction[4];
  int generator[4];
  int pong[4]; /* this is the index where it hits at the arrival node */
  int overpassisodd;
  short int color;
  short int colorodd;
};

struct embedding {
  int k;  /* number of trivalent nodes */
  int n;  /* number of crossings (four-valent nodes) */
  int choice;
  int numhcomponents;  /* number of connected components that contain trivalent nodes */
  int numrings;
  int orientation;
  int *connections;
  struct emb_node *nodes;
};

struct dual_region {
  int id;
  int valency;
  struct dual_region **ping; /* these are the adjacent regions, ordered counterclockwise */
  int *wedgeij;  /* pointer to the vector of wedges */
  struct dual_region *next;
  // int pingpong[];  /* valency is the size of this variable size portion */
};

struct dualembedding {
  int numregions;
  int v;
  int e;
  int numwedges;
  int *wedgeij;
  struct dual_region *regions;
};

/* prototypes */

struct sketch *embedding2sketch (struct embedding *emb);
struct dualembedding *embedding2dual (struct embedding *emb);
void printdualembedding (struct dualembedding *dual, struct embedding *emb);
void printembedding (struct embedding *emb);
void print_dual_type (struct dualembedding *dual, struct embedding *emb);
void freedualembedding (struct dualembedding *dual);
void freedualregions (struct dual_region *region);
struct embedding *readembedding (FILE *file);
struct presentation *wirtingerfromembedding (struct embedding *emb);
struct presentation *ccasloop (struct embedding *emb);
int isexcluded (int color, int loop);
struct presentation *wirtingerfromembeddingraw (struct embedding *emb);
int emb_color (struct embedding *emb);
int emb_color4 (struct embedding *emb);
int emb_remove_dup_rules (struct presentation *p);
int emb_meridians_longitudes (struct embedding *emb, struct presentation *p, int cc);
int emb_meridians_longitudes_torus (struct embedding *emb, struct presentation *p, int cc);
int numunderpasses_on_spanning_tree (int i, int *node_flood, int *underpasses);
int underpasses_on_arc (int i_and_k, int *var, struct embedding *emb);
int emb_orient (struct embedding *emb);
void printembrules (struct embedding *emb, struct dualembedding *dual);
struct vecofintlist *embeddingtoloiv (struct embedding *emb);
void freeembedding (struct embedding *emb);
int embedding_connectedness (struct dualembedding *dual, struct embedding *emb);
