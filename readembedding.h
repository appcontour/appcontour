struct emb_node {
  int id;
  int valency;
  struct emb_node *next;
  int ping[4];
  int direction[4];
  int generator[4];
  int pong[4]; /* this is the index where it hits at the arrival node */
  int overpassisodd;
  int color;
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

/* prototypes */

struct sketch *embeddingtosketch (FILE *file);
struct embedding *readembedding (FILE *file);
struct presentation *wirtingerfromembedding (struct embedding *emb);
int emb_color (struct embedding *emb);
int emb_remove_dup_rules (struct presentation *p);
int emb_meridians_longitudes (struct embedding *emb, struct presentation *p);
int numunderpasses_on_spanning_tree (int i, int *node_flood, int *underpasses);
int underpasses_on_arc (int i_and_k, int *var, struct embedding *emb);
int emb_orient (struct embedding *emb);
struct vecofintlist *embeddingtoloiv (struct embedding *emb);
