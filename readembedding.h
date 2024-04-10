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
  struct emb_node *nodes;
};

/* prototypes */

struct embedding *readembedding_low (FILE *file);
struct presentation *wirtingerfromembedding (struct embedding *emb);
int emb_color (struct embedding *emb, int *connections);
int emb_remove_dup_rules (struct presentation *p);
int emb_meridians_longitudes (struct embedding *emb, int *connections, struct presentation *p);
int underpasses_on_spanning_tree (int i, int *node_flood, int *underpasses);
int underpasses_from_to (int i, int j, int *var);
