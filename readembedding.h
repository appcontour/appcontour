struct emb_node {
  int id;
  int valency;
  struct emb_node *next;
  int ping[4];
  int direction[4];
  int generator[4];
  int pong[4]; /* this is the index where it hits at the arrival node */
  int overpassisodd;
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
