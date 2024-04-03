struct emb_node {
  int id;
  int valency;
  struct emb_node *next;
  int adj[4];
};

struct embedding {
  int k;  /* number of trivalent nodes */
  int n;  /* number of crossings (four-valent nodes) */
  int choice;
  struct emb_node *nodes;
  struct emb_node *crossings;
};
