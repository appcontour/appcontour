#define MAXDTCODELEN 200

struct vecofintlist {
  int len;
  int dim;
  struct vecofintlist *next;
  int *handedness;
  int vec[];
};

#define SIZEOFLOIV(n) (sizeof(struct vecofintlist) + n*sizeof(int))

/* prototypes */

struct vecofintlist *readnakedvecofintlist (FILE *file);
struct vecofintlist *readvecofintlist (FILE *file);
void freeloiv (struct vecofintlist *loiv);
struct sketch *realize_dtcode (int numnodes, int *vecofint, int *gregionsign);
struct sketch *realize_gausscode (int nlabels, int *vecofint);
int reconstruct_sign (int which, int *gregionsign);
int nextlabel (int label);
int prevlabel (int label);
int inherit (int label);
int isinpath (int i, int start, int end);
int propagate (int sign, int label);
void display_arcs_from_arcs (struct sketch *s);
void display_arcs_from_nodes (struct sketch *s);
void display_regions (struct sketch *s);
void display_regions_from_arcs (struct sketch *s);
void display_regions_from_nodes (struct sketch *s);
int isconsistent (void);
int tour_of_region (int label, int velocity);
void walk_left (int *labelpt, int *velocitypt);

