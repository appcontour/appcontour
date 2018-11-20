#define MAXDTCODELEN 200

struct vecofintlist {
  int len;
  int dim;
  struct vecofintlist *next;
  int *handedness;
  int vec[];
};

#define SIZEOFLOIV(n) (sizeof(struct vecofintlist) + (n)*sizeof(int))

/* prototypes */

struct vecofintlist *readdtcode2loiv (FILE *file, int *isgausscodept);
struct vecofintlist *readnakedvecofintlist (FILE *file);
struct vecofintlist *readvecofintlist (FILE *file);
void freeloiv (struct vecofintlist *loiv);
struct sketch *readgausscodeloiv (struct vecofintlist *loiv);
struct sketch *readlinkfromtable (char *linkname);
void chg_underpass (struct vecofintlist *loiv, int nodenum);
struct vecofintlist *read_gausscode_from_string (char *gc);
struct vecofintlist *gausscode_link_to_knot (struct vecofintlist *loiv);
void gausscode2dtcode (struct vecofintlist *loiv, int *vecofint);
struct sketch *orientedgauss2sketch (struct vecofintlist *loiv);
void inherit_gauss2gauss (struct vecofintlist *loiv_knot, struct vecofintlist *loiv_link, int *dt_realization);
struct sketch *realize_dtcode (int numnodes, int *vecofint, int *gregionsign);
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

