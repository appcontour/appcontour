#define V_REGULAR 1
#define V_CUSP 2
#define V_CROSS 3
#define V_FIXED 4

struct morseevent {
  int type;
  int ori;
  int ori2;
  int cusps;
  int cusps2;
  struct morseevent *next;
};

struct polyline {
  struct vertex *vertex;
  struct line *line;
  int numvertices;
  double *gradx;
  double *grady;
  struct vertex **specnodes;
  int specnodesnum;
  double h;
};

struct vertex {
  int tag;
  int type;
  double x;
  double y;
  struct vertex *next;
  struct line *line[];
};

struct line {
  int tag;
  int orientation;
  int cusps;
  struct vertex *a;
  struct vertex *b;
  struct line *next;
};

/* prototypes */

char *grinit (int *argcpt, char *argv[]);
void grsetresponsivity (double incrtime);
int grmain (void);
double evolve (struct polyline *contour, double incrtime);
struct line *nextp (struct line *l, struct vertex *p);
struct line *prevp (struct line *l, struct vertex *p);

