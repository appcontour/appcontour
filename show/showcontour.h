#define V_REGULAR 1
#define V_CUSP 2
#define V_CROSS 3
#define V_FIXED 4

#define EVENT_REDISTRIBUTENODES 1
#define EVENT_REPULSIVEENERGY 2
#define EVENT_KICKOUT 3

#define CHECK_GRADIENT 1

struct timerevent {
  int event;
  double time;
  struct timerevent *next;
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
  double baseh;
  double time;
  double rentime;
  double rgradtime;
  double renergy;
  double *rgradx;
  double *rgrady;
  int rgraddim;
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
  struct vertex *a;
  struct vertex *b;
  struct arc *arc;
  struct line *next;
};

struct arc {
  //int type;
  int orientation;
  int cusps;
  int numsegments;
  int cuspsinserted;
  struct line *first;
  struct line *last;
  struct line *loop;
  struct arc *next;
  int refcount;
};

/* prototypes */

char *grinit (int *argcpt, char *argv[]);
void toggle_motion (int toggle);
void grsetresponsivity (double incrtime);
int grmain (void);
double evolve (struct polyline *contour, double incrtime);
struct line *nextp (struct line *l, struct vertex *p);
struct line *prevp (struct line *l, struct vertex *p);
double normsq (double *x, int dim);
#ifdef CHECK_GRADIENT
void check_gradient (struct polyline *contour);
#endif

