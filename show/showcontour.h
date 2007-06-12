#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define GO_NULL 1
#define GO_GTK 2
#define GO_GLUT 3
#define GO_XFIG 4

#define V_REGULAR 1
#define V_CUSP 2
#define V_FIXED 4
/* the following can coexist, check for crossing
 * should be done with "xxx->type & V_CROSS"
 * (bitwise and)
 */
#define V_CROSS 8
#define V_TYPES 15         /* mask type information */
#define V_CROSS_NWSE 16
#define V_CROSS_NESW 32

#define EVENT_REDISTRIBUTENODES 1
#define EVENT_REPULSIVEENERGY 2
#define EVENT_KICKOUT 3

#define MINF (-10000)

//#define CHECK_GRADIENT 1

extern char *xfigproblem;

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

/*
 * please note, in case of CROSSING, we have 4 lines coming out from
 * a node; they are numbered from 0 to 3 in such a way that
 * 3 - i give the opposite line (the one that continues directly
 * across the node)
 */

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
  int d;
  struct vertex *a;
  struct vertex *b;
  struct earc *earc;
  struct rarc *rarc;
  struct line *next;
};

struct earc {
  //int type;
  int orientation;
  int cusps;
  int numsegments;
  int cuspsinserted;
  int *d;
  struct line *first;
  struct line *last;
  struct line *loop;
  //struct arc *next;
  //struct arc *parent;
  int refcount;
};

struct rarc {
  int numsegments;
  int d;
  struct line *first;
  struct line *last;
  struct line *loop;
};

/* prototypes */

char *grinit (int *argcpt, char *argv[]);
#ifdef HAVE_GTK
char *gtk_grinit (int *argcpt, char *argv[]);
#endif
#ifdef HAVE_GLUT
char *glut_grinit (int *argcpt, char *argv[]);
#endif

void toggle_motion (int toggle);
void glut_toggle_motion (int toggle);
void gtk_toggle_motion (int toggle);
void grsetresponsivity (double incrtime);
int grmain (void);
int glut_grmain (void);
int gtk_grmain (void);
double evolve (struct polyline *contour, double incrtime, double incsimtime);
/* incrtime is a realtime increment, incsimtime is an increment in the
 * simulation time
 */
void redistributenodes (struct polyline *contour);
void init_rarc (struct polyline *contour);
struct line *nextp (struct line *l, struct vertex *p);
struct line *prevp (struct line *l, struct vertex *p);
double normsq (double *x, int dim);
void removenode (struct polyline *contour, struct vertex *v);
void removeline (struct polyline *contour, struct line *l);
#ifdef CHECK_GRADIENT
void check_gradient (struct polyline *contour);
#endif

