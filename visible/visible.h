#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>

// #define WITHCONNECTIONS 1

struct mdesc {
  struct line *lines;
};

struct line {
  int tag;
  struct event *events;
  struct line *next;
};

struct event {
  int tag;
  int type;
  int orientation;
  int external_on_right;
  int external_on_left;
  int external_inside;   /* only for max, min, tj */
#ifdef WITHCONNECTIONS
  struct event *a1;
  struct event *a2;
  struct event *a3;
#endif
  struct line *line;
  struct event *next;
};

struct dangarc {
  int d;
  int orientation;
  struct dangarc *next;
};

/*
 * per eventi VERT: a1:up, a2:down
 * per eventi MAX/MIN: a1:left, a2:right
 * per eventi EPTOP: a2:down
 * per eventi EPBOT: a1:up
 * per eventi TJ: a3 e' il tratto che si nasconde
 *             a1 e a2 puntano rispettivamente verso l'alto e verso il basso
 */

/* --------------------------------------------- */

#define EVENT_VERT 1
#define EVENT_MAX 2
#define EVENT_MIN 3
#define EVENT_EPTOP 4
#define EVENT_EPBOT 5
#define EVENT_TJNW 6
#define EVENT_TJNE 7
#define EVENT_TJSW 8
#define EVENT_TJSE 9
#define EVENT_CUSP 10
#define EVENT_XNE 11
#define EVENT_XNW 12
#define EVENT_XSE 13
#define EVENT_XSW 14
#define ORIENT_LEFT 1
#define ORIENT_RIGHT -1
#define ORIENT_UP 1
#define ORIENT_DOWN -1

/* --------------------------------------------- */

int build_contour (struct mdesc *contour);
struct dangarc *skip_invis (struct dangarc *dang);
void print_event (int type, int ori, int d1, int d2,
                  struct dangarc *start, struct dangarc *ndang);
void print_dang_list (struct dangarc *list, struct dangarc *target);
int compute_f_value (struct dangarc *list, struct dangarc *darc);

void usage (int argc, char *argv[]);
struct mdesc *read_vis_contour (FILE *file);
int get_token (FILE *file);
int get_tokenl (FILE *file);
void unget_token (int);
struct line *get_line (FILE *file);
char *getword (FILE *file, char *word);

#ifdef WITHCONNECTIONS
void connect_events_to_each_other (struct mdesc *);
#endif
void extend_orientations (struct mdesc *);
int check_if_all_oriented (struct mdesc *);
void orientation_from_external_region (struct mdesc *contour);
int inherit_orientation_lines (struct line *l);
void extend_external_region (struct mdesc *);
int extend_ext_region_line (struct line *);
int extend_ext_region_lines (struct line *);
int count_dangling_down (struct line *l);
int count_dangling_up (struct line *l);

int check_conditions (struct mdesc *);

void dump (struct mdesc *);
