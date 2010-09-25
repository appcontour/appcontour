#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>

#define INVALIDINT (-9999)

struct mdesc {
  struct line *lines;
  struct line *lastline;
};

struct line {
  int tag;
  struct event *events;
  struct line *next;
  struct line *prev;
};

struct event {
  int tag;
  int type;
  int orientation;
  int orientation2;
  int huffman;
  int huffman2;
  int cuspsign;
  struct line *line;
  struct event *next;
  struct patch *cusps;
  struct patch *cusps2;
};

struct patch {
  int dvalues;
  int *d;
  struct patchend *start;
  struct patchend *end;
};

struct patchend {
  int isstart;
  struct patch *patch;
  struct patchend *right;
  struct patchend *left;
};

/* --------------------------------------------- */

#define EVENT_VERT 1
#define EVENT_MAX 2
#define EVENT_MIN 3
#define EVENT_CUSP 4
#define EVENT_CROSS 5
#define ORIENT_LEFT 1
#define ORIENT_RIGHT -1
#define ORIENT_UP 1
#define ORIENT_DOWN -1

/* --------------------------------------------- */

void usage (int argc, char *argv[]);
struct mdesc *read_contour (FILE *file);
void write_contour (struct mdesc *);
int get_token (FILE *file);
int get_tokenl (FILE *file);
void unget_token (int);
struct line *get_line (FILE *file);
char *getword (FILE *file, char *word);

void extend_orientations (struct mdesc *);
int check_if_all_oriented (struct mdesc *);
int inherit_orientation_lines (struct line *l);
int count_dangling_down (struct line *l);
int count_dangling_up (struct line *l);

void collect_simple_arcs (struct mdesc *);

void dump (struct mdesc *);
