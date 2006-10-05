/*
 * compile with "cc -o glutcontour glutcontour.c parser.o -lglut"
 *
 * requires freeglut-devel package (on Fedora core)
 *
 * usage:  "./glutcontour <example.morse"
 */

#include <stdio.h>
#include <stdlib.h>
#include "parser.h"
#include "GL/freeglut.h"

#define BUFSIZE 1000

#define TYPE_TRAN 1
#define TYPE_TOP 2
#define TYPE_BOT 3
#define TYPE_CROSS 4

#define V_REGULAR 1
#define V_CUSP 2
#define V_CROSS 3

#define TOP_LENGTH 0.25

struct morseevent {
  int type;
  int ori;
  int ori2;
  int cusps;
  int cusps2;
  struct morseevent *next;
};

struct morsedesc {
  int numrows;
  int rowlimit;
  int maxrowlen;
  int *rowsize;
  int **row;
  int **ori;
  int **cusps;
};

struct polyline {
  struct vertex *vertex;
  struct line *line;
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
  struct line *next;
};

int loadmorse (FILE *file, struct morsedesc *mdesc);
struct polyline *buildpolyline (struct morsedesc *mdesc);
struct vertex *newvertex (struct polyline *contour, 
          double x, double y, int type);
struct line *newline (struct polyline *contour, 
                      struct vertex *a, 
                      struct vertex *b);

static struct polyline *contour;

void
display (void)
{
  struct line *line;
  struct vertex *a, *b;

  glClear(GL_COLOR_BUFFER_BIT);
  //glBegin(GL_LINELOOP);  per una poligonale chiusa...
  glBegin(GL_LINES);
    glColor3f(1.0, 1.0, 1.0);  /* white */
    for (line = contour->line; line; line = line->next)
    {
      a = line->a;
      b = line->b;
      glVertex2d(a->x, a->y);
      glVertex2d(b->x, b->y);
    }
  glEnd();
  glFlush();  /* Single buffered, so needs a flush. */
}

int
main (int argc, char *argv[])
{
  struct morsedesc mdesc;
  struct line *line;
  struct vertex *a, *b;
  int rowsize[BUFSIZE];
  int *row[BUFSIZE];
  int *ori[BUFSIZE];
  int *cusps[BUFSIZE];
  int numrows, i, j;

  mdesc.rowsize = rowsize;
  mdesc.row = row;
  mdesc.ori = ori;
  mdesc.cusps = cusps;
  mdesc.rowlimit = BUFSIZE;

  numrows = loadmorse (stdin, &mdesc);
  if (numrows == 0) exit (1);

  mdesc.maxrowlen = 0;
  for (i = 0; i < numrows; i++)
    if (rowsize[i] > mdesc.maxrowlen) mdesc.maxrowlen = rowsize[i];

  //printf ("rows: %d\n", numrows);
  //for (j = 0; j < numrows; j++)
  //{
  //  for (i = 0; i < rowsize[j]; i++)
  //  {
  //    printf ("%d", row[j][i]);
  //  }
  //  printf ("\n");
  //}

  contour = buildpolyline (&mdesc);

  //for (line = contour->line; line; line = line->next)
  //{
  //  a = line->a;
  //  b = line->b;
  //  printf ("line from (%lf, %lf) to (%lf, %lf)\n",
  //          a->x, a->y, b->x, b->y);
  //}
  glutInit(&argc, argv);
  glutCreateWindow("single triangle");
  glutDisplayFunc(display);
  //glutReshapeFunc(reshape);
  glutMainLoop();
}

struct line *
newline (struct polyline *contour, struct vertex *a, struct vertex *b)
{
  struct line *line;

  line = (struct line *) malloc (sizeof (struct line));
  line->a = a;
  line->b = b;
  line->orientation = 0;
  line->next = contour->line;
  contour->line = line;
}

struct vertex *
newvertex (struct polyline *contour, double x, double y, int type)
{
  struct vertex *v;
  int numarcs, i;

  numarcs = 2;
  if (type == V_CROSS) numarcs = 4;
  v = (struct vertex *) malloc (sizeof (struct vertex) 
                         + numarcs*sizeof (struct line *));
  v->x = x;
  v->y = y;
  v->next = contour->vertex;
  v->type = type;
  for (i = 0; i < numarcs; i++) v->line[i] = 0;
  contour->vertex = v;
  return (v);
}

struct polyline *
buildpolyline (struct morsedesc *mdesc)
{
  double dx, dy, x, y, maxx, maxy;
  struct vertex *danglingnodes[BUFSIZE];
  struct vertex *prevdanglingnodes[BUFSIZE];
  struct vertex *v1, *v2, *v3, *v4, *v5, *v;
  struct line *line;
  struct polyline *contour;
  int numdnodes = 0;
  int i, j, k, prevdangnodes, dangind, prevdangind;
  int numrows = 0, numcols = 0;

  contour = (struct polyline *) malloc (sizeof (struct polyline));
  contour->vertex = 0;
  contour->line = 0;

  //dx = 2.0/(mdesc->maxrowlen + 2);
  //dy = 2.0/(mdesc->numrows + 2);
  dx = dy = 1.0;    /* aggiusto alla fine */
  maxx = maxy = 0.0;
  prevdangnodes = dangind = 0;
  y = 0.0;
  for (i = 0; i < mdesc->numrows; i++)
  {
    y += dy;
    x = 0.0;
    numrows++;
    prevdangind = dangind = 0;
    for (j = 0; j < mdesc->rowsize[i]; j++)
    {
      x += dx;
      if (j + 1 > numcols) numcols = j + 1;
      switch (mdesc->row[i][j])
      {
        case TYPE_TOP:
          v1 = newvertex (contour, x - TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
          v2 = newvertex (contour, x + TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
          line = newline (contour, v1, v2);
          v1->line[0] = line;
          v2->line[0] = line;
          danglingnodes[dangind++] = v1;
          danglingnodes[dangind++] = v2;
          break;

        case TYPE_BOT:
          v1 = newvertex (contour, x - TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
          v2 = newvertex (contour, x + TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
          line = newline (contour, v1, v2);
          v1->line[0] = line;
          v2->line[0] = line;
          newline (contour, prevdanglingnodes[prevdangind++], v1);
          newline (contour, prevdanglingnodes[prevdangind++], v2);
          break;

        case TYPE_CROSS:
          v1 = newvertex (contour, x, y, V_CROSS);
          v2 = newvertex (contour, x - TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
          v3 = newvertex (contour, x + TOP_LENGTH*dx, y - TOP_LENGTH*dy,
                V_REGULAR);
          v4 = newvertex (contour, x - TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
          v5 = newvertex (contour, x + TOP_LENGTH*dx, y + TOP_LENGTH*dy,
                V_REGULAR);
          line = newline (contour, v1, v2);
          v1->line[0] = line;
          v2->line[0] = line;
          line = newline (contour, v1, v3);
          v1->line[1] = line;
          v2->line[0] = line;
          line = newline (contour, v1, v4);
          v1->line[2] = line;
          v2->line[0] = line;
          line = newline (contour, v1, v5);
          v1->line[3] = line;
          v2->line[0] = line;
          newline (contour, prevdanglingnodes[prevdangind++], v2);
          newline (contour, prevdanglingnodes[prevdangind++], v3);
          danglingnodes[dangind++] = v4;
          danglingnodes[dangind++] = v5;
          break;

        case TYPE_TRAN:
	  v1 = newvertex (contour, x, y - TOP_LENGTH*dy, V_REGULAR);
	  v2 = newvertex (contour, x, y + TOP_LENGTH*dy, V_REGULAR);
          line = newline (contour, v1, v2);
          v1->line[0] = line;
          v2->line[0] = line;
          line = newline (contour, prevdanglingnodes[prevdangind++], v1);
          danglingnodes[dangind++] = v2;
          break;

        default:
          printf ("caso non gestito: %d\n", mdesc->row[i][j]);
          break;
      }
    }
    if (prevdangnodes != prevdangind)
    {
      printf ("dangling nodes: prev = %d, current = %d\n",
              prevdangnodes, dangind);
    }
    for (k = 0; k < dangind; k++)
    {
      prevdanglingnodes[k] = danglingnodes[k];
    }
    prevdangnodes = dangind;
  }

  dx = 2.0/(mdesc->maxrowlen + 2);
  dy = 2.0/(mdesc->numrows + 2);

  for (v = contour->vertex; v; v = v->next)
  {
    v->x = dx*v->x - 1.0;
    v->y = 1.0 - dy*v->y;
  }

  printf ("in buildpolyline\n");
  return (contour);
}

int
loadmorse (FILE *file, struct morsedesc *mdesc)
{
  char buf[BUFSIZE];
  int tok, i, j, endinput, numrows, maxrowlen;
  int **row, **ori, **cusps, *rowsize, dim;
  char ch;

  rowsize = mdesc->rowsize;
  row = mdesc->row;
  ori = mdesc->ori;
  cusps = mdesc->cusps;
  dim = mdesc->rowlimit;

  tok = gettoken (file);
  if (tok != TOK_MORSE) return (0);
  tok = gettoken (file);
  if (tok != TOK_LBRACE) return (0);

  endinput = 0;
  i = j = 0;
  while ((tok = gettokens (file)) != TOK_RBRACE)
  {
    if (i >= BUFSIZE - 3) exit (1);
    switch (tok)
    {
      case KEY_HAT:
      tok = KEY_A;
      break;

      case KEY_U:
      case KEY_UNDERSCORE:
      tok = KEY_V;
      break;

      case KEY_SLASH:
      case KEY_BSLASH:
      case KEY_BACKQUOTE:
      case TOK_LPAREN:
      case TOK_RPAREN:
      case KEY_PIPE:
      tok = KEY_I;
      break;
    }

    switch (tok)
    {
      case KEY_I:
        buf[i++] = TYPE_TRAN;
        break;

      case KEY_A: 
        buf[i++] = TYPE_TOP;
        break;

      case KEY_V:
        buf[i++] = TYPE_BOT;
        break;

      case KEY_X:
        buf[i++] = TYPE_CROSS;
        break;

      case TOK_SEMICOLON:
        rowsize[j] = i;
        row[j] = (int *) malloc (i * sizeof(int));
        for (i = 0; i < rowsize[j]; i++) row[j][i] = buf[i];
        j++;
        if (j > dim - 3) return (0);
        i = 0;
        break;
    }
  }
  numrows = j;
  mdesc->numrows = numrows;
  return (numrows);
}

