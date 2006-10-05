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

#define TOP_LENGTH 0.25

struct morsedesc {
  int numrows;
  int rowlimit;
  int maxrowlen;
  int *rowsize;
  int **row;
};

struct polyline {
  struct vertex *vertex;
  struct line *line;
};

struct vertex {
  int tag;
  double x;
  double y;
  struct vertex *next;
};

struct line {
  int tag;
  struct vertex *a;
  struct vertex *b;
  struct line *next;
};

int loadmorse (FILE *file, struct morsedesc *mdesc);
struct polyline *buildpolyline (struct morsedesc *mdesc);
struct vertex *newvertex (struct polyline *contour, double x, double y);
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
  int numrows, i, j;

  mdesc.rowsize = rowsize;
  mdesc.row = row;
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
  line->next = contour->line;
  contour->line = line;
}

struct vertex *
newvertex (struct polyline *contour, double x, double y)
{
  struct vertex *v;

  v = (struct vertex *) malloc (sizeof (struct vertex));
  v->x = x;
  v->y = y;
  v->next = contour->vertex;
  contour->vertex = v;
  return (v);
}

struct polyline *
buildpolyline (struct morsedesc *mdesc)
{
  double dx, dy, x, y;
  struct vertex *danglingnodes[BUFSIZE];
  struct vertex *prevdanglingnodes[BUFSIZE];
  struct vertex *v1, *v2, *v3, *v4, *v5;
  struct line *line;
  struct polyline *contour;
  int numdnodes = 0;
  int i, j, k, prevdangnodes, dangind, prevdangind;

  contour = (struct polyline *) malloc (sizeof (struct polyline));
  contour->vertex = 0;
  contour->line = 0;

  dx = 2.0/(mdesc->maxrowlen + 2);
  dy = 2.0/(mdesc->numrows + 2);

  prevdangnodes = dangind = 0;
  y = 1.0;
  for (i = 0; i < mdesc->numrows; i++)
  {
    y -= dy;
    x = -1.0;
    prevdangind = dangind = 0;
    for (j = 0; j < mdesc->rowsize[i]; j++)
    {
      x += dx;
      switch (mdesc->row[i][j])
      {
        case TYPE_TOP:
          v1 = newvertex (contour, x - TOP_LENGTH*dx, y - TOP_LENGTH*dy);
          v2 = newvertex (contour, x + TOP_LENGTH*dx, y - TOP_LENGTH*dy);
          line = newline (contour, v1, v2);
          danglingnodes[dangind++] = v1;
          danglingnodes[dangind++] = v2;
          break;

        case TYPE_BOT:
          v1 = newvertex (contour, x - TOP_LENGTH*dx, y + TOP_LENGTH*dy);
          v2 = newvertex (contour, x + TOP_LENGTH*dx, y + TOP_LENGTH*dy);
          newline (contour, v1, v2);
          newline (contour, prevdanglingnodes[prevdangind++], v1);
          newline (contour, prevdanglingnodes[prevdangind++], v2);
          break;

        case TYPE_CROSS:
          v1 = newvertex (contour, x, y);
          v2 = newvertex (contour, x - TOP_LENGTH*dx, y + TOP_LENGTH*dy);
          v3 = newvertex (contour, x + TOP_LENGTH*dx, y + TOP_LENGTH*dy);
          v4 = newvertex (contour, x - TOP_LENGTH*dx, y - TOP_LENGTH*dy);
          v5 = newvertex (contour, x + TOP_LENGTH*dx, y - TOP_LENGTH*dy);
          line = newline (contour, v1, v2);
          line = newline (contour, v1, v3);
          line = newline (contour, v1, v4);
          line = newline (contour, v1, v5);
          newline (contour, prevdanglingnodes[prevdangind++], v2);
          newline (contour, prevdanglingnodes[prevdangind++], v3);
          danglingnodes[dangind++] = v4;
          danglingnodes[dangind++] = v5;
          break;

        case TYPE_TRAN:
	  v1 = newvertex (contour, x, y);
          line = newline (contour, prevdanglingnodes[prevdangind++], v1);
          danglingnodes[dangind++] = v1;
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

  printf ("in buildpolyline\n");
  return (contour);
}

int
loadmorse (FILE *file, struct morsedesc *mdesc)
{
  char buf[BUFSIZE];
  int tok, i, j, endinput, numrows, maxrowlen;
  int **row, *rowsize, dim;
  char ch;

  rowsize = mdesc->rowsize;
  row = mdesc->row;
  dim = mdesc->rowlimit;

  tok = gettoken (file);
  if (tok != TOK_MORSE) return (0);
  tok = gettoken (file);
  if (tok != TOK_LBRACE) return (0);

  endinput = 0;
  i = j = 0;
  while ((ch = fgetc (file)) != EOF)
  {
    if (i >= BUFSIZE - 3) exit (1);
    switch (ch)
    {
      case '|':
      case '\\':
      case '/':
      case ')':
      case '(':
        buf[i++] = TYPE_TRAN;
        break;

      case '^': 
        buf[i++] = TYPE_TOP;
        break;

      case 'U':
        buf[i++] = TYPE_BOT;
        break;

      case 'X':
        buf[i++] = TYPE_CROSS;
        break;

      case ';':
        rowsize[j] = i;
        row[j] = (int *) malloc (i * sizeof(int));
        for (i = 0; i < rowsize[j]; i++) row[j][i] = buf[i];
        j++;
        if (j > dim - 3) return (0);
        i = 0;
        break;

      case '}':
        numrows = j;
        break;
    }
    if (endinput) break;
  }
  mdesc->numrows = numrows;
  return (numrows);
}

