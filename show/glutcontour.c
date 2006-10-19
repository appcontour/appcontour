#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "GL/glut.h"
#include "showcontour.h"

extern struct polyline *contour;

static double incrtime = 0.25;
static int motion = 1;
//static int motion = 0;     /* for now start with no motion */

#define SIZE 0.01

void idle (void);

void
display (void)
{
  struct line *line;
  struct vertex *a, *b, *v;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;
  double xb, yb, gray;
  char dbuf[80];
  int i;

  maxx = -1000.0;
  maxy = -1000.0;
  minx =  1000.0;
  miny =  1000.0;
  for (v = contour->vertex; v; v = v->next)
  {
    if (maxx < v->x) maxx = v->x;
    if (maxy < v->y) maxy = v->y;
    if (minx > v->x) minx = v->x;
    if (miny > v->y) miny = v->y;
  }
  xmed = (maxx + minx)/2.0;
  ymed = (maxy + miny)/2.0;
  zoomx = 1.9/(maxx - minx);
  zoomy = 1.9/(maxy - miny);
  zoom = zoomx;
  if (zoom > zoomy) zoom = zoomy;
  glClear(GL_COLOR_BUFFER_BIT);
  //glBegin(GL_LINELOOP);  per una poligonale chiusa...
  glBegin(GL_LINES);
//    glColor3f(1.0, 1.0, 1.0);  /* white */
    for (line = contour->line; line; line = line->next)
    {
      gray = 1.0/(line->d + 1);
      glColor3f(gray, gray, gray);  /* white */
      a = line->a;
      b = line->b;
      glVertex2d((a->x - xmed)*zoom, (a->y - ymed)*zoom);
      glVertex2d((b->x - xmed)*zoom, (b->y - ymed)*zoom);
      //if (line->a->type != V_REGULAR) printf ("line->d = %d\n", line->d);
      //if (line->arc->loop == line) printf ("loop: %d\n", line->d);
    }
  glEnd();
  for (line = contour->line; line; line = line->next)
  {
    if (line->a->type != V_REGULAR || line->arc->loop == line) 
    {
      snprintf (dbuf, 70, " %d", line->d);
      xb = line->b->x;
      yb = line->b->y;
      glRasterPos2d((xb - xmed)*zoom, (yb - ymed)*zoom);
      for (i = 0; i < strlen (dbuf); i++)
        glutBitmapCharacter (GLUT_BITMAP_HELVETICA_12, dbuf[i]);
    }
  }
//  glGetDoublev (GL_POINT_SIZE, &ptsize);
//  printf ("point size: %lf\n", ptsize);
  glColor3f(1.0, 0.0, 0.0);  /* red */
  for (v = contour->vertex; v; v = v->next)
  {
    switch (v->type)
    {
      case V_CROSS:
      case V_CUSP:
        //glBegin(GL_POLYGON);
        glBegin(GL_LINE_LOOP);
          glVertex2d((v->x - xmed)*zoom + SIZE, (v->y - ymed)*zoom + SIZE);
          glVertex2d((v->x - xmed)*zoom + SIZE, (v->y - ymed)*zoom - SIZE);
          glVertex2d((v->x - xmed)*zoom - SIZE, (v->y - ymed)*zoom - SIZE);
          glVertex2d((v->x - xmed)*zoom - SIZE, (v->y - ymed)*zoom + SIZE);
        glEnd();
        break;

    }
  }
  glFlush();  /* Single buffered, so needs a flush. */
}

void
toggle_motion (int toggle)
{
  if (toggle) motion = 1 - motion;
  if (motion) {
    glutIdleFunc (idle);
  } else {
    glutIdleFunc (NULL);
  }
}

void
menu (int value)
{
  switch (value) {
  case 1:
    toggle_motion (1);
    break;

  case 2:
    incrtime /= sqrt(2.0);
    break;

  case 3:
    incrtime *= sqrt(2.0);
    break;

  case 666:
    exit (0);
  }
}

void
visible (int state)
{
  if (state == GLUT_VISIBLE) {
    toggle_motion (0);
  }
}

void
idle (void)
{
  double time;
  char buf[100];

  time = evolve (contour, incrtime);
  snprintf (buf, 98, "showconfig, time=%lf", time);
  glutSetWindowTitle(buf);
  glutPostRedisplay();
}

char *
grinit (int *argcpt, char *argv[])
{
  static char ident[]="glut";
  int goon = 1;
  int i, j;

  while (goon)
  {
    goon = 0;
    for (i = 1; i < *argcpt; i++)
    {
      if (strcmp (argv[i], "--nomotion") == 0)
      {
        motion = 0;
        goon = 1;
        (*argcpt)--;
        for (j = i; j < *argcpt; j++)
        {
          argv[j] = argv[j+1];
        }
      }
    }
  }

  glutInit(argcpt, argv);
  return (ident);
}

void
grsetresponsivity (double lincrtime)
{
  incrtime = lincrtime;
}

int
grmain (void)
{
  glutCreateWindow("showcontour");
  glutDisplayFunc(display);
  glutVisibilityFunc(visible);
  glutCreateMenu (menu);
  glutAddMenuEntry ("Toggle motion", 1);
  glutAddMenuEntry ("Increase responsivity", 2);
  glutAddMenuEntry ("Decrease responsivity", 3);
  glutAddMenuEntry ("Quit", 666);
  glutAttachMenu (GLUT_RIGHT_BUTTON);
  //glEnable (GL_POINT_SMOOTH);
  glutMainLoop();
  return (0);
}

