#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "GL/glut.h"
#include "showcontour.h"

extern struct polyline *contour;

static double incrtime = 0.25;
static int motion = 1;

#define SIZE 0.015

void idle (void);

void
display (void)
{
  struct line *line;
  struct vertex *a, *b, *v;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;

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
    glColor3f(1.0, 1.0, 1.0);  /* white */
    for (line = contour->line; line; line = line->next)
    {
      a = line->a;
      b = line->b;
      glVertex2d((a->x - xmed)*zoom, (a->y - ymed)*zoom);
      glVertex2d((b->x - xmed)*zoom, (b->y - ymed)*zoom);
    }
  glEnd();
//  glGetDoublev (GL_POINT_SIZE, &ptsize);
//  printf ("point size: %lf\n", ptsize);
    glColor3f(1.0, 0.0, 0.0);  /* white */
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
menu (int value)
{
  switch (value) {
  case 1:
    motion = 1 - motion;
    if (motion) {
      glutIdleFunc (idle);
    } else {
      glutIdleFunc (NULL);
    }
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
    if (motion)
      glutIdleFunc (idle);
  } else {
    glutIdleFunc (NULL);
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

