#include <stdlib.h>
#include <math.h>
#include "GL/glut.h"
#include "showcontour.h"

extern struct polyline *contour;

static double incrtime = 0.25;
static int motion = 1;

void idle (void);

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
  evolve (contour, incrtime);
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
  glutMainLoop();
  return (0);
}

