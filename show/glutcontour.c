#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "GL/glut.h"
#include "showcontour.h"
#include "xfigexport.h"

extern struct polyline *contour;

static double incrtime = 0.25;
extern int motion;
extern int steps;
extern char *title;
static int isfullscreen = 0;

#define SIZE 0.01

#define MENU_TOGGLE_MOTION 1
#define MENU_SINGLE_STEP 2
#define MENU_REFINE 3
#define MENU_DEREFINE 4
#define MENU_INC_RESP 10
#define MENU_DEC_RESP 11
#define MENU_XFIG_EXPORT 12
#define MENU_QUIT 100

void glut_idle (void);
void keyboard (unsigned char key, int x, int y);

void
display (void)
{
  struct line *line;
  struct vertex *a, *b, *v;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;
  double xb, yb, gray, vx, vy, vnorm;
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
    if ((line->a->type == V_REGULAR && line->a->line[0]->a->type != V_REGULAR) || line->rarc->loop == line) 
    {
      if (line->d > 0)
      {
        snprintf (dbuf, 70, "%d", line->d);
        xb = (line->b->x + line->a->x)/2.0;
        yb = (line->b->y + line->a->y)/2.0;
        vx = (line->b->x - line->a->x);
        vy = (line->b->y - line->a->y);
        vnorm = sqrt (vx*vx + vy*vy);
        xb -= vy*0.02/vnorm;
        yb += vx*0.02/vnorm;
        glRasterPos2d((xb - xmed)*zoom, (yb - ymed)*zoom);
        for (i = 0; i < strlen (dbuf); i++)
          glutBitmapCharacter (GLUT_BITMAP_HELVETICA_12, dbuf[i]);
      }
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
glut_toggle_motion (int toggle)
{
  if (toggle) motion = 1 - motion;
  if (motion) {
    glutIdleFunc (glut_idle);
  } else {
    glutIdleFunc (NULL);
  }
}

void
menu (int value)
{
  double time;
  //FILE *exportfile;

  switch (value) {
  case MENU_TOGGLE_MOTION:
    glut_toggle_motion (1);
    break;

  case MENU_SINGLE_STEP:
    motion = 1;
    steps = 1;
    glut_toggle_motion (0);
    break;

  case MENU_REFINE:
    contour->h /= sqrt(2.0);
    redistributenodes (contour);
    time = evolve (contour, 0.1);
    redistributenodes (contour);
    break;

  case MENU_DEREFINE:
    contour->h *= sqrt(2.0);
    redistributenodes (contour);
    time = evolve (contour, 0.1);
    redistributenodes (contour);
    break;

  case MENU_INC_RESP:
    incrtime /= sqrt(2.0);
    break;

  case MENU_DEC_RESP:
    incrtime *= sqrt(2.0);
    break;

  case MENU_XFIG_EXPORT:
    //exportfile = fopen ("contour.fig", "w");
    xfig_export0 (contour, xfigproblem, &grflags);
    //fclose (exportfile);
    break;

  case MENU_QUIT:
    exit (0);
  }
}

void
visible (int state)
{
  if (state == GLUT_VISIBLE) {
    glut_toggle_motion (0);
  }
}

void
glut_idle (void)
{
  double time;
  char buf[100];

  steps--;
  if (steps <= 0) {glut_toggle_motion(1); steps = 10000;}
  time = evolve (contour, incrtime);
  snprintf (buf, 98, "showcontour, time=%lf", time);
  if (title) glutSetWindowTitle(title); else glutSetWindowTitle(buf);
  glutPostRedisplay();
}

void
mykeyboard (unsigned char key, int x, int y)
{
  static int posx, posy, sizex, sizey;

  switch (key)
  {
    case 'q':
    case 'Q':
      menu (MENU_QUIT);
      break;
    case 'p':
    case 'P':
    case ' ':
      menu (MENU_TOGGLE_MOTION);
      break;
    case 's':
    case 'S':
      menu (MENU_SINGLE_STEP);
      break;
    case 'f':
    case 'F':
      isfullscreen = 1 - isfullscreen;
      if (isfullscreen)
      {
        posx = glutGet (GLUT_WINDOW_X);
        posy = glutGet (GLUT_WINDOW_Y);
        sizex = glutGet (GLUT_WINDOW_WIDTH);
        sizey = glutGet (GLUT_WINDOW_HEIGHT);
        glutFullScreen ();
      } else {
        glutPositionWindow (posx - 4, posy - 20);
        glutReshapeWindow (sizex, sizey);
        /* ugly... it seems that we must take into account
         * the window decoration, which can depend on the
         * window manager used
         */
      }
      break;
    default:
      break;
  }
}

char *
glut_grinit (int *argcpt, char *argv[])
{
  static char ident[]="glut";

  //grparser (argcpt, argv);
  glutInit(argcpt, argv);
  return (ident);
}

void
grsetresponsivity (double lincrtime)
{
  incrtime = lincrtime;
}

int
glut_grmain (void)
{
  glutCreateWindow("showcontour");
  glutDisplayFunc(display);
  glutVisibilityFunc(visible);
  glutCreateMenu (menu);
  glutAddMenuEntry ("Toggle motion (p)", MENU_TOGGLE_MOTION);
  glutAddMenuEntry ("Single step (s)", MENU_SINGLE_STEP);
  glutAddMenuEntry ("Refine nodes", MENU_REFINE);
  glutAddMenuEntry ("Derefine nodes", MENU_DEREFINE);
  glutAddMenuEntry ("Increase responsivity", MENU_INC_RESP);
  glutAddMenuEntry ("Decrease responsivity", MENU_DEC_RESP);
  glutAddMenuEntry ("xfig export", MENU_XFIG_EXPORT);
  glutAddMenuEntry ("Quit (q)", MENU_QUIT);
  glutAttachMenu (GLUT_RIGHT_BUTTON);
  glutKeyboardFunc (mykeyboard);
  //glEnable (GL_POINT_SMOOTH);
  glutMainLoop();
  return (0);
}

