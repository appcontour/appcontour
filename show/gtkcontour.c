#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gtk/gtk.h>
#include <strings.h>
#include "grcommon.h"
#include "showcontour.h"

static void menuitem_response (GtkWidget *,gint ,GtkWidget *);
static void menuitem_stop (GtkWidget  *);
static void menuitem_singlestep (GtkWidget  *);
static void gtk_graf_expose (GtkWidget *, GdkPixmap *, double );

extern struct polyline *contour;
static double incrtime = 0.25;
extern int motion;
gboolean  nograf = 1;
extern int steps;
extern char *title;

#define MENU_TOGGLE_MOTION 1
#define MENU_SINGLE_STEP 2
#define MENU_REFINE 3
#define MENU_DEREFINE 4
#define MENU_INC_RESP 10
#define MENU_DEC_RESP 11
#define MENU_QUIT 100

static GdkPixmap *pixmap = NULL;
//GtkWidget *graf;

static GtkItemFactoryEntry menu_items[] = {
  { "/_Menu",         NULL,         NULL, 0, "<Branch>" },
  { "/Menu/Toggle",   "p", menuitem_stop, 1, "<Item>" },
  { "/Menu/Single step",   "s", menuitem_singlestep, 0, "<Item>"},
  { "/Menu/Refine nodes",  NULL , menuitem_response, MENU_REFINE, "<Item>"},
  { "/Menu/Derefine nodes",  NULL , menuitem_response, MENU_DEREFINE, "<Item>"},
  { "/Menu/Increase responsivity",  NULL , menuitem_response, MENU_INC_RESP, "<Item>"},
  { "/Menu/Decrease responsivity",  NULL , menuitem_response, MENU_DEC_RESP, "<Item>"},
  { "/Menu/Quit",     "q", gtk_main_quit, 0, NULL },
};

void get_main_menu( GtkWidget  *window,
                    GtkWidget **menubar,
                    GtkWidget *graf )
{
  GtkItemFactory *item_factory;
  GtkAccelGroup *accel_group;
  gint nmenu_items = sizeof (menu_items) / sizeof (menu_items[0]);

  accel_group = gtk_accel_group_new ();

  /* This function initializes the item factory.
     Param 1: The type of menu - can be GTK_TYPE_MENU_BAR, GTK_TYPE_MENU,
              or GTK_TYPE_OPTION_MENU.
     Param 2: The path of the menu.
     Param 3: A pointer to a gtk_accel_group.  The item factory sets up
              the accelerator table while generating menus.
  */

  item_factory = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<main>",
                                       accel_group);

  /* This function generates the menu items. Pass the item factory,
     the number of items in the array, the array itself, and any
     callback data for the the menu items. */
  gtk_item_factory_create_items (item_factory, nmenu_items, menu_items, graf);

  /* Attach the new accelerator group to the window. */
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);

  if (menubar)
    /* Finally, return the actual menu bar created by the item factory. */
    *menubar = gtk_item_factory_get_widget (item_factory, "<main>");
}

static gint configure_event( GtkWidget         *widget,
                             GdkEventConfigure *event )
{
  if (pixmap)
    gdk_pixmap_unref(pixmap);

  pixmap = gdk_pixmap_new(widget->window,
                          widget->allocation.width,
                          widget->allocation.height,
                          -1);
  gtk_graf_expose(widget, pixmap, 0.);

  return TRUE;
}

static gint expose_event( GtkWidget      *widget,
                          GdkEventExpose *event )
{
  gdk_draw_pixmap(widget->window,
                  widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
                  pixmap,
                  event->area.x, event->area.y,
                  event->area.x, event->area.y,
                  event->area.width, event->area.height);

  return FALSE;
}

static void
gtk_graf_expose (GtkWidget *widget, GdkPixmap *pixmap, double time)
{
 
 
  struct line *line;
  struct vertex *a, *b, *v;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;
  double xb, yb, vx, vy, vnorm;
  char buf[100];
  GdkFont *fixed_font;
  GdkColor color_black;
  GdkColor color_gray1;
  GdkColor color_gray2;
  GdkColor color_white;
  GdkColor color_line;
  GdkColor color_point;
  GdkColormap *colormap;
  fixed_font = gdk_font_load ("-misc-fixed-medium-r-*-*-*-140-*-*-*-*-*-*");

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
  zoomx = 0.7/(maxx - minx);
  zoomy = 0.7/(maxy - miny);
  zoom = zoomx;
  if (zoom > zoomy) zoom = zoomy;


  gdk_draw_rectangle (pixmap,
                      widget->style->white_gc,
                      TRUE,
                      0, 0,
                      widget->allocation.width,
                      widget->allocation.height);

  colormap=gdk_colormap_get_system();
  color_point.red = 0xffff;
  color_point.green = 0;
  color_point.blue = 0;
  if (!gdk_color_alloc(colormap, &color_point)) {
    g_error("couldn't allocate color");
  }
  color_black.red = 0;
  color_black.green = 0;
  color_black.blue = 0;
  if (!gdk_color_alloc(colormap, &color_black)) {
    g_error("couldn't allocate color");
  }
  color_white.red = 65350;
  color_white.green = 65350;
  color_white.blue = 65350;
  if (!gdk_color_alloc(colormap, &color_white)) {
    g_error("couldn't allocate color");
  }
  color_gray1.red=49000;
  color_gray1.green=49000;
  color_gray1.blue=49000;
  if (!gdk_color_alloc(colormap, &color_gray1)) {
    g_error("couldn't allocate color");
  }
  color_gray2.red=32000;
  color_gray2.green=32000;
  color_gray2.blue=32000;
  if (!gdk_color_alloc(colormap, &color_gray2)) {
    g_error("couldn't allocate color");
  }

    for (line = contour->line; line; line = line->next)
    {
      a = line->a;
      b = line->b;

      switch (line->d) {
      case 0:
	color_line=color_black;
        break;
      case 1:
	color_line=color_gray2;
        break;
      case 2:
	color_line=color_gray2;
        break;
      default:
	color_line=color_gray1;
      }

      gdk_gc_set_foreground(widget->style->fg_gc[widget->state],&color_line);
      gdk_gc_set_line_attributes( widget->style->fg_gc[widget->state],2,GDK_LINE_SOLID,GDK_CAP_ROUND,GDK_JOIN_BEVEL);
      gdk_draw_line (pixmap,
                     widget->style->fg_gc[widget->state],
                     widget->allocation.width*(0.5+(a->x - xmed)*zoom) ,
                     widget->allocation.height*(0.5-(a->y - ymed)*zoom) ,
                     widget->allocation.width*(0.5+(b->x - xmed)*zoom) ,
                     widget->allocation.height*(0.5-(b->y - ymed)*zoom) );
      if (a->type == V_CUSP || a->type == V_CROSS)
      {
        gdk_gc_set_foreground(widget->style->fg_gc[widget->state],&color_point);
        gdk_draw_rectangle (pixmap,
                     widget->style->fg_gc[widget->state], 0,
                     widget->allocation.width*(0.5+(a->x - xmed)*zoom)-1.5 ,
                     widget->allocation.height*(0.5-(a->y - ymed)*zoom)-1.5,
		     3,3);
        gdk_gc_set_foreground (widget->style->fg_gc[widget->state],&color_line);
      }
    }

  for (line = contour->line; line; line = line->next)
  {
    if ((line->a->type == V_REGULAR && line->a->line[0]->a->type != V_REGULAR) || line->rarc->loop == line)
    {
      if (line->d > 0)
      {
        snprintf (buf, 70, "%d", line->d);
        xb = (line->b->x + line->a->x)/2.0;
        yb = (line->b->y + line->a->y)/2.0;
        vx = (line->b->x - line->a->x);
        vy = (line->b->y - line->a->y);
        vnorm = sqrt (vx*vx + vy*vy);
        xb -= vy*0.02/vnorm;
        yb += vx*0.02/vnorm;
        gdk_draw_string(pixmap,fixed_font,widget->style->fg_gc[widget->state],
		widget->allocation.width*(0.5+(xb - xmed)*zoom),
		widget->allocation.height*(0.5-(yb - ymed)*zoom), buf);
      }
    }
  }

  snprintf (buf, 98, "showconfig, time=%lf", time);
  gdk_draw_string(pixmap,fixed_font,widget->style->fg_gc[widget->state],5,50,buf);

  return;
}

char
*gtk_grinit (int *argcpt, char *argv[])
{
  static char ident[]="gtk";
//  grparser (argcpt, argv);
  nograf=motion;   
  gtk_init(argcpt, &argv);
  return (ident);
}

gboolean
idle (gpointer data)
{
  GtkWidget *graf =data;
  double time;
  GdkPixmap *local_pixmap;
  local_pixmap=pixmap;

  time = evolve (contour, incrtime);
  gtk_graf_expose(graf,local_pixmap,time);
  pixmap=local_pixmap;
  gdk_draw_pixmap(graf->window,
                  graf->style->fg_gc[GTK_WIDGET_STATE (graf)],
                  pixmap,
                  0,0,
                  0,0,
                  graf->allocation.width,graf->allocation.height);


  steps--;
  if (steps <= 0) {nograf = 0; steps = 10000;}

  return (nograf);
}

int
gtk_grmain ()
{
  GtkWidget *window;
  GtkWidget *graf;
  GtkWidget *vbox;
  GtkWidget *menubar;

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  
  gtk_signal_connect_object (GTK_OBJECT (window), "destroy",
                      GTK_SIGNAL_FUNC (gtk_main_quit), NULL);
  gtk_window_set_title (GTK_WINDOW(window), title);
  gtk_widget_set_usize (GTK_WIDGET(window), 500,500);
  gtk_container_border_width (GTK_CONTAINER (window), 1);
  graf=gtk_drawing_area_new ();

  vbox = gtk_vbox_new (FALSE, 0);
  gtk_container_add (GTK_CONTAINER (window), vbox);
  gtk_widget_show (vbox);

  get_main_menu (window, &menubar, graf);
  gtk_box_pack_start (GTK_BOX (vbox), menubar, FALSE, FALSE, 1);
  gtk_widget_show (menubar);

  if (nograf == 1)
    g_idle_add(idle, graf);

  gtk_box_pack_end (GTK_BOX (vbox), graf, TRUE, TRUE, 1);
   
  gtk_widget_show (graf);
  gtk_signal_connect (GTK_OBJECT(graf),"expose_event",
                      (GtkSignalFunc) expose_event, NULL);
  gtk_signal_connect (GTK_OBJECT(graf),"configure_event",
                      (GtkSignalFunc) configure_event, NULL);
  gtk_widget_show (window);
  gtk_main ();

  return 0;
}

void gtk_toggle_motion(int i )
{
      nograf = 0;
}

static void menuitem_singlestep(GtkWidget  *graf )
{
  steps = 1;
  g_idle_add(idle,graf);
}
static void menuitem_stop( GtkWidget  *graf)
{
    if (nograf)
      nograf = 0;
    else
    {
      nograf = 1;
      g_idle_add(idle,graf);

    } 
}
    

static void menuitem_response(GtkWidget *graf, gint value, GtkWidget *menu_item)
{
  double time;

  switch (value) {

  case MENU_INC_RESP :
    incrtime /= sqrt(2.0);
    break;

  case MENU_DEC_RESP :
    incrtime *= sqrt(2.0);
    break;

  case MENU_DEREFINE :
    contour->h *= sqrt(2.0);
    redistributenodes (contour);
    time = evolve (contour, 0.1);
    redistributenodes (contour);
    break;

  case MENU_REFINE :
    contour->h /= sqrt(2.0);
    redistributenodes (contour);
    time = evolve (contour, 0.1);
    redistributenodes (contour);
    break;
   

  case 666:
    exit (0);
  }
}
