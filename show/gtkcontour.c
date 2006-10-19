#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gtk/gtk.h>
#include <strings.h>
#include "showcontour.h"

static void menuitem_response (gint );
static void menuitem_stop (GtkWidget  *);

extern struct polyline *contour;
static double incrtime = 0.25;
static gboolean  nograf = 1;


static gint
gtk_graf_expose (GtkWidget *widget, double time)
{

  struct line *line;
  struct vertex *a, *b, *v;
  double maxx, maxy, minx, miny, xmed, ymed, zoomx, zoomy, zoom;
  char buf[100];
  GdkFont *fixed_font;
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


  g_return_val_if_fail (widget != NULL, FALSE);

  gdk_window_clear_area (widget->window,
                         0, 0,
                         widget->allocation.width,
                         widget->allocation.height);

    for (line = contour->line; line; line = line->next)
    {
      a = line->a;
      b = line->b;
      gdk_draw_line (widget->window,
                     widget->style->fg_gc[widget->state],
                     widget->allocation.width*(0.5+(a->x - xmed)*zoom) ,
                     widget->allocation.height*(0.5+(a->y - ymed)*zoom) ,
                     widget->allocation.width*(0.5+(b->x - xmed)*zoom) ,
                     widget->allocation.height*(0.5+(b->y - ymed)*zoom) );
    }
  snprintf (buf, 98, "showconfig, time=%lf", time);
  gdk_draw_string(widget->window,fixed_font,widget->style->fg_gc[widget->state],5,50,buf);
  return TRUE;
}

char
*grinit (int *argcpt, char *argv[])
{
  static char ident[]="gtk";
  gtk_init(argcpt, &argv);
  return (ident);
}

gboolean
idle (gpointer data)
{
  GtkWidget *graf =data;
  double time;

  time = evolve (contour, incrtime);
  gtk_graf_expose(graf,time);

  return (nograf);
}

int
grmain ()
{
  GtkWidget *window;
  GtkWidget *graf;
  GtkWidget *vbox;
  int *i,*j;
  int k,kk;
  i=&k;
  j=&kk;

    GtkWidget *menu;
    GtkWidget *menu_bar;
    GtkWidget *root_menu;
    GtkWidget *menu_items;
    char buf[128];

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  
  gtk_signal_connect_object (GTK_OBJECT (window), "destroy",
                      GTK_SIGNAL_FUNC (gtk_main_quit), NULL);
  gtk_window_set_title (GTK_WINDOW(window), "Grafico");
  gtk_widget_set_usize (GTK_WIDGET(window), 500,500);
  gtk_container_border_width (GTK_CONTAINER (window), 1);
  *i=2;
  *j=3;
  graf=gtk_frame_new(NULL);

    menu = gtk_menu_new ();
    sprintf (buf, "Toggle motion");
    menu_items = gtk_menu_item_new_with_label (buf);
    gtk_widget_show (menu_items);
    gtk_menu_append (GTK_MENU (menu), menu_items);
    gtk_signal_connect_object  (GTK_OBJECT (menu_items), "activate",
                 GTK_SIGNAL_FUNC(menuitem_stop), (gpointer) graf);
    sprintf (buf, "Increase responsivity");
    menu_items = gtk_menu_item_new_with_label (buf);
    gtk_widget_show (menu_items);
    gtk_menu_append (GTK_MENU (menu), menu_items);
    gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                 GTK_SIGNAL_FUNC(menuitem_response), (gpointer) j);
    sprintf (buf, "Decrease responsivity");
    menu_items = gtk_menu_item_new_with_label (buf);
    gtk_widget_show (menu_items);
    gtk_menu_append (GTK_MENU (menu), menu_items);
    gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
               GTK_SIGNAL_FUNC(menuitem_response), (gpointer) i);
    sprintf (buf, "Quit");
    menu_items = gtk_menu_item_new_with_label (buf);
    gtk_widget_show (menu_items);
    gtk_menu_append (GTK_MENU (menu), menu_items);
    gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
               GTK_SIGNAL_FUNC(gtk_main_quit),NULL);

    root_menu = gtk_menu_item_new_with_label ("Menu");
    gtk_menu_item_set_submenu (GTK_MENU_ITEM (root_menu), menu);

    vbox = gtk_vbox_new (FALSE, 0);
    gtk_container_add (GTK_CONTAINER (window), vbox);
  
    menu_bar = gtk_menu_bar_new ();
    gtk_menu_bar_append (GTK_MENU_BAR (menu_bar), root_menu);
    gtk_box_pack_start (GTK_BOX (vbox), menu_bar, FALSE, FALSE, 0);
    
  g_idle_add(idle, graf);
  gtk_box_pack_end (GTK_BOX (vbox), graf, TRUE, TRUE, 1);
   
  gtk_widget_show (vbox);
  gtk_widget_show (menu_bar);
  gtk_widget_show (root_menu);
  gtk_widget_show (menu);
  gtk_widget_show (graf);
  gtk_widget_show (window);
  gtk_main ();

  return 0;
}

void toggle_motion(int i )
{
      nograf = 0;
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
    

static void menuitem_response(gint i)
{

  switch (i) {

  case 2 :
    incrtime /= sqrt(2.0);
    break;

  case 3 :
    incrtime *= sqrt(2.0);
    break;

  case 666:
    exit (0);
  }
}
