#include "gtkcontour.h"
#include "parser.h"
  
/*
static void gtk_add_column(GtkWidget *widget)
{
  if ( ix+45 >= width_max)
  {
    gtk_widget_set_size_request(widget,ix+95,widget->allocation.height);
    width_max=ix+50;
  }
}
*/

static inline int max(int a, int b) {
  return a > b ? a : b;
}

void set_pixmapp_iniz (GtkWidget *widget, int * xpm)
{
   int i = 0;
   set_pixmapp (widget, xpm, &i);
}
/* Cambia la pixmapp delle icone */
void set_pixmapp (GtkWidget *widget, int * xpm, int *orientamento)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  switch (*orientamento) {
  case 0:
    sprintf(xpm_filename,"tasto%i.xpm",*xpm);
    break;
  case 1:
    sprintf(xpm_filename,"tasto%ir.xpm",*xpm);
    break;
  case 2:
    sprintf(xpm_filename,"tasto%il.xpm",*xpm);
    break;
  case 3:
    sprintf(xpm_filename,"tasto%iu.xpm",*xpm);
    break;
  case 4:
    sprintf(xpm_filename,"tasto%id.xpm",*xpm);
    break;
  case 5:
    sprintf(xpm_filename,"tasto%i-2.xpm",*xpm);
    break;
  case 6:
    sprintf(xpm_filename,"tasto%i-1.xpm",*xpm);
    break;
  }
 
  style = gtk_widget_get_style( widget );
  pixmapwid = gtk_image_new_from_pixmap( pixmap, mask );
  gtk_widget_show( pixmapwid );
  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);

  tipodatoattivo=*xpm;
}

void disegnatratti(GtkWidget *widget, struct elemento *datoloc, int y)
{
  int lux;
  int ldx;
//  struct elemento *dato;

  if (datoloc->datosotto1 != NULL)
  {
    if (datoloc->tipodato == 1)
      lux = datoloc->posizionex+20;
    else
      lux = datoloc->posizionex;

    if (datoloc->datosotto1->tipodato == 1)
      ldx = datoloc->datosotto1->posizionex+20;
    else
    {
      if (datoloc->tipodato == 1)
      {
    /*    dato = datoloc->datodopo;
        while (dato != NULL && dato->datosotto1 == NULL)
          dato = dato->datodopo;
 
        if (datoloc->datosotto1->tipodato ==1)
          ldx = datoloc->datosotto1->posizionex+20;
        else if(dato != NULL && dato->datosotto1 == datoloc->datosotto1)
          ldx = datoloc->datosotto1->posizionex;
        else
          ldx = datoloc->datosotto1->posizionex+40;
   */ 
       if (datoloc->datosotto1->tipodato ==1)
          ldx = datoloc->datosotto1->posizionex+20;
       else if(datoloc->datosotto1->posizionex > datoloc->posizionex)
          ldx = datoloc->datosotto1->posizionex;
       else
          ldx = datoloc->datosotto1->posizionex+40;
     }
      else
      {
      if ( datoloc->datosotto2 == datoloc->datosotto1 )
        ldx = datoloc->datosotto1->posizionex;
      else
        ldx = datoloc->datosotto1->posizionex+40;
      }
    } 

    gdk_draw_line(pixmap,widget->style->black_gc,lux,y+30,ldx,y+50);
  }

  if (datoloc->datosotto2 != NULL)
  {
    lux = datoloc->posizionex + 40 ;

    if (datoloc->datosotto2->tipodato == 1)
      ldx = datoloc->datosotto2->posizionex+20;
    else
    {
      if (datoloc->datosotto2 == datoloc->datosotto1 )
        ldx = datoloc->datosotto2->posizionex+40;
      else
        ldx = datoloc->datosotto2->posizionex;
    }

    gdk_draw_line(pixmap,widget->style->black_gc,lux,y+30,ldx,y+50);
  }
//  gdk_gc_set_line_attributes(widget->style->black_gc,1,GDK_LINE_SOLID,
//                                  GDK_CAP_NOT_LAST,
//                                  GDK_JOIN_MITER);
}

/* ReDraw */
static void redraw_brush( GtkWidget *widget, GtkWidget *widget1)
{
  struct elemento *datoloc;
  struct riga *riga;
  int x,y,tipodato,i,xm,sommaarchiaperti=0,code=1;
  PangoLayout * pangolayout;
  gchar *stampa;
  GdkPixmap *pixmaploc;

  stampa=(gchar *) malloc(3 *sizeof(gchar));
  gdk_draw_rectangle (pixmap, widget->style->white_gc, TRUE, 0,0,
                      widget->allocation.width,widget->allocation.height);
  for (i=50;i<widget->allocation.height;i=i+50)
    gdk_draw_line(pixmap,widget->style->black_gc,0,i,50,i);

  gdk_draw_line(pixmap,widget->style->black_gc,50,0,50,widget->allocation.height);

  riga = primariga;

  tipodato=tipodatoattivo;
  while (riga != NULL)
  {
    sommaarchiaperti = sommaarchiaperti + riga->archiinf;
    datoloc = riga->dato;
    y = riga->posizione;
    
    if (y >= widget->allocation.height-50 ) {
      widget->allocation.height = widget->allocation.height +50;
      gtk_widget_set_size_request(widget,widget->allocation.width,
                                  widget->allocation.height);
      pixmaploc = gdk_pixmap_new(widget->window, widget->allocation.width,
                          widget->allocation.height, -1);
      gdk_draw_rectangle (pixmaploc, widget->style->white_gc, TRUE, 0, 0,
                      widget->allocation.width,widget->allocation.height);
      gdk_draw_drawable (pixmaploc, widget->style->white_gc, pixmap,
                   0,0,0,0, -1 ,-1);
      for (i=50;i<widget->allocation.height;i=i+50)
        gdk_draw_line(pixmaploc,widget->style->black_gc,0,i,50,i);

      gdk_draw_line(pixmaploc,widget->style->black_gc,50,0,50,widget->allocation.height);

      pixmap = pixmaploc;
    }

    gdk_draw_rectangle (pixmap, widget->style->white_gc, TRUE,
                      15,y+2,16,20);
    gdk_draw_rectangle (pixmap, widget->style->white_gc, TRUE,
                      15,y+24,16,20);
    sprintf(stampa,"%i",riga->archisup);
    pangolayout=gtk_widget_create_pango_layout(widget,stampa);
    gdk_draw_layout(pixmap,widget->style->black_gc,15,y+8, pangolayout);
    sprintf(stampa,"%i",riga->archiinf);
    pangolayout=gtk_widget_create_pango_layout(widget,stampa);
    gdk_draw_layout(pixmap,widget->style->black_gc,15,y+30, pangolayout);

    while (datoloc != NULL) {
      x = datoloc->posizionex;
      if (datoloc->datosotto2 != NULL)  
        xm = datoloc->datosotto2->posizionex+50;
      else
        xm = x;

      if (xm >= widget->allocation.width-50 ) {
        widget->allocation.width = xm +50;
        gtk_widget_set_size_request(widget,widget->allocation.width,
                                  widget->allocation.height);
        width_max=widget->allocation.width;
        pixmaploc = gdk_pixmap_new(widget->window, widget->allocation.width,
                          widget->allocation.height, -1);
        gdk_draw_rectangle (pixmaploc, widget->style->white_gc, TRUE, 0, 0,
                      widget->allocation.width,widget->allocation.height);
        gdk_draw_drawable (pixmaploc, widget->style->white_gc, pixmap,
                   0,0,0,0, -1 ,-1);
        pixmap = pixmaploc;
      }

      set_pixmapp(widget,&datoloc->tipodato,&datoloc->orientamento);
      disegnatratti(widget,datoloc,y);
      gdk_draw_drawable (pixmap, widget->style->black_gc, pixmapp,
                   0,0, x , y, 40,40);
      if (datoloc->profondita != NULL) {
        sprintf(stampa,"%s",datoloc->profondita);
        pangolayout=gtk_widget_create_pango_layout(widget,stampa);
        gdk_draw_layout(pixmap,widget->style->black_gc,x+40,y+10,pangolayout);
      }
  
      gtk_widget_queue_draw_area (widget, 155,y-5, width_max,50);

      datoloc = datoloc->datodopo;
    }
    riga = riga->rigadopo;
  }
  tipodatoattivo=tipodato;

  gdk_draw_drawable (widget->window,widget->style->black_gc,pixmap,0,0,0,0,-1,-1);

#ifdef HAVE_CONTOUR
  if ( pixmapsem == NULL)
    pixmapsem = gdk_pixmap_new(widget1->window, widget1->allocation.width,
                          widget1->allocation.height, -1);

  gdk_draw_rectangle (pixmapsem, widget1->style->bg_gc[4], TRUE, 0,0,
                      widget1->allocation.width,widget1->allocation.height);
  if (sommaarchiaperti != 0)
    pixmapp = gdk_pixmap_create_from_xpm( widget1->window, &mask,&style->bg[GTK_STATE_NORMAL],
                                         "sem_rosso.png");
  else {
    code = verifica(0);
    if (code == 0)
      pixmapp = gdk_pixmap_create_from_xpm( widget1->window, &mask,&style->bg[GTK_STATE_NORMAL],
                                         "sem_verde.png");
    else
      //pixmapp = gdk_pixmap_create_from_xpm( widget1->window, &mask,&style->bg[GTK_STATE_NORMAL],
      pixmapp = gdk_pixmap_create_from_xpm( widget1->window, &mask,&style->bg[0],
                                         "sem_giallo.png");
  }
    
  gdk_draw_drawable (pixmapsem, widget1->style->black_gc, pixmapp,0,0,0,0,50,50);
  
  if (code == 0) {
    code = verifica(1);
    if (code == 0)
      pixmapp = gdk_pixmap_create_from_xpm( widget1->window, &mask,&style->bg[GTK_STATE_NORMAL],
                                         "sem_verde.png");
    else
      pixmapp = gdk_pixmap_create_from_xpm( widget1->window, &mask,&style->bg[GTK_STATE_NORMAL],
                                         "sem_giallo.png");

    gdk_draw_drawable (pixmapsem, widget1->style->black_gc, pixmapp,0,0,0,53,50,50);
  }
  gdk_draw_drawable (widget1->window,widget1->style->black_gc,pixmapsem,0,0,0,0,-1,-1);
#endif

}

/* Draw a icon on the screen */
/*
static void draw_brush( GtkWidget *widget)
{
  int i=0;
  if (pixmapp == NULL)
    set_pixmapp(widget,&i,&i);

  gdk_draw_drawable (pixmap,widget->style->black_gc,
                   pixmapp,0,0,ix,iy,-1,-1);
  gtk_widget_queue_draw_area (widget,ix,iy,-1,-1);
  gdk_draw_drawable (widget->window,widget->style->black_gc,pixmapp,
		   0,0,ix,iy,-1,-1);
  gtk_add_column(widget);
  ix=ix+50;
}
*/

static void enter_callback( GtkWidget *widget, GtkWidget *entry )
{
  gchar *entry_text;
  entry_text=(gchar *) malloc(12 * sizeof(gchar));
  entry_text = gtk_editable_get_chars (GTK_EDITABLE (entry),0,-1);
  if (strlen(entry_text) != 0 ) {
    datoattivo->profondita = (gchar *) malloc((strlen(entry_text)+1)*sizeof(gchar));
    sprintf(datoattivo->profondita,",%s",entry_text);
  }
  else
    datoattivo->profondita = NULL;

  gtk_widget_destroy(gtk_widget_get_toplevel (widget));
  gtk_main_quit();
}

static void enter_callback2( GtkWidget *widget, struct entries *entry )
{
  gchar *entry_text, *entry_text1;
  entry_text=(gchar *) malloc(12 * sizeof(gchar));
  entry_text = gtk_editable_get_chars (GTK_EDITABLE (entry->entry),0,-1);
  datoattivo->profondita = (gchar *) malloc((strlen(entry_text)+1)*sizeof(gchar));
  sprintf(datoattivo->profondita,"%s",entry_text);
  entry_text = gtk_editable_get_chars (GTK_EDITABLE (entry->entry1),0,-1);
  if (strlen(entry_text) != 0 ) {
    entry_text1 = datoattivo->profondita;
    datoattivo->profondita = (gchar *) malloc((strlen(entry_text1)+strlen(entry_text)+3)*sizeof(gchar));
    sprintf(datoattivo->profondita,"[%s]%s",entry_text1,entry_text);
  }
  gtk_widget_destroy(gtk_widget_get_toplevel (widget));
  gtk_main_quit();
}

void richiede_profondita()
{
    GtkWidget *window;
    GtkWidget *hbox;
    GtkWidget *entry;
    GtkWidget *label;

    /* create a new window */
    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_decorated (GTK_WINDOW (window),FALSE);
    gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_MOUSE);

    hbox = gtk_hbox_new (FALSE, 0);
    gtk_container_add (GTK_CONTAINER (window), hbox);
    gtk_widget_show (hbox);

    entry = gtk_entry_new ();
    gtk_entry_set_max_length (GTK_ENTRY (entry), 10);
    g_signal_connect (G_OBJECT (entry), "activate",
                      G_CALLBACK (enter_callback),
                      (gpointer) entry);
//    gtk_entry_set_text (GTK_ENTRY (entry), "0");
    gtk_box_pack_end (GTK_BOX (hbox), entry, TRUE, TRUE, 0);
    gtk_widget_show (entry);

    label=gtk_label_new("Inserisci profondita'");
    gtk_box_pack_end (GTK_BOX (hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);
    entry = gtk_entry_new ();

    gtk_widget_show (window);
    gtk_main();
}


void richiede_profondita_2rami(struct entries *wid)
{
    GtkWidget *window;
    GtkWidget *table;
    GtkWidget *label;
    struct entries *entries_local;

    entries_local = (struct entries *) malloc (sizeof(struct entries));
    /* create a new window */
    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_decorated (GTK_WINDOW (window),FALSE);
    gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_MOUSE);

    table = gtk_table_new(2, 2, TRUE);
    gtk_container_add (GTK_CONTAINER (window), table);
    gtk_widget_show(table);
    label=gtk_label_new("Inserisci profondita' ramo di destra");
    gtk_widget_show(label);
    gtk_table_attach(GTK_TABLE(table), label, 0, 1, 0, 1,
         (GtkAttachOptions) (GTK_EXPAND), (GtkAttachOptions) (0), 0, 0);
    label=gtk_label_new("Inserisci profondita' ramo di sinistra");
    gtk_widget_show(label);
    gtk_table_attach(GTK_TABLE(table), label, 0, 1, 1, 2,
         (GtkAttachOptions) (GTK_EXPAND), (GtkAttachOptions) (0), 0, 0);

    entries_local->entry = gtk_entry_new ();
    gtk_widget_show (entries_local->entry);
    gtk_entry_set_max_length (GTK_ENTRY (entries_local->entry), 10);
    gtk_table_attach(GTK_TABLE(table), entries_local->entry, 1, 2, 0, 1,
         (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 5, 5);
    entries_local->entry1 = gtk_entry_new ();
    gtk_widget_show (entries_local->entry1);
    gtk_entry_set_max_length (GTK_ENTRY (entries_local->entry1), 10);
    gtk_table_attach(GTK_TABLE(table), entries_local->entry1, 1, 2, 1, 2,
         (GtkAttachOptions) (GTK_FILL), (GtkAttachOptions) (0), 5, 5);
    g_signal_connect (G_OBJECT (entries_local->entry), "activate",
                      G_CALLBACK (enter_callback2),
                      (gpointer) entries_local);
    g_signal_connect (G_OBJECT (entries_local->entry1), "activate",
                      G_CALLBACK (enter_callback2),
                      (gpointer) entries_local);

    gtk_widget_show (window);
    gtk_main();
    redraw_brush(wid->entry,wid->entry1);
}

void menuitem_cancellaorientamento(GtkWidget *widget, struct entries *wid)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto%i.xpm",datoattivo->tipodato);

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
//  draw_brush(wid);
  datoattivo->orientamento=0;
  g_free(datoattivo->profondita);
  datoattivo->profondita = NULL;
  redraw_brush(wid->entry,wid->entry1);
}

void menuitem_response(struct entries  *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto%ir.xpm",datoattivo->tipodato);

  pixmapp = gdk_pixmap_create_from_xpm( widget->entry->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  datoattivo->orientamento=1;
  richiede_profondita();
  redraw_brush(widget->entry,widget->entry1);
}

void menuitem_response1(struct entries  *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto%il.xpm",datoattivo->tipodato);

  pixmapp = gdk_pixmap_create_from_xpm( widget->entry->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);

  datoattivo->orientamento=2;
  richiede_profondita();
  redraw_brush(widget->entry,widget->entry1);
}

void menuitem_response2(struct entries *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto1u.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->entry->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);

  datoattivo->orientamento=3;
  richiede_profondita();
  redraw_brush(widget->entry,widget->entry1);
}

void menuitem_response3( struct entries *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto1d.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->entry->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);

  datoattivo->orientamento=4;
  richiede_profondita();
  redraw_brush(widget->entry,widget->entry1);
}

void menuitem_response4(struct entries *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto2-2.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->entry->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);

  datoattivo->orientamento=5;
  redraw_brush(widget->entry,widget->entry1);

}

void menuitem_response5( struct entries *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto2-1.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->entry->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);

  datoattivo->orientamento=6;
  redraw_brush(widget->entry,widget->entry1);
}

void aggiorna_posizionex_su_righe_prec(int posizione, struct elemento *dato)
{
  struct riga *riga;
  int posizioneiniz;
  struct elemento *datoqui; 
  
  riga = rigaattiva->rigaprima;
  posizioneiniz = dato->posizionex;
  if (posizioneiniz >= posizione+50)
    return;

  aggiorna_posizionex_su_riga(posizione,dato);

  while (riga != NULL)
  {
    datoqui = riga->dato;
    while (datoqui != NULL && datoqui->posizionex < posizioneiniz)
      datoqui = datoqui->datodopo;

    if ( datoqui != NULL )
      aggiorna_posizionex_su_riga(posizione+datoqui->posizionex-posizioneiniz,datoqui);
      
    riga = riga->rigaprima;
  }
}

void aggiorna_posizionex_su_righe_dopo(int posizione, struct elemento *dato)
{
  struct riga *riga;
  int posizioneiniz;
  struct elemento *datoqui; 
  
  riga = rigaattiva->rigadopo->rigadopo;
  posizioneiniz = dato->posizionex;
  if (posizioneiniz >= posizione+50)
    return;

  aggiorna_posizionex_su_riga(posizione,dato);

  while (riga != NULL)
  {
    datoqui = riga->dato;
    while (datoqui != NULL && datoqui->posizionex < posizioneiniz)
      datoqui = datoqui->datodopo;

    if ( datoqui != NULL )
      aggiorna_posizionex_su_riga(posizione+datoqui->posizionex-posizioneiniz,datoqui);
      
    riga = riga->rigadopo;
  }
}

void aggiorna_posizionex_su_riga(int posizione, struct elemento *dato)
{
   int incremento;

   incremento = posizione + 50 - dato->posizionex;
 
   if (dato->posizionex_old != 0) 
     dato->posizionex_old = dato->posizionex;

   dato->posizionex = dato->posizionex + incremento;
   while (dato->datodopo != NULL)
   {
     dato->datodopo->posizionex_old = dato->datodopo->posizionex;
     dato->datodopo->posizionex = dato->datodopo->posizionex + incremento;
     dato = dato->datodopo;
   }
}

int cerco_dati_dopo(struct elemento *datodopo, struct elemento *datoloc, int tipodato)
{
  int datidopotrovati=0;
  if (datodopo != NULL) {
    datidopotrovati = 1;
    datodopo->archiapertiu = datodopo->archiapertiu -1;
    datoloc->archiapertid = datoloc->archiapertid -1;
    if ( datoloc->datosotto1 == NULL ) {
      datoloc->datosotto1 = datodopo;
      if ( datoloc->archiapertid  == 1) {
        while (datodopo != NULL && datodopo->archiapertiu == 0)
          datodopo = datodopo->datodopo;

        if ( datodopo != NULL) {
          datodopo->archiapertiu = datodopo->archiapertiu -1;
          datoloc->datosotto2 = datodopo;
          datoloc->archiapertid = datoloc->archiapertid -1;
          datidopotrovati = 2;
        }
        if (datidopotrovati == 2) {
          if ( datoloc->datosotto1 != datoloc->datosotto2) {
            aggiorna_posizionex_su_righe_prec(datoloc->datosotto1->posizionex,
                                            datoloc);
            aggiorna_posizionex_su_righe_dopo(datoloc->posizionex,datoloc->datosotto2);
          } else {
            if (datoloc->posizionex > datoloc->datosotto1->posizionex)
              aggiorna_posizionex_su_righe_dopo(datoloc->posizionex-50,datoloc->datosotto1);
            else
              aggiorna_posizionex_su_righe_prec(datoloc->datosotto1->posizionex-50,
                                            datoloc);
          }
        } else
           aggiorna_posizionex_su_righe_prec(datoloc->datosotto1->posizionex,datoloc);
      } else {
        if (datoloc->datosotto1->tipodato != 1) {
          if (datoloc->datosotto1->archiapertiu == 0 )
            aggiorna_posizionex_su_righe_prec(datoloc->datosotto1->posizionex,
                                            datoloc);
          else
            aggiorna_posizionex_su_righe_dopo(datoloc->posizionex,datoloc->datosotto1);
        } else {
          if (datoloc->posizionex < datoloc->datosotto1->posizionex)
            aggiorna_posizionex_su_righe_prec(datoloc->datosotto1->posizionex-50,
                                          datoloc);
          else 
            aggiorna_posizionex_su_righe_dopo(datoloc->posizionex-50,
                                          datoloc->datosotto1);
        }   
      }
    } else {
      datoloc->datosotto2 = datodopo;
      aggiorna_posizionex_su_righe_dopo(datoloc->posizionex,datoloc->datosotto2);
    }

  }
  return (datidopotrovati);
}

void resetta_posizione_old(void)
{
  struct riga *riga;
  struct elemento *dato;
  
  riga = primariga;
  while (riga != NULL) {
    dato = riga->dato;
    while ( dato != NULL) {
      dato->posizionex_old = dato->posizionex;
      dato = dato->datodopo;
    }
    riga = riga->rigadopo;
  }
}
static gint button_press_event( GtkWidget *widget, GdkEventButton *event , GtkWidget *drawing_area)
{
  struct elemento *datoloc;
  struct elemento *datoprima;
  struct entries *draw;
  GtkWidget *menu;
  GtkWidget *menu_items;
  char buf[128];
  int x,y;
  struct riga *cercoriga;
  draw = (struct entries *) malloc (sizeof(struct entries));

  draw->entry = widget;
  draw->entry1 = drawing_area;

  if (riga_canc_punt != NULL) {
    cancella_puntatori();
    riga_canc_punt = NULL;
  }

  resetta_posizione_old();

  x=(int) event->x;
  y=(int) event->y;
    
  cercoriga=rigaattiva;
  while (cercoriga!= NULL && !((y-cercoriga->posizione)>0 && (y-cercoriga->posizione)<50))
  {
    if (y < cercoriga->posizione)
      cercoriga=cercoriga->rigaprima;
    else
      cercoriga=cercoriga->rigadopo;
  }
 
  if (cercoriga == NULL) {
    datoloc = NULL;
    menu = gtk_menu_new ();
    sprintf(buf,"attenzione posizionarsi su una riga");
    menu_items = gtk_menu_item_new_with_label (buf);
    gtk_menu_append (GTK_MENU (menu), menu_items);
    gtk_widget_show (menu_items);
    gtk_widget_show (menu);
    gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                       event->button, event->time);
  }
  else {
    datoloc=cercoriga->dato;
    rigaattiva = cercoriga;
    iy=rigaattiva->posizione;

    if (event->button == 1 && tipodatoattivo >= 0) { 
    if ( ( tipodatoattivo > 0 &&
      (rigaattiva->rigaprima == NULL || 
      (rigaattiva->rigaprima != NULL && rigaattiva->rigaprima->archiinf == 0) ||
      (rigaattiva->rigaprima != NULL && rigaattiva->rigaprima->archiinf == 1 && tipodatoattivo > 1))) ) {
      menu = gtk_menu_new ();
      sprintf(buf,"attenzione non ci sono abbastanza archi aperti nella riga precedente");
      menu_items = gtk_menu_item_new_with_label (buf);
      gtk_menu_append (GTK_MENU (menu), menu_items);
      gtk_widget_show (menu_items);
      gtk_widget_show (menu);
      gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                       event->button, event->time);
    }
    else
    {
      gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy-50, 16,44);
      gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy+1, 16,44);
      gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy+50, 16,44);
      gtk_widget_queue_draw_area (widget,15,iy-50, 25,iy+95);

      if (datoloc == NULL) {
        datoloc=(struct elemento *) malloc(sizeof(struct elemento));
        datoloc->datodopo = NULL;
        rigaattiva->dato=datoloc;
        ix = 155;
        datoloc->datosotto1 = NULL;
        datoloc->datosotto2 = NULL;
      }
      else {
        if ( x < datoloc->posizionex ) {
// inserisco il dato all'inizio 
          datoprima = (struct elemento *) malloc(sizeof(struct elemento));
          datoprima->datodopo = datoloc;
          datoloc = datoprima;
          cercoriga->dato = datoloc;
          ix = 155;
          if (ix == datoloc->datodopo->posizionex) {
            datoprima = datoloc->datodopo;
            while ( datoprima != NULL) {
              datoprima->posizionex_old = datoprima->posizionex;
              datoprima->posizionex = datoprima->posizionex +50;
              datoprima = datoprima->datodopo;
            }
          }
        }
        else {
          while (datoloc->datodopo != NULL && !(x > datoloc->posizionex && x < datoloc->datodopo->posizionex ))
            datoloc=datoloc->datodopo;
        
          if (datoloc->datodopo == NULL) {
// inserisco il dato alla fine
            datoloc->datodopo=(struct elemento *) malloc(sizeof(struct elemento));
            ix = datoloc->posizionex + 50;
            datoloc=datoloc->datodopo;
            datoloc->datodopo=NULL;
          }
          else {
// inserisco il nuovo dato nella fila
            datoprima = (struct elemento *) malloc(sizeof(struct elemento));
            datoprima->datodopo = datoloc->datodopo;
            datoloc->datodopo = datoprima;
            ix = datoloc->posizionex +50;
            datoloc = datoprima;
            if (ix == datoloc->datodopo->posizionex) {
              datoprima = datoloc->datodopo;
              while ( datoprima != NULL) {
                datoprima->posizionex_old = datoprima->posizionex;
                datoprima->posizionex = datoprima->posizionex +50;
                datoprima = datoprima->datodopo;
              }
            }
          }
        }
      }

      datoloc->tipodato=tipodatoattivo;
      datoloc->datosotto1=NULL;
      datoloc->datosotto2=NULL;
      datoloc->orientamento=0;
      datoloc->profondita = NULL;
      datoloc->posizionex = ix;
      datoloc->posizionex_old = 0;

      switch (tipodatoattivo)
      {
        case 0:
          rigaattiva->archiinf=rigaattiva->archiinf+2;
          datoloc->archiapertid=2;
          datoloc->archiapertiu=0;
          break;
        case 1:
          rigaattiva->archisup=rigaattiva->archisup+1;  
          rigaattiva->archiinf=rigaattiva->archiinf+1;  
          datoloc->archiapertid=1;
          datoloc->archiapertiu=1;
          break;
        case 2:
          rigaattiva->archisup=rigaattiva->archisup+2;  
          rigaattiva->archiinf=rigaattiva->archiinf+2;  
          datoloc->archiapertid=2;
          datoloc->archiapertiu=2;
          break;
        case 3:
          rigaattiva->archisup=rigaattiva->archisup+2;  
          datoloc->archiapertid=0;
          datoloc->archiapertiu=2;
          break;
      }

      cercoriga = rigaattiva;
      sistemo_posizione();
      rigaattiva = cercoriga;
 
      redraw_brush (widget,drawing_area);
    }
  }
  else if (event->button ==3 || tipodatoattivo == -1)
  {
    datoprima = NULL;
    while (datoloc!=NULL && !((x-datoloc->posizionex)>0 && (x-datoloc->posizionex)<50))
    {
      datoprima = datoloc;
      datoloc=datoloc->datodopo;
    }

    if (datoloc == NULL)
    {
      if (tipodatoattivo == -1)
      {
        menu = gtk_menu_new ();
        sprintf(buf,"non ci sono elementi da cancellare");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_widget_show (menu_items);
        gtk_widget_show (menu);
        gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                     event->button, event->time);
      }
    }
    else
    {
     datoattivo = datoloc;
     iy=cercoriga->posizione;
     ix=datoloc->posizionex;
     
     if (tipodatoattivo == -1) {
       gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy-50, 16,44);
       gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy+1, 16,44);
       gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy+50, 16,44);
       gtk_widget_queue_draw_area (widget, 15,iy-50, 25,iy+95);
  
       cancello_elemento(datoloc,datoprima,cercoriga);
       redraw_brush(widget,drawing_area);
     }
     else {
      menu = gtk_menu_new ();
      if (datoloc->tipodato == 0 || datoloc->tipodato == 3) {
        sprintf(buf,"orientamento verso destra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response), draw);
        gtk_widget_show (menu_items);
        sprintf(buf,"orientamento verso sinistra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response1), draw);
        gtk_widget_show (menu_items);
        sprintf(buf,"cancella orientamento");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_cancellaorientamento), draw);
        gtk_widget_show (menu_items);
        gtk_widget_show (menu);
        gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                          event->button, event->time);
      } 
      else if (datoloc->tipodato == 1) {
        sprintf(buf,"orientamento verso l'alto");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response2), draw);
        gtk_widget_show (menu_items);
        sprintf(buf,"orientamento verso il basso");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response3), draw);
        gtk_widget_show (menu_items);
        sprintf(buf,"cancella orientamento");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_cancellaorientamento), draw);
        gtk_widget_show (menu_items);
        gtk_widget_show (menu);
        gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                          event->button, event->time);
      }
      else {
        sprintf(buf,"sopra ramo basso-sinistra verso alto-destra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response4),draw);
        gtk_widget_show (menu_items);
        sprintf(buf,"sopra ramo alto-sinistra verso basso-destra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response5), draw);
        gtk_widget_show (menu_items);
        sprintf(buf,"inserisci profondita' dei rami");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect_swapped (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (richiede_profondita_2rami), draw);
//                  GTK_SIGNAL_FUNC (menuitem_response6), widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"cancella orientamento");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        g_signal_connect (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_cancellaorientamento), draw);
        gtk_widget_show (menu_items);
        gtk_widget_show (menu);
        gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                          event->button, event->time);
      }
     }
    }
    }
  }
  return TRUE;

}

void cancello_elemento(struct elemento *datoloc, struct elemento *datoprima, struct riga *cercoriga)
{
  if (datoprima != NULL)
    datoprima->datodopo = datoloc->datodopo;
  else
    cercoriga->dato = datoloc->datodopo;

  if (cercoriga->rigaprima != NULL)
    datoprima = cercoriga->rigaprima->dato;
  else
    datoprima = NULL;

  while (datoprima != NULL) {
    if (datoprima->datosotto1 == datoloc) {
          datoprima->datosotto1 = NULL;
          datoprima->archiapertid = datoprima->archiapertid + 1;
    }
        
    if (datoprima->datosotto2 == datoloc) {
          datoprima->datosotto2 = NULL;
          datoprima->archiapertid = datoprima->archiapertid + 1;
    }

    datoprima = datoprima->datodopo;
  }
  
  int archiusatid = 0;
  if (datoloc->datosotto1 != NULL){
    datoloc->datosotto1->archiapertiu = datoloc->datosotto1->archiapertiu +1;
    archiusatid++;  
  }
 
  if (datoloc->datosotto2 != NULL){
    datoloc->datosotto2->archiapertiu = datoloc->datosotto2->archiapertiu +1;
    archiusatid++;  
  }
 
  switch (datoloc->tipodato) {
    case 0:
      cercoriga->archiinf=max(cercoriga->archiinf-2+archiusatid,0);
      break;
    case 1:
      cercoriga->archisup=cercoriga->archisup-1;  
      cercoriga->rigaprima->archiinf=cercoriga->rigaprima->archiinf+1;  
      cercoriga->archiinf=max(cercoriga->archiinf-1+archiusatid,0);  
      break;
    case 2:
      cercoriga->archisup=cercoriga->archisup-2;  
      cercoriga->rigaprima->archiinf=cercoriga->rigaprima->archiinf+2;  
      cercoriga->archiinf=max(cercoriga->archiinf-2+archiusatid,0);  
      break;
    case 3:
      cercoriga->archisup=cercoriga->archisup-2;  
      cercoriga->rigaprima->archiinf=cercoriga->rigaprima->archiinf+2;  
      break;
  }
}

GtkWidget *xpm_label_box (gchar *xpm_filename, gchar *label_text)
{
    GtkWidget *box;
    GtkWidget *label;
    GtkWidget *image;

    /* Create box for image and label */
    box = gtk_hbox_new (FALSE, 0);
    gtk_container_set_border_width (GTK_CONTAINER (box), 2);

    /* Now on to the image stuff */
    image = gtk_image_new_from_file (xpm_filename);

    /* Create a label for the button */
    label = gtk_label_new (label_text);

    /* Pack the image and label into the box */
    gtk_box_pack_start (GTK_BOX (box), image, FALSE, FALSE, 3);
    gtk_box_pack_start (GTK_BOX (box), label, FALSE, FALSE, 3);

    gtk_widget_show (image);
    gtk_widget_show (label);

    return box;
}
/* Create a new backing pixmap of the appropriate size */
static gint configure_event( GtkWidget         *widget, GdkEventConfigure *event )
{
  GdkPixmap *pixmaploc;
  PangoLayout * pangolayout;
  gchar *stampa;
  gint i;

  stampa=(gchar *) malloc(sizeof(gchar));
  pixmaploc = gdk_pixmap_new(widget->window, widget->allocation.width,
			  widget->allocation.height, -1);
  gdk_draw_rectangle (pixmaploc, widget->style->white_gc, TRUE, 0, 0,
  		      widget->allocation.width, widget->allocation.height);

  for (i=50;i<widget->allocation.height;i=i+50)
    gdk_draw_line(pixmaploc,widget->style->black_gc,0,i,50,i);

  gdk_draw_line(pixmaploc,widget->style->black_gc,50,0,50,widget->allocation.height);

  if (pixmap)
  {
    gdk_draw_drawable (pixmaploc, widget->style->white_gc, pixmap,
		   0 , 0, 0 , 0, -1 ,-1);
  }
      
  sprintf(stampa,"%i",rigaattiva->archisup);
  pangolayout=gtk_widget_create_pango_layout(widget,stampa);
  gdk_draw_layout(pixmaploc,widget->style->black_gc,15,iy+8,pangolayout);
  sprintf(stampa,"%i",rigaattiva->archiinf);
  pangolayout=gtk_widget_create_pango_layout(widget,stampa);
  gdk_draw_layout(pixmaploc,widget->style->black_gc,15,iy+30,pangolayout);
      
  pixmap=pixmaploc;
  return TRUE;
}

/* Redraw the screen from the backing pixmap */
static gint expose_event( GtkWidget *widget, GdkEventExpose *event)
{
  gdk_draw_drawable(widget->window,
		  widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		  pixmap,
		  event->area.x, event->area.y,
		  event->area.x, event->area.y,
		  event->area.width, event->area.height);

  return FALSE;
}
static gint expose_event1( GtkWidget *widget, GdkEventExpose *event)
{
  if (pixmapsem != NULL )
    gdk_draw_drawable(widget->window,
		  widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		  pixmapsem,
		  event->area.x, event->area.y,
		  event->area.x, event->area.y,
		  event->area.width, event->area.height);

  return FALSE;
}

//static void gtk_add_drawing_line( GtkWidget *widget, GdkEventExpose *event, GtkWidget *button, struct entries *draw)
static void add_line( GtkWidget *button, struct entries *draw)
{ 
//  GtkWidget *menu;
//  menu = gtk_menu_new ();
//  GtkWidget *menu_items;
//  char buf[128];
  
//  menu = gtk_menu_new ();
//  sprintf(buf,"Clicca dove vuoi inserire la nuova riga");
//  menu_items = gtk_menu_item_new_with_label (buf);
//  gtk_menu_append (GTK_MENU (menu), menu_items);
//  gtk_widget_show (menu_items);
//  gtk_widget_show (menu);
//  gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
//                       1,0);
  id = g_signal_connect(G_OBJECT(draw->entry),"button_press_event",
                      G_CALLBACK(ricavo_posizione),draw);
}

static void del_last_element(void)
{
  struct riga *riga;
  struct elemento *dato,*datoprima;

  riga = primariga;
  while (riga != NULL) {
    dato = riga->dato;
    datoprima = NULL;
    while ( dato != NULL) {
      if (dato->posizionex_old == 0 )
        cancello_elemento(dato,datoprima,riga);
       
      dato->posizionex = dato->posizionex_old;
      datoprima = dato;
      dato = dato->datodopo;
    }
    riga = riga->rigadopo;
  }
}

static void del_last_op(GtkWidget *button, struct entries *draw)
{
  struct riga *riga;
  struct elemento *dato;

  if (riga_canc_punt != NULL || (riga_canc_punt == NULL && primariga->dato == NULL)) {

    riga = primariga;
    while (riga != NULL) {
      dato = riga->dato;
      while ( dato != NULL) {
        dato->posizionex_old = dato->posizionex;
        dato = dato->datodopo;
      }
      riga = riga->rigadopo;
    }
    del_add_line();
  }
  else
    del_last_element();
  redraw_brush(draw->entry,draw->entry1);
}
static void del_add_line(void)
{
  if (riga_canc_punt == NULL) { 
//    if (primariga->dato == NULL) {
      riga_canc_punt = primariga;
      riga_canc_punt->posizione = riga_canc_punt->posizione -50;
      primariga = primariga->rigadopo;
    }
//    else
//      return;
//  }
  else {
  riga_canc_punt->rigadopo = riga_canc_punt->rigadopo->rigadopo;
  }

  rigaattiva = riga_canc_punt->rigadopo;
  while (riga_canc_punt->rigadopo != NULL) {
    riga_canc_punt = riga_canc_punt->rigadopo;
    riga_canc_punt->posizione = riga_canc_punt->posizione -50;
  }
  riga_canc_punt = NULL;
}

static void ricavo_posizione( GtkWidget *widget, GdkEventButton *event, struct entries *draw)
{
  int x,y;
  struct riga *nuovariga;
  struct elemento *dato;
 
  x=(int) event->x;
  y=(int) event->y;
 
  ix=155;
  iy=rigaattiva->posizione+50;
  if (y > iy) {
    while (rigaattiva->rigadopo != NULL && y >= iy+50) {
      rigaattiva=rigaattiva->rigadopo;
      iy=iy+50;
    }
  }
  else {
    if (rigaattiva->rigaprima != NULL)
      rigaattiva=rigaattiva->rigaprima;

    while (rigaattiva->rigaprima != NULL && y <= iy-50)
    {
      rigaattiva=rigaattiva->rigaprima;
      iy=iy-50;
    }
    iy=iy-50;
  }
  nuovariga=(struct riga *) malloc(sizeof(struct riga));
  nuovariga->dato = NULL;
  if (primariga == NULL) {
    primariga = nuovariga;
    iy = 55;
    primariga->rigaprima = NULL;
    rigaattiva = primariga;
  }
  else if (primariga->posizione +50 >= y){
    nuovariga->rigadopo = primariga;
    primariga->rigaprima = nuovariga;
    primariga = nuovariga;
    rigaattiva = primariga;
    primariga->rigaprima = NULL;
    iy = 55;
    }else {
      nuovariga->rigaprima = rigaattiva;
      if (rigaattiva->rigadopo == NULL){
        rigaattiva->rigadopo = nuovariga;
        rigaattiva->rigadopo->rigadopo = NULL;
      }
      else{
        nuovariga->rigadopo = rigaattiva->rigadopo;
        rigaattiva->rigadopo->rigaprima = nuovariga;
        rigaattiva->rigadopo = nuovariga;
        riga_canc_punt = rigaattiva;
    } 
    rigaattiva=rigaattiva->rigadopo;
  }
  rigaattiva->archisup=0;
  rigaattiva->archiinf=0;
  rigaattiva->posizione=iy;
  nuovariga = rigaattiva->rigadopo;
  while (nuovariga != NULL)
  {
    nuovariga->posizione=nuovariga->posizione+50;
    nuovariga = nuovariga->rigadopo;
  }
  
  redraw_brush(draw->entry,draw->entry1);
  g_signal_handler_disconnect(G_OBJECT(widget),id);
}

//
// cancella tutti i puntatori della riga indicata
//
void cancella_puntatori()
{ 
  struct elemento *dato;

  dato=riga_canc_punt->dato;
  while (dato != NULL) {
    if (dato->datosotto1 != NULL) {
      dato->archiapertid = dato->archiapertid +1;
      dato->datosotto1->archiapertiu = dato->datosotto1->archiapertiu +1;
      riga_canc_punt->archiinf = riga_canc_punt->archiinf +1;
    }
    if (dato->datosotto2 != NULL) {  
      dato->archiapertid = dato->archiapertid +1;
      dato->datosotto2->archiapertiu = dato->datosotto2->archiapertiu +1;
      riga_canc_punt->archiinf = riga_canc_punt->archiinf +1;
    }

    dato->datosotto1 = NULL;
    dato->datosotto2 = NULL;
    dato = dato->datodopo;
  } 
}

void file_ok_sel (GtkWidget *w, struct file_gest *file)
//void file_ok_sel( struct file_gest *file, GdkEventExpose *event, GtkWidget *w)
{
  const gchar *nomefile;
  GtkWidget *fs;
  FILE *fileout;
  int fdes;

  fs= file->filew;

  nomefile=gtk_file_selection_get_filename (GTK_FILE_SELECTION (fs));

  if (file->salvalegge == NULL) {
    fileout = fopen(nomefile,"w");
    fdes=fileno(fileout);
    salvadati(fdes);
    fclose(fileout);
  }
  else
  {
    leggidati(nomefile);
    redraw_brush (file->salvalegge->entry,file->salvalegge->entry1);
  }
}

struct elemento * alloca_elemento( struct elemento *datoloc, int tipo, int archiu, int archid)
{

  if (datoloc == NULL)
  {
    rigaattiva->dato = (struct elemento *) malloc (sizeof(struct elemento));
    datoloc = rigaattiva->dato;
  }
  else
  {
    datoloc->datodopo = (struct elemento *) malloc (sizeof(struct elemento));
    datoloc = datoloc->datodopo;
  }
  datoloc->datodopo = NULL;
  datoloc->datosotto1 = NULL;
  datoloc->datosotto2 = NULL;
  datoloc->orientamento = 0;
  datoloc->profondita = NULL;
  rigaattiva->archisup = rigaattiva->archisup + archiu;
  rigaattiva->archiinf = rigaattiva->archiinf + archid;
  datoloc->tipodato = tipo;
  datoloc->archiapertiu = archiu;
  datoloc->archiapertid = archid;
  ix = ix + 50;
  datoloc->posizionex  = ix;
  datoloc->posizionex_old  = 0;
  return (datoloc); 
} 

void stampa(struct riga *riga) 
{
  struct elemento *dato;

  while (riga != NULL)
  {
    dato = riga->dato;
    while (dato != NULL)
    {
      g_print("dato tipo %i posizione = % i \n",dato->tipodato,dato->posizionex);
      dato = dato->datodopo;
    }
    g_print("nuova riga\n");
    riga = riga->rigadopo;
  }
  g_print("\n");
  g_print("fine elementi \n");
  g_print("\n");
}

void sistemo_posizione(void)
{
  struct elemento * datoloc;
  struct elemento * datodopo;
//  int archichiusi;
  
  rigaattiva = primariga;
  while (rigaattiva != NULL) {
    datoloc = rigaattiva->dato;
    while (datoloc != NULL ) {
     if (datoloc->archiapertid != 0 ) {
      datodopo=NULL;
      if ( rigaattiva->rigadopo != NULL && datoloc->tipodato != 3) {
        datodopo=rigaattiva->rigadopo->dato;
        while (datodopo != NULL && datodopo->archiapertiu == 0)
          datodopo = datodopo->datodopo;
      } 

      rigaattiva->archiinf=rigaattiva->archiinf - cerco_dati_dopo(datodopo,datoloc,datoloc->tipodato);
     }
     datoloc=datoloc->datodopo;
    }
    rigaattiva = rigaattiva->rigadopo;
  }
}

int getarcinfo (FILE *file,gchar *buffer)
{
  int tok, orientamento;
  int require_rbr = 1;
  gchar *buffpt;

  orientamento = 0;
  buffpt = buffer;

  tok = gettokens (file);
  if (tok == ISNUMBER || tok == KEY_CUSP || tok == KEY_LEFT ||
      tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN || tok == TOK_COMMA)
  {
    ungettoken (tok);
    tok = TOK_LBRACKET;
    require_rbr = 0;
  }
  if (tok != TOK_LBRACKET)
  {
    ungettoken (tok);
    return(0);
  }
  tok = gettokens (file);
  if (tok == TOK_RBRACKET) {
    return(200);
  }

  switch (tok) {
    case KEY_LEFT:
      orientamento = 2; 
      tok = gettokens (file);
      break;
    case KEY_RIGHT:
      orientamento = 1; 
      tok = gettokens (file);
      break;
    case KEY_UP:
      orientamento = 3; 
      tok = gettokens (file);
      break;
    case KEY_DOWN:
      orientamento = 4; 
      tok = gettokens (file);
      break;
  }
   
  if (tok == TOK_COMMA || tok == ISNUMBER || tok == KEY_CUSP)
  {
    if (tok == ISNUMBER || tok == KEY_CUSP) ungettoken (tok);
    while ((tok = gettokens (file)) == ISNUMBER ||
            tok == TOK_PLUS || tok == TOK_MINUS || tok == KEY_CUSP)
    {
      switch (tok) {
        case ISNUMBER:
          *buffpt++ = gettokenchar();
        break;
        case TOK_PLUS:
          *buffpt++ = '+';
        break;
        case TOK_MINUS:
          *buffpt++ = '-';
        break;
       case KEY_CUSP:
          *buffpt++ = 'c';
        break;
      }
    }
  }
  *buffpt = '\0';
  if (require_rbr == 0)
  {
    ungettoken (tok);
    tok = TOK_RBRACKET;
  }
  if (tok != TOK_RBRACKET)
  {
    fprintf (stderr, "Error: right paren expected: %d\n", tok);
  }
  return (orientamento);
}

void leggidati(const gchar *nomefile)
{
  FILE *filein;
//  struct riga *riga;
  struct elemento *datoloc;
  gchar *buffer, *supp;
  int tok;

  buffer = (gchar *) malloc( 100 * sizeof (gchar));
  filein=fopen(nomefile,"r");
  tok = gettoken(filein);
 
  if ( tok != TOK_MORSE) { 
    g_print("il file selezionato non contiene una descrizione morse \n");
    fclose(filein);
    return;
  }

  tok = gettoken(filein);

  if ( tok != TOK_LBRACE){
    g_print("il file selezionato non contiene una descrizione morse \n");
    fclose(filein);
    return;
  }
  
    free(primariga);

    ix=155;
    iy=55;
    rigaattiva=(struct riga * ) malloc(sizeof(struct riga));
    rigaattiva->archisup=0;
    rigaattiva->archiinf=0;
    rigaattiva->dato=NULL;
    rigaattiva->rigaprima=NULL;
    rigaattiva->rigadopo=NULL;
    rigaattiva->posizione=iy;
    primariga=rigaattiva;
    datoloc = NULL;
    while ((tok=gettokens(filein)) != TOK_RBRACE ) {
      buffer[0] = '\0';
      if ( tok == TOK_SEMICOLON) {
        rigaattiva->rigadopo=(struct riga * ) malloc(sizeof(struct riga));
        rigaattiva->rigadopo->rigaprima = rigaattiva;
        rigaattiva = rigaattiva->rigadopo;
        iy=iy+50;
        ix=155;
        rigaattiva->archisup=0;
        rigaattiva->archiinf=0;
        rigaattiva->dato=NULL;
        rigaattiva->rigadopo=NULL;
        rigaattiva->posizione=iy;
        datoloc = NULL;
      }
      else {
        switch (tok)
        {
          case TOK_LPAREN:
          case TOK_RPAREN:
          case KEY_I:
          case KEY_SLASH:
          case KEY_PIPE:
          case KEY_BSLASH:
            datoloc = alloca_elemento(datoloc,1,1,1);
            datoloc->orientamento = getarcinfo (filein,buffer);
            if (strlen(buffer) != 0){
              datoloc->profondita = (gchar *) malloc(strlen(buffer) * sizeof(gchar));
              sprintf(datoloc->profondita,"%s",buffer);
            }
            break;
          case KEY_X:
            datoloc = alloca_elemento(datoloc,2,2,2);
            tok = gettokens (filein);
            switch (tok){
              case KEY_NWSE:
               datoloc->orientamento = 6;
               break;
              case KEY_NESW:
               datoloc->orientamento = 5;
               break;
              default :
                ungettoken(tok);
               break;
            }
            if ((getarcinfo (filein,buffer) >= 200) || strlen(buffer) != 0) {
              datoloc->profondita = (gchar *) malloc((strlen(buffer)+2) * sizeof(gchar));
              sprintf(datoloc->profondita,"[%s]",buffer);
            }
            buffer[0] = '\0';
            getarcinfo (filein,buffer);
            if (strlen(buffer) != 0) {
              supp = datoloc->profondita;
              datoloc->profondita = (gchar *) malloc((strlen(buffer) + strlen(supp) +2) * sizeof(gchar));
              sprintf(datoloc->profondita,"%s[%s]",supp,buffer);
            }
            break;
          case KEY_U:
          case KEY_V:
            datoloc = alloca_elemento(datoloc,3,2,0);
            datoloc->orientamento = getarcinfo (filein,buffer);
            if (strlen(buffer) != 0){
              datoloc->profondita = (gchar *) malloc(strlen(buffer) * sizeof(gchar));
              sprintf(datoloc->profondita,"%s",buffer);
            }
            break;
          case KEY_A:
          case KEY_HAT:
            datoloc = alloca_elemento(datoloc,0,0,2);
            datoloc->orientamento = getarcinfo (filein,buffer);
            if (strlen(buffer) != 0){
              datoloc->profondita = (gchar *) malloc(strlen(buffer) * sizeof(gchar));
              sprintf(datoloc->profondita,"%s",buffer);
            }
            break;
          default:
            printf("token ignorato %d \n", tok);
         } 
      }
   }
  fclose(filein);
  if (rigaattiva->dato == NULL)
    rigaattiva->rigaprima->rigadopo = NULL;

  sistemo_posizione();
  rigaattiva = primariga;
}

void salvadati(int fdes)
{
  struct riga *rigaloc;
  struct elemento *datoloc;
  gchar *buf;

  buf = (gchar *) malloc(10 * sizeof(gchar));

  sprintf(buf,"morse { \n");
  write(fdes,buf,strlen(buf));

  rigaloc=rigaattiva;
  while (rigaloc->rigaprima != NULL)
    rigaloc=rigaloc->rigaprima;
  while (rigaloc != NULL)
  {
    datoloc=rigaloc->dato;
    while (datoloc != NULL)
    {
      switch (datoloc->tipodato)
      {
      case 0:
        write(fdes," ^",2);
        if (datoloc->orientamento == 1)
        {
          write(fdes,"r",1);
          if (datoloc->profondita != NULL ) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        else if (datoloc->orientamento == 2)
        {
          write(fdes,"l",1);
          if (datoloc->profondita != NULL) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        break;
      case 1:
        write(fdes," |",2);
        if (datoloc->orientamento == 3)
        {
          write(fdes,"u",1);
          if (datoloc->profondita != NULL) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        else if (datoloc->orientamento == 4)
        {
          write(fdes,"d",1);
          if (datoloc->profondita != NULL) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        break;
      case 2:
        write(fdes," X",2);
        if (datoloc->orientamento == 6) 
          write(fdes,"`",1);
        if (datoloc->orientamento == 5)
          write(fdes,"'",1);
        if (datoloc->profondita != NULL) {
          sprintf(buf,"%s",datoloc->profondita);
          write(fdes,buf,strlen(buf));
        }
        break;
      case 3:
        write(fdes," U",2);
        if (datoloc->orientamento == 1)
        {
          write(fdes,"r",1);
          if (datoloc->profondita != NULL) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        else if (datoloc->orientamento == 2)
        {
          write(fdes,"l",1);
          if (datoloc->profondita != NULL) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        break;
      }
      datoloc=datoloc->datodopo;
    } 
    write(fdes," ; \n",4);
    rigaloc=rigaloc->rigadopo;
  }
  write(fdes," } \n",4);
}
/*
static void verify( GtkWidget *wid, GdkEventExpose *event, GtkWidget *button)
{
  int code;
  GtkWidget *menu;
  GtkWidget *menu_items;
  char buf[128];

  code = verifica(0);
  if ( code == 0)
    sprintf(buf,"il grafico rappresenta un contorno apparente");
  else
    sprintf(buf,"il grafico NON rappresenta un contorno apparente");

  menu = gtk_menu_new ();
  menu_items = gtk_menu_item_new_with_label (buf);
  gtk_menu_append (GTK_MENU (menu), menu_items);
  gtk_widget_show (menu_items);
  gtk_widget_show (menu);
  gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                   0,0);
}
*/

int verifica( int scelta )
{
  pid_t cpid;
  int retcode;
  int status;
  int code;
  int fdes[2];
  char *comando;
//  extern char **environ;
  extern int errno;
//  int i=0;

  if ( scelta == 0 )
    comando = PATH_CONTOUR;
  else
    comando = PATH_CONTOUR_1;

  retcode = pipe (fdes);
  cpid = fork ();
  if (cpid < 0) exit (1000);

  if (cpid == 0)
  {
      close (fdes[1]);
      dup2 (fdes[0], 0);
      retcode = system(comando);
      exit (WEXITSTATUS (retcode));
  } else {
    close (fdes[0]);
    salvadati(fdes[1]);
    close (fdes[1]);
    waitpid (cpid, &status, 0);
    code = WEXITSTATUS (status);
  }
return (code);
}

static void saveload( GtkWidget *button, struct entries *wid)
{
//  GtkWidget *text;
  GtkWidget *filew;
  struct file_gest *file;

  file=(struct file_gest * ) malloc(sizeof(struct file_gest));
  file->salvalegge = wid;

  /* Crea un nuovo widget di selezione file */
  filew = gtk_file_selection_new ("File selection");
  file->filew = filew;

  g_signal_connect (G_OBJECT (filew), "destroy",
                    G_CALLBACK (gtk_widget_destroy), G_OBJECT (filew));
  /* Connette ok_button alla funzione file_ok_sel */
  g_signal_connect(G_OBJECT (GTK_FILE_SELECTION (filew)->ok_button),
                      "clicked", G_CALLBACK(file_ok_sel), (gpointer) file);
  g_signal_connect_swapped(G_OBJECT (GTK_FILE_SELECTION (filew)->ok_button),
                    "clicked", G_CALLBACK(gtk_widget_destroy), G_OBJECT (filew));

  /* Connette cancel_button alla funzione di distruzione del widget */
  g_signal_connect_swapped (G_OBJECT (GTK_FILE_SELECTION (filew)->cancel_button),
                             "clicked", G_CALLBACK(gtk_widget_destroy),
                             G_OBJECT (filew));

  /* Preassegnamo un nome di file, come se stessimo dando un valore per difetto in
  dialogo di tipo `` salva con nome '' */
  gtk_file_selection_set_filename (GTK_FILE_SELECTION(filew),
                                   "morse.morse");
  gtk_widget_show(filew);
}

void quit ()
{
  gtk_exit (0);
}

int main( int argc, char *argv[] )
{
  GtkWidget *window;
  GtkWidget *window_scrol;
  GtkWidget *drawing_area;
  GtkWidget *drawing1_area;
  GtkWidget *separator;
  GtkWidget *vbox,*hbox;
//  GtkObject *adjustment;
  GtkWidget *bbox;

  GtkWidget *button,*box1;
  int tipodato;
  int tipodato1;
  int tipodato2;
  int tipodato3;
  int tipodato4;
  struct entries *draw;

  struct riga *riga;
//  struct elemento *elemento;

  gtk_init (&argc, &argv);

  draw = (struct entries *) malloc (sizeof(struct entries));
  riga=(struct riga * ) malloc(sizeof(struct riga));
  riga->archisup=0;
  riga->archiinf=0;
  riga->dato=NULL;
  riga->rigaprima=NULL;
  riga->rigadopo=NULL;
  riga->posizione=iy;
  rigaattiva=riga;
  primariga = riga;
  
  tipodato=-1;
  tipodato1=0;
  tipodato2=1;
  tipodato3=2;
  tipodato4=3;

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_default_size(GTK_WINDOW(window),900,300);
  gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
  gtk_signal_connect (GTK_OBJECT (window), "destroy",
		      GTK_SIGNAL_FUNC (quit), NULL);

  style = gtk_widget_get_style( window );

  hbox = gtk_hbox_new (FALSE, 0);
  gtk_container_add (GTK_CONTAINER (window), hbox);
  gtk_widget_show (hbox);

  vbox = gtk_vbox_new (FALSE, 0);
  gtk_box_pack_start (GTK_BOX (hbox), vbox,TRUE, TRUE, 0);
  gtk_widget_show (vbox);

  /* Create the drawing areas */
  hbox = gtk_hbox_new (FALSE, 0);
  gtk_box_pack_start (GTK_BOX (vbox),hbox, TRUE, TRUE, 0);
  gtk_widget_show (hbox);

  drawing1_area = gtk_drawing_area_new ();
  gtk_widget_set_size_request (drawing1_area, 50,  150);
  gtk_widget_show (drawing1_area);
  /* Signals used to handle backing pixmap */
  g_signal_connect (GTK_OBJECT (drawing1_area), "expose_event",
    		      (GtkSignalFunc) expose_event1, NULL);

  gtk_box_pack_start (GTK_BOX (hbox), drawing1_area, FALSE, FALSE, 0);

  window_scrol = gtk_scrolled_window_new (NULL,NULL);
  gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (window_scrol),
                                    GTK_POLICY_ALWAYS, GTK_POLICY_ALWAYS);

  draw->entry1 = drawing1_area;

  drawing_area = gtk_drawing_area_new ();
  gtk_widget_set_size_request (drawing_area, 300, 150);
  gtk_widget_show (drawing_area);

  gtk_scrolled_window_add_with_viewport (
                   GTK_SCROLLED_WINDOW (window_scrol), drawing_area);
  gtk_box_pack_start (GTK_BOX (hbox), window_scrol, TRUE, TRUE, 0);
//  gtk_box_pack_start (GTK_BOX (vbox), window_scrol, TRUE, TRUE, 0);
  gtk_widget_show (window_scrol);

  /* Signals used to handle backing pixmap */
  g_signal_connect (GTK_OBJECT (drawing_area), "expose_event",
    		      (GtkSignalFunc) expose_event, NULL);
  g_signal_connect (GTK_OBJECT(drawing_area),"configure_event",
  		      (GtkSignalFunc) configure_event, NULL);

  /* Event signals */
  g_signal_connect_after (GTK_OBJECT (drawing_area), "button_press_event",
		      (GtkSignalFunc) button_press_event, drawing1_area );

  gtk_widget_set_events (drawing_area, GDK_EXPOSURE_MASK
			 | GDK_LEAVE_NOTIFY_MASK
  			 | GDK_BUTTON_PRESS_MASK);
 
  draw->entry = drawing_area;

  button = gtk_button_new_with_label ("Aggiungi una riga");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  g_signal_connect (G_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (add_line),
			     draw);
  gtk_widget_show (button);

  button = gtk_button_new_with_label ("Annulla ultimo inserimento");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  g_signal_connect (G_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (del_last_op),
			     draw);
  gtk_widget_show (button);

  button = gtk_button_new_with_label ("Save");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  g_signal_connect(G_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (saveload),
			     NULL);
  gtk_widget_show (button);

  button = gtk_button_new_with_label ("Load");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  g_signal_connect(G_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (saveload),
			     draw);
  gtk_widget_show (button);

//#ifdef HAVE_CONTOUR
//  button = gtk_button_new_with_label ("Verify");
//  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
//
//  gtk_signal_connect_object(GTK_OBJECT (button), "clicked",
//			     GTK_SIGNAL_FUNC (verify),
//			     NULL);
//  gtk_widget_show (button);
//#endif

  button = gtk_button_new_with_label ("Quit");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  gtk_signal_connect_object (GTK_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (gtk_widget_destroy),
			     GTK_OBJECT (window));
  gtk_widget_show (button);

  separator = gtk_vseparator_new ();
  gtk_box_pack_start (GTK_BOX (hbox), separator, FALSE, TRUE, 5);
  gtk_widget_show (separator);

  vbox = gtk_vbox_new (FALSE, 0);
  gtk_box_pack_start (GTK_BOX (hbox), vbox,FALSE,FALSE, 0);
  gtk_widget_show (vbox);

  bbox = gtk_vbutton_box_new ();
  gtk_box_pack_start (GTK_BOX (vbox),bbox, FALSE, FALSE, 0);
  gtk_widget_show (bbox);
  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto0.xpm", "punto di massimo");
  g_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato1);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto1.xpm", "punto di attraversamento");
  g_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato2);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto2.xpm", "punto di incrocio");
  g_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato3);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto3.xpm", "punto di minimo");
  g_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato4);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto-1.xpm", "cancella elemento");
  g_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);


  gtk_widget_show (window);

  gtk_main ();

  return 0;
}
