#include "gtkcontour.h"

static void gtk_add_column(GtkWidget *widget)
{
  if ( ix+45 >= width_max)
  {
    gtk_widget_set_size_request(widget,ix+95,widget->allocation.height);
    width_max=ix+50;
  }
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
  pixmapwid = gtk_pixmap_new( pixmap, mask );
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
static void redraw_brush( GtkWidget *widget)
{
  struct elemento *datoloc;
  struct riga *riga;
  int x,y,tipodato,i,xm;
  PangoLayout * pangolayout;
  gchar *stampa;
  GdkPixmap *pixmaploc;

  stampa=(gchar *) malloc(3 *sizeof(gchar));
  gdk_draw_rectangle (pixmap, widget->style->white_gc, TRUE, 100,0,
                      widget->allocation.width-100,widget->allocation.height);
  riga = primariga;

  tipodato=tipodatoattivo;
  while (riga != NULL)
  {
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
      gtk_widget_queue_draw_area (widget, 155,y-5, width_max,50);

      datoloc = datoloc->datodopo;
    }
    riga = riga->rigadopo;
  }
  tipodatoattivo=tipodato;
  gdk_draw_drawable (widget->window,widget->style->black_gc,pixmap,0,0,0,0,-1,-1);
}

/* Draw a icon on the screen */
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

static void enter_callback( GtkWidget *widget, GtkWidget *entry )
{
  gchar *entry_text;
  entry_text=(gchar *) malloc(12 * sizeof(gchar));
  entry_text = gtk_editable_get_chars (GTK_EDITABLE (entry),0,-1);
  datoattivo->profondita = (gchar *) malloc((strlen(entry_text)+1)*sizeof(gchar));
  sprintf(datoattivo->profondita,",%s",entry_text);
  gtk_widget_destroy(gtk_widget_get_toplevel (widget));
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
    gtk_entry_set_text (GTK_ENTRY (entry), "0");
    gtk_box_pack_end (GTK_BOX (hbox), entry, TRUE, TRUE, 0);
    gtk_widget_show (entry);

    label=gtk_label_new("Inserisci profondita'");
    gtk_box_pack_end (GTK_BOX (hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);
    entry = gtk_entry_new ();

    gtk_widget_show (window);
}

void menuitem_cancellaorientamento(GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto%i.xpm",datoattivo->tipodato);

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=0;
  datoattivo->profondita = "";
}

void menuitem_response(GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto%ir.xpm",datoattivo->tipodato);

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=1;
  richiede_profondita();
}

void menuitem_response1(GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto%il.xpm",datoattivo->tipodato);

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=2;
  richiede_profondita();
}

void menuitem_response2( GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto1u.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=3;
  richiede_profondita();
}

void menuitem_response3( GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto1d.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=4;
  richiede_profondita();
}

void menuitem_response4( GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto2-2.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=5;
}

void menuitem_response5( GtkWidget *wid, GtkWidget *widget)
{
  gchar *xpm_filename;
  xpm_filename=(gchar *) malloc(10 * sizeof(gchar));
  sprintf(xpm_filename,"tasto2-1.xpm");

  pixmapp = gdk_pixmap_create_from_xpm( widget->window, &mask,
                                         &style->bg[GTK_STATE_NORMAL],
                                         xpm_filename);
  draw_brush(wid);

  datoattivo->orientamento=6;
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

   dato->posizionex = dato->posizionex + incremento;
   while (dato->datodopo != NULL)
   {
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


static gint button_press_event( GtkWidget *widget, GdkEventButton *event )
{
  struct elemento *datoloc;
  struct elemento *datoprima;
//  struct elemento *datodopo;
//  struct elemento *datoprimaprima;
  GtkWidget *menu;
  GtkWidget *menu_items;
  char buf[128];
  int x,y;
//  int datidopotrovati;
  struct riga *cercoriga;
    
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
 
  if (cercoriga != NULL)
    datoloc=cercoriga->dato;
  else
    datoloc = NULL;

  if (cercoriga == NULL)
  {
      menu = gtk_menu_new ();
      sprintf(buf,"attenzione posizionarsi su una riga");
      menu_items = gtk_menu_item_new_with_label (buf);
      gtk_menu_append (GTK_MENU (menu), menu_items);
      gtk_widget_show (menu_items);
      gtk_widget_show (menu);
      gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                       event->button, event->time);
  }
  else
  {
  rigaattiva = cercoriga;
  iy=rigaattiva->posizione;

  if (event->button == 1 && tipodatoattivo >= 0)
  { 
    if ( ( tipodatoattivo > 0 &&
      (rigaattiva->rigaprima == NULL || 
      (rigaattiva->rigaprima != NULL && rigaattiva->rigaprima->archiinf == 0) ||
      (rigaattiva->rigaprima != NULL && rigaattiva->rigaprima->archiinf == 1 && tipodatoattivo > 1))) )
    {
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

      if (datoloc == NULL)
      {
        datoloc=(struct elemento *) malloc(sizeof(struct elemento));
        datoloc->datodopo = NULL;
        rigaattiva->dato=datoloc;
        ix = 155;
        datoloc->datosotto1 = NULL;
        datoloc->datosotto2 = NULL;
      }
      else
      {
        while (datoloc->datodopo != NULL && !(x > datoloc->posizionex && x < datoloc->datodopo->posizionex ))
          datoloc=datoloc->datodopo;
        
        if (datoloc->datodopo == NULL) 
        {
// inserisco il dato alla fine
            datoloc->datodopo=(struct elemento *) malloc(sizeof(struct elemento));
            ix = datoloc->posizionex + 50;
            datoloc=datoloc->datodopo;
            datoloc->datodopo=NULL;
        }
        else
        {
// inserisco il nuovo dato nella fila
            datoprima = (struct elemento *) malloc(sizeof(struct elemento));
            datoprima->datodopo = datoloc->datodopo;
            datoloc->datodopo = datoprima;
            ix = datoloc->posizionex +50;
            datoloc = datoprima;
            if (ix == datoloc->datodopo->posizionex) {
              datoprima = datoloc->datodopo;
              while ( datoprima != NULL) {
                datoprima->posizionex = datoprima->posizionex +50;
                datoprima = datoprima->datodopo;
              }
            }

        }
      }

      datoloc->tipodato=tipodatoattivo;
      datoloc->datosotto1=NULL;
      datoloc->datosotto2=NULL;
      datoloc->orientamento=0;
      datoloc->posizionex = ix;

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
 
      redraw_brush (widget);
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
     
     if (tipodatoattivo == -1)
     {

        gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy-50, 16,44);
        gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy+1, 16,44);
        gdk_draw_rectangle (pixmap, widget->style->white_gc,
                      TRUE, 15,iy+50, 16,44);
        gtk_widget_queue_draw_area (widget, 15,iy-50, 25,iy+95);
  
      if (datoprima != NULL)
        datoprima->datodopo = datoloc->datodopo;
      else
        cercoriga->dato = datoloc->datodopo;

      if (cercoriga->rigaprima != NULL)
        datoprima = cercoriga->rigaprima->dato;
      else
        datoprima = NULL;

      while (datoprima != NULL)
      {
        if (datoprima->datosotto1 == datoloc)
        {
          datoprima->datosotto1 = NULL;
          datoprima->archiapertid = datoprima->archiapertid + 1;
        }
        
        if (datoprima->datosotto2 == datoloc)
        {
          datoprima->datosotto2 = NULL;
          datoprima->archiapertid = datoprima->archiapertid + 1;
        }

        datoprima = datoprima->datodopo;
      }
  
      if (datoloc->datosotto1 != NULL)
        datoloc->datosotto1->archiapertiu = datoloc->datosotto1->archiapertiu +1;
 
      if (datoloc->datosotto2 != NULL)
        datoloc->datosotto2->archiapertiu = datoloc->datosotto2->archiapertiu +1;
 
      switch (datoloc->tipodato)
      {
        case 0:
          cercoriga->archiinf=max(cercoriga->archiinf-2,0);
          break;
        case 1:
          cercoriga->archisup=max(cercoriga->archisup-1,0);  
          cercoriga->rigaprima->archiinf=cercoriga->rigaprima->archiinf+1;  
          cercoriga->archiinf=max(cercoriga->archiinf-1,0);  
          break;
        case 2:
          cercoriga->archisup=max(cercoriga->archisup-2,0);  
          cercoriga->rigaprima->archiinf=cercoriga->rigaprima->archiinf+2;  
          cercoriga->archiinf=max(cercoriga->archiinf-2,0);  
        break;
      case 3:
        cercoriga->archisup=max(cercoriga->archisup-2,0);  
        cercoriga->rigaprima->archiinf=cercoriga->rigaprima->archiinf+2;  
        break;
      }

     redraw_brush(widget);
     }
     else
     {
      menu = gtk_menu_new ();
      if (datoloc->tipodato == 0 || datoloc->tipodato == 3)
      {
        sprintf(buf,"orientamento verso destra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response), widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"orientamento verso sinistra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response1), widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"cancella orientamento");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_cancellaorientamento), widget);
        gtk_widget_show (menu_items);
        gtk_widget_show (menu);
        gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                          event->button, event->time);
      }
      else if (datoloc->tipodato == 1)
      {
        sprintf(buf,"orientamento verso l'alto");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response2), widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"orientamento verso il basso");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response3), widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"cancella orientamento");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_cancellaorientamento), widget);
        gtk_widget_show (menu_items);
        gtk_widget_show (menu);
        gtk_menu_popup( GTK_MENU(menu), NULL, NULL, NULL, NULL,
                          event->button, event->time);
      }
      else
      {
        sprintf(buf,"sopra ramo basso-sinistra verso alto-destra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response4),widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"sopra ramo alto-sinistra verso basso-destra");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_response5), widget);
        gtk_widget_show (menu_items);
        sprintf(buf,"cancella orinamento");
        menu_items = gtk_menu_item_new_with_label (buf);
        gtk_menu_append (GTK_MENU (menu), menu_items);
        gtk_signal_connect_object (GTK_OBJECT (menu_items), "activate",
                  GTK_SIGNAL_FUNC (menuitem_cancellaorientamento), widget);
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
static gint expose_event( GtkWidget      *widget, GdkEventExpose *event )
{
  gdk_draw_drawable(widget->window,
		  widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		  pixmap,
		  event->area.x, event->area.y,
		  event->area.x, event->area.y,
		  event->area.width, event->area.height);

  return FALSE;
}

static void gtk_add_drawing_line( GtkWidget *widget, GdkEventExpose *event, GtkWidget *button)
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
  
  id=g_signal_connect(G_OBJECT(widget),"button_press_event",
                      G_CALLBACK(ricavo_posizione),NULL);
}

static void ricavo_posizione(GtkWidget *widget, GdkEventButton *event)
{
  int x,y;
  struct riga *nuovariga;
  struct elemento *dato;
 
  x=(int) event->x;
  y=(int) event->y;
 
  ix=155;
  iy=rigaattiva->posizione+50;
  if (y > iy)
  {
    while (rigaattiva->rigadopo != NULL && y >= iy+50)
    {
      rigaattiva=rigaattiva->rigadopo;
      iy=iy+50;
    }
  }
  else
  {
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
  nuovariga->rigaprima = rigaattiva;
  if (rigaattiva->rigadopo == NULL)
  {
    rigaattiva->rigadopo = nuovariga;
    rigaattiva->rigadopo->rigadopo = NULL;
  }
  else
  {
    nuovariga->rigadopo=rigaattiva->rigadopo;
    rigaattiva->rigadopo = nuovariga;
    dato=rigaattiva->dato;
    while (dato != NULL)
    {
      if (dato->datosotto1 != NULL)
      {
        dato->archiapertid = dato->archiapertid +1;
        dato->datosotto1->archiapertiu = dato->datosotto1->archiapertiu +1;
        rigaattiva->archiinf = rigaattiva->archiinf +1;
      }
      if (dato->datosotto2 != NULL)
      {  
        dato->archiapertid = dato->archiapertid +1;
        dato->datosotto2->archiapertiu = dato->datosotto2->archiapertiu +1;
        rigaattiva->archiinf = rigaattiva->archiinf +1;
      }

      dato->datosotto1 = NULL;
      dato->datosotto2 = NULL;
      dato = dato->datodopo;
    }
  }
  rigaattiva->rigadopo->rigaprima=rigaattiva;
  rigaattiva=rigaattiva->rigadopo;
  rigaattiva->archisup=0;
  rigaattiva->archiinf=0;
  rigaattiva->posizione=iy;
  nuovariga = rigaattiva->rigadopo;
  while (nuovariga != NULL)
  {
    nuovariga->posizione=nuovariga->posizione+50;
    nuovariga = nuovariga->rigadopo;
  }
  
  redraw_brush(widget);
  g_signal_handler_disconnect(G_OBJECT(widget),id);
}

//void file_ok_sel (GtkWidget *w, struct file_gest *file)
void file_ok_sel( struct file_gest *file, GdkEventExpose *event, GtkWidget *w)
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
    redraw_brush (file->salvalegge);
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
  datoloc->datosotto1=NULL;
  datoloc->datosotto2=NULL;
  datoloc->orientamento = 0;
  rigaattiva->archisup = rigaattiva->archisup + archiu;
  rigaattiva->archiinf = rigaattiva->archiinf + archid;
  datoloc->tipodato = tipo;
  datoloc->archiapertiu = archiu;
  datoloc->archiapertid = archid;
  ix = ix + 50;
  datoloc->posizionex  = ix;
  return datoloc;
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

  rigaattiva = primariga;
  while (rigaattiva != NULL)
  {
    datoloc = rigaattiva->dato;
    while (datoloc != NULL )
    {
     if (datoloc->archiapertid != 0 ) 
     {
      datodopo=NULL;
      if ( rigaattiva->rigadopo != NULL && datoloc->tipodato != 3)
      {
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

void leggidati(const gchar *nomefile)
{
  FILE *filein;
//  struct riga *riga;
  struct elemento *datoloc;
//  struct elemento *datodopo;
  gchar valoreletto[80];
//  gchar valorelavoro[5];
  int indice;
  char *err;

  filein=fopen(nomefile,"r");
  err=fgets(valoreletto,80,filein);
  
  while (strncmp(valoreletto,"morse",5) && err != NULL)
    err=fgets(valoreletto,80,filein);

  if (err == NULL)
    g_print("il file selezionato non contiene una descrizione morse \n");
  else
  {
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
    while ((indice=getc(filein)) != '}' )
    {
      if ( indice == ';')
      {
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
      else
      {
        if (indice == '(' || indice == ')' || indice =='\\' || indice == '|' || indice == '/')
          datoloc = alloca_elemento(datoloc,1,1,1);
        else if ( indice == 'X')
          datoloc = alloca_elemento(datoloc,2,2,2);
        else if ( indice == 'U')
          datoloc = alloca_elemento(datoloc,3,2,0);
        else if ( indice == '^')
          datoloc = alloca_elemento(datoloc,0,0,2);
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
          if (*datoloc->profondita != 0) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        else if (datoloc->orientamento == 2)
        {
          write(fdes,"l",1);
          if (*datoloc->profondita != 0) {
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
          if (*datoloc->profondita != 0) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        else if (datoloc->orientamento == 4)
        {
          write(fdes,"d",1);
          if (*datoloc->profondita != 0) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        break;
      case 2:
        write(fdes," X",2);
        break;
      case 3:
        write(fdes," U",2);
        if (datoloc->orientamento == 1)
        {
          write(fdes,"r",1);
          if (*datoloc->profondita != 0) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        else if (datoloc->orientamento == 2)
        {
          write(fdes,"l",1);
          if (*datoloc->profondita != 0) {
            sprintf(buf,"%s",datoloc->profondita);
            write(fdes,buf,strlen(buf));
          }
        }
        break;
      }
      datoloc=datoloc->datodopo;
    } 
    write(fdes,"; \n",3);
    rigaloc=rigaloc->rigadopo;
  }
  write(fdes," } \n",4);
}

static void verify( GtkWidget *wid, GdkEventExpose *event, GtkWidget *button)
{
  int code;
  GtkWidget *menu;
  GtkWidget *menu_items;
  char buf[128];

  code = verifica();
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

int verifica( void )
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

  comando = PATH_CONTOUR;

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

static void saveload( GtkWidget *wid, GdkEventExpose *event, GtkWidget *button)
{
//  GtkWidget *text;
  GtkWidget *filew;
  struct file_gest *file;

  file=(struct file_gest * ) malloc(sizeof(struct file_gest));
  file->salvalegge = wid;

  /* Crea un nuovo widget di selezione file */
  filew = gtk_file_selection_new ("File selection");
  file->filew = filew;

  gtk_signal_connect (GTK_OBJECT (filew), "destroy",
                      (GtkSignalFunc) gtk_widget_destroy, &filew);
  /* Connette ok_button alla funzione file_ok_sel */
  gtk_signal_connect_object(GTK_OBJECT (GTK_FILE_SELECTION (filew)->ok_button),
                      "clicked", (GtkSignalFunc) file_ok_sel, file);
  gtk_signal_connect_object (GTK_OBJECT (GTK_FILE_SELECTION (filew)->ok_button),
                             "clicked", (GtkSignalFunc) gtk_widget_destroy,
                             GTK_OBJECT (filew));

  /* Connette cancel_button alla funzione di distruzione del widget */
  gtk_signal_connect_object (GTK_OBJECT (GTK_FILE_SELECTION (filew)->cancel_button),
                             "clicked", (GtkSignalFunc) gtk_widget_destroy,
                             GTK_OBJECT (filew));

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

  struct riga *riga;
//  struct elemento *elemento;

  gtk_init (&argc, &argv);

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


  hbox = gtk_hbox_new (FALSE, 0);
  gtk_container_add (GTK_CONTAINER (window), hbox);
  gtk_widget_show (hbox);

  vbox = gtk_vbox_new (FALSE, 0);
  gtk_box_pack_start (GTK_BOX (hbox), vbox,TRUE, TRUE, 0);
  gtk_widget_show (vbox);

  gtk_signal_connect (GTK_OBJECT (window), "destroy",
		      GTK_SIGNAL_FUNC (quit), NULL);

  /* Create the drawing area */
  window_scrol = gtk_scrolled_window_new (NULL,NULL);
  gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (window_scrol),
                                    GTK_POLICY_ALWAYS, GTK_POLICY_ALWAYS);

  drawing_area = gtk_drawing_area_new ();
  gtk_widget_set_size_request (drawing_area, 300, 150);
  gtk_widget_show (drawing_area);

  gtk_scrolled_window_add_with_viewport (
                   GTK_SCROLLED_WINDOW (window_scrol), drawing_area);
  gtk_box_pack_start (GTK_BOX (vbox), window_scrol, TRUE, TRUE, 0);
  gtk_widget_show (window_scrol);

  /* Signals used to handle backing pixmap */
  gtk_signal_connect (GTK_OBJECT (drawing_area), "expose_event",
    		      (GtkSignalFunc) expose_event, NULL);
  gtk_signal_connect (GTK_OBJECT(drawing_area),"configure_event",
  		      (GtkSignalFunc) configure_event, NULL);

  /* Event signals */
  gtk_signal_connect_after (GTK_OBJECT (drawing_area), "button_press_event",
		      (GtkSignalFunc) button_press_event, NULL);

  gtk_widget_set_events (drawing_area, GDK_EXPOSURE_MASK
			 | GDK_LEAVE_NOTIFY_MASK
  			 | GDK_BUTTON_PRESS_MASK);

  button = gtk_button_new_with_label ("Aggiungi una riga");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  gtk_signal_connect_object (GTK_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (gtk_add_drawing_line),
			     GTK_OBJECT (drawing_area));
  gtk_widget_show (button);

  button = gtk_button_new_with_label ("Save");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  gtk_signal_connect_object(GTK_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (saveload),
			     NULL);
  gtk_widget_show (button);

  button = gtk_button_new_with_label ("Load");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  gtk_signal_connect_object(GTK_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (saveload),
			     GTK_OBJECT (drawing_area));
  gtk_widget_show (button);

#ifdef HAVE_CONTOUR
  button = gtk_button_new_with_label ("Verify");
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);

  gtk_signal_connect_object(GTK_OBJECT (button), "clicked",
			     GTK_SIGNAL_FUNC (verify),
			     NULL);
  gtk_widget_show (button);
#endif

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
  gtk_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato1);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto1.xpm", "punto di attraversamento");
  gtk_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato2);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto2.xpm", "punto di incrocio");
  gtk_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato3);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto3.xpm", "punto di minimo");
  gtk_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato4);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);

  button = gtk_button_new();
  gtk_container_add (GTK_CONTAINER (bbox), button);
  box1 = xpm_label_box("tasto-1.xpm", "cancella elemento");
  gtk_signal_connect (GTK_OBJECT (button), "clicked",
                        GTK_SIGNAL_FUNC (set_pixmapp_iniz), &tipodato);
  gtk_widget_show(box1);
  gtk_container_add (GTK_CONTAINER (button), box1);
  gtk_widget_show (button);


  gtk_widget_show (window);

  gtk_main ();

  return 0;
}
