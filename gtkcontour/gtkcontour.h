#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtk/gtk.h>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <gdk/gdkpango.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <errno.h>

/* strutture locali */

struct entries
{
  GtkWidget *entry;
  GtkWidget *entry1;
};

struct file_gest
{
  GtkWidget *filew;
  GtkWidget *salvalegge;
};

struct elemento
{
// tipodato=0 massimo; =1 attraversamento; =2 incrocio; =3 minimo);
  int tipodato;
  int archiapertiu;
  int archiapertid;
  struct elemento *datodopo;
  struct elemento *datosotto1;
  struct elemento *datosotto2;
  int orientamento;
  gchar *profondita;
  int posizionex;
};

struct riga
{
  int archisup;
  int archiinf;
  int posizione;
  struct elemento *dato;
  struct riga *rigadopo;
  struct riga *rigaprima;
};

/* variabili globali */

gulong id;
static GdkPixmap *pixmap = NULL;
static GdkPixmap *pixmapp = NULL;
static GdkBitmap *mask;
static GtkStyle *style;
static GtkWidget *pixmapwid;
static gint width_max=300;
static gint ix=155;
static gint iy=55;
static gint tipodatoattivo;
static struct riga *rigaattiva;
static struct riga *primariga;
static struct elemento *datoattivo;

/* strutture di funzioni */

static void gtk_add_column(GtkWidget *);
//static gint set_pixmapp (GtkWidget *, int *);
void set_pixmapp_iniz (GtkWidget *, int *);
void set_pixmapp (GtkWidget *, int *, int *);
void disegnatratti(GtkWidget *, struct elemento *, int);
static void redraw_brush( GtkWidget *);
static void draw_brush( GtkWidget *);
static void enter_callback( GtkWidget *,GtkWidget *);
static void enter_callback2( GtkWidget *,struct entries *);
void richiede_profondita();
void richiede_profondita_2rami(GtkWidget *);
void menuitem_cancellaorientamento(GtkWidget *, GtkWidget *);
void menuitem_response(GtkWidget *, GtkWidget *);
void menuitem_response1(GtkWidget *, GtkWidget *);
void menuitem_response2(GtkWidget *, GtkWidget *);
void menuitem_response3(GtkWidget *, GtkWidget *);
void menuitem_response4(GtkWidget *, GtkWidget *);
void menuitem_response5(GtkWidget *, GtkWidget *);
//void menuitem_response6(GtkWidget *, GtkWidget *);
void aggiorna_posizionex_su_righe_prec(int, struct elemento *);
void aggiorna_posizionex_su_righe_dopo(int, struct elemento *);
void aggiorna_posizionex_su_riga(int, struct elemento *);
int cerco_dati_dopo(struct elemento *, struct elemento *, int);
static gint button_press_event( GtkWidget *, GdkEventButton *);
GtkWidget *xpm_label_box (gchar *, gchar *);
static gint configure_event( GtkWidget *, GdkEventConfigure *);
static gint expose_event( GtkWidget *, GdkEventExpose *);
static void gtk_add_drawing_line(GtkWidget *, GdkEventExpose *, GtkWidget *);
static void ricavo_posizione(GtkWidget *, GdkEventButton *);
void file_ok_sel( struct file_gest *, GdkEventExpose *, GtkWidget *);
struct elemento * alloca_elemento( struct elemento *, int , int , int );
void stampa(struct riga *);
void sistemo_posizione(void);
void leggidati(const gchar *);
void salvadati(int );
char * trovacomando(void);
static void verify( GtkWidget *, GdkEventExpose *, GtkWidget *);
int verifica( void);
static void saveload( GtkWidget *, GdkEventExpose *, GtkWidget *);
void quit ();

static inline int max(int a, int b) {
  return a > b ? a : b;
}
