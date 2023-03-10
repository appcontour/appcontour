#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "contour.h"
#include "parser.h"

extern int debug;

static int has_huffman_labelling = 0;

struct sketch *
readsketch (FILE *file)
{
  int tok, res;
  struct sketch *sketch;
  int arcnum = 0;
  int regionnum = 0;

  /* leggi una descrizione in formato descrizione regioni */
  tok = gettoken (file);
  if (tok != TOK_LBRACE)
  {
    fprintf (stderr, "Error: left brace expected\n");
    return (0);
  }
  sketch = newsketch ();
  has_huffman_labelling = 0;
  while ((tok = gettoken (file)) != TOK_RBRACE)
  {
    //ungettoken (tok);
    switch (tok)
    {
      case TOK_ARC:
      if (debug) printf ("leggo la descrizione di un arco\n");
      res = readsketch_arc (arcnum++, sketch, file);
      break;

      case TOK_REGION:
      if (debug) printf ("leggo la descrizione di una regione\n");
      res = readsketch_region (regionnum++, sketch, file);
      break;

      default:
      if (debug) printf ("token invalido: %d\n", tok);
      freesketch (sketch);
      return (0);
      break;
    }
    if (res == 0)
    {
      freesketch (sketch);
      printf ("errore di lettura dal file\n");
      return (0);
    }
    if (debug) printsketch (sketch);
  }
  sketch->huffman_labelling = has_huffman_labelling;
  if (sketch->regions->next == 0) {
    fprintf (stderr, "Warning: empty sketch!\n");
    sketch->isempty = sketch->huffman_labelling = 1;
  }
  if (debug) printsketch (sketch);
  postprocesssketch (sketch);
  return (sketch);
}

int
readsketch_arc (int arcid, struct sketch *sketch, FILE *file)
{
  struct arc *arc;
  int tok, i, cusps_no_d, depthsdim, buf[100];

  arc = newarc (sketch);
  tok = gettoken (file);
  if (tok != ISNUMBER)
  {
    printf ("arc number expected\n");
    return (0);
  }
  arc->tag = gettokennumber ();
  if (arc->next && arc->tag > arc->next->tag)
  {
    sketch->arcs = arc->next;
    arc->next = 0;
    insert_arc_in_list (arc, sketch->arcs);
  }
  if (debug) printf ("leggo la descrizione dell'arco %d\n", arc->tag);

  if ((tok = gettoken (file)) != TOK_COLON)
  {printf ("colon expected\n"); return (0);}
  tok = gettoken (file);
  if (tok != TOK_LPAREN && tok != TOK_LBRACKET)
  {printf ("'(' or '[' expected\n"); return (0);}

  depthsdim = 0;
  cusps_no_d = 0;
  while ((tok = gettoken (file)) == ISNUMBER || tok == KEY_CUSP)
  {
    if (tok == KEY_CUSP) {cusps_no_d++; continue;}
    if (depthsdim > 98) {printf ("too many d values for arc %d\n", arc->tag);
                 return (0);}
    buf[depthsdim++] = gettokennumber ();
    has_huffman_labelling = 1;
  }
  assert (cusps_no_d == 0 || depthsdim == 0);
  if (debug) printf ("letti %d valori di d\n", depthsdim);
  arc->depths = (int *) malloc (depthsdim * sizeof (int));
  arc->depthsdim = depthsdim;
  arc->cusps = depthsdim - 1;
  if (arc->cusps < 0) arc->cusps = 0;
  if (cusps_no_d > 0) arc->cusps = cusps_no_d;
  for (i = 0; i < depthsdim; i++) arc->depths[i] = buf[i];

  if (tok != TOK_RPAREN && tok != TOK_RBRACKET)
  {printf ("')' or ']' expected\n"); return (0);}
  if ((tok = gettoken (file)) != TOK_SEMICOLON)
  {printf ("';' expected\n"); return (0);}
  if (debug) printsketch (sketch);
  return (1);
}

int
readsketch_region (int regionid, struct sketch *sketch, FILE *file)
{
  struct region *region;
  //struct borderlist *bl;
  int tok;
  static int printwarning = 1;

  region = newregion (sketch);
  assert (region->border == 0);
  tok = gettoken (file);
  if (tok != ISNUMBER)
  {
    printf ("region number expected\n");
    return (0);
  }
  region->tag = gettokennumber ();
  if (sketch->regions->next == region) /* region was inserted in second place */
  {
    /* swap first two positions */
    sketch->regions->next = region->next;
    region->next = sketch->regions;
    sketch->regions = region;
  }
  assert (sketch->regions == region);
  if (region->next && region->tag > region->next->tag)
  {
    sketch->regions = region->next;
    region->next = 0;
    insert_region_in_list (region, sketch->regions);
  }
  if (debug) printf ("leggo la descrizione della regione %d\n", region->tag);

  if ((tok = gettoken (file)) == TOK_LPAREN)
  {
    if (printwarning)
      fprintf (stderr, "warning: definition of f is skipped\n");
    printwarning = 0;
    if (gettoken (file) != KEY_F)
    {fprintf (stderr, "'f' expected\n"); return (0);}
    if (gettoken (file) != TOK_EQUAL)
    {fprintf (stderr, "'=' expected\n"); return (0);}
    tok = gettoken (file);
    if (tok == TOK_MINUS) tok = gettoken (file);
    if (tok != ISNUMBER)
    {fprintf (stderr, "number expected\n"); return (0);}
    if (gettoken (file) != TOK_RPAREN)
    {fprintf (stderr, "')' expected\n"); return (0);}
    tok = gettoken (file);
  }
  if (tok != TOK_COLON)
  {printf ("colon expected\n"); return (0);}
  while ((tok = gettoken (file)) == TOK_LPAREN)
  {
    //bl = readsketch_bl (region, sketch, file);
    readsketch_bl (region, sketch, file);
    //bl->next = region->border;
    //region->border = bl;
    if (debug) printsketch (sketch);
  }

  if (tok != TOK_SEMICOLON)
  {printf ("';' expected: %d\n", tok); return (0);}

  return (1);
}

struct borderlist *
readsketch_bl (struct region *r, struct sketch *sketch, FILE *file)
{
  struct border *b, *blast = 0;
  struct borderlist *bl;
  struct arc *arc;
  int tok, atag;

  /* la parentesi aperta e' gia stata letta! */

  if (debug) printf ("componente connessa\n");
  bl = newborderlist_tail (r);
  while ((tok = gettoken (file)) != TOK_RPAREN)
  {
    if (tok != TOK_PLUS && tok != TOK_MINUS)
    {ungettoken (tok); tok = TOK_PLUS; printf ("warning: + assumed\n");}
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    if (tok == TOK_MINUS) b->orientation = -1;

    if (tolower(mygetchar(file)) != 'a')
    {printf ("'a' expected\n"); return (0);}

    if (gettoken (file) != ISNUMBER)
    {printf ("number expected\n"); return (0);}

    atag = gettokennumber ();
    if (debug) printf ("atag = %d\n", atag);
    for (arc = sketch->arcs; arc; arc = arc->next)
    {
      if (arc->tag == atag)
      {
        b->info = arc;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
  }
  if (debug)
  {
    if (bl->sponda) printborder (bl->sponda, bl->region);
      else printf (" () ");
    printf ("\n");
  }
  return (bl);
}

/* local functions */

void
insert_arc_in_list (struct arc *arc, struct arc *list)
{
  struct arc *a;

  assert (arc->tag > list->tag);
  for (a = list; a; a = a->next)
  {
    if (a->next == 0 || a->next->tag > arc->tag)
    {
      arc->next = a->next;
      a->next = arc;
      return;
    }
  }
  assert (0);
  return;
}

void
insert_region_in_list (struct region *region, struct region *list)
{
  struct region *r;

  assert (region->tag > list->tag);
  for (r = list; r; r = r->next)
  {
    if (r->next == 0 || r->next->tag > region->tag)
    {
      region->next = r->next;
      r->next = region;
      return;
    }
  }
  assert (0);
  return;
}
