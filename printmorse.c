#include <stdio.h>
#include <assert.h>
#include "contour.h"

extern int debug;

static int crbefore = 0;
static int crafter = 0;

struct morse {
  struct border *left;
  struct border *right;
  struct morse *prev;
  struct morse *next;
  struct morse *inf;
};

/* local prototypes */

void morse_islands (struct borderlist *bl);
void nested_morse_islands (struct morse *left, struct borderlist *bl);
void enter_arc_event (struct border *b, struct morse *before);
struct morse *newmblock (struct morse *before);
void do_morse (struct morse *minf);
void advance_morse_line (struct morse *minf);
void enter_region (struct morse *left, int type);
void nested_morse (struct morse *mb1, struct morse *mb2, int arcbelow);
int isonleft (struct morse *mb1, struct morse *mb2);
void morse_close (struct morse *left, struct morse *right);
struct morse *segui_bordo (struct border *bin, struct morse *mbin, int type);
int countleft (struct morse *m);
int countright (struct morse *m);
void outrowp (int type, struct morse *left, struct border *bll, struct border *brl);
void outrow (int type, struct border *bll, struct border *brl);
void outinfodepths (struct arc *arc);

#define TYPE_TOP 1
#define TYPE_CROSS 2
#define TYPE_BOT 3

void
printmorse (struct sketch *sketch)
{
  struct arc *arc;

  if (debug) printf ("Entering in printmorse\n");
  printf ("morse {\n");
  morse_islands (sketch->regions->border->next);
  printf ("}\n");
  /* now reset all arc tags */
  for (arc = sketch->arcs; arc; arc = arc->next)
  {
    if (arc->tag < 0) arc->tag *= -1;
      else fprintf (stderr, "Error: forgotten arc %d\n", arc->tag);
  }
  if (debug) printf ("  printmorse finished\n");
}

void
morse_islands (struct borderlist *bl)
{
  struct morse minf;
  minf.inf = 0;

  for (; bl; bl = bl->next)
  {
    minf.next = minf.prev = &minf;
    enter_arc_event (bl->sponda, &minf);  
    /* volendo posso decidere di inserire tutto il bordo... */
    do_morse (&minf);
  }
}

void
enter_arc_event (struct border *b, struct morse *before)
{
  struct morse *left, *right;
  struct border *bin;

  if (debug) printf ("New arc event (%d), this is a '^'\n",
    b->info->tag);
  bin = gettransborder (b);
  /* devo creare due nuovi blocchi di morse */
  left = newmblock (before);
  left->left = b;
  left->right = bin;
  right = newmblock (left);
  right->left = bin;
  right->right = b;
  outrowp (TYPE_TOP, left, b, 0);

  enter_region (left, TYPE_TOP);
  return;
}

struct morse *
newmblock (struct morse *before)
{
  struct morse *block;

  block = (struct morse *) malloc (sizeof (struct morse));
  block->prev = before;
  block->next = before->next;
  before->next = block;
  block->next->prev = block;
  if (before->inf) block->inf = before->inf;
    else block->inf = before;
  return (block);
}

void
do_morse (struct morse *minf)
{
  /* move the morse line untill empty */
  while (minf->next->inf)
  {
    advance_morse_line (minf);
  }
}

void
advance_morse_line (struct morse *minf)
{
  struct morse *left, *right;
  struct border *bgen, *bll, *brl;

  assert (minf->next->inf);

  /* lavoro sul primo intervallo */
  left = minf->next;
  right = left->next;
  assert (right->inf);

  /* ora sfrutto il fatto che non ci sono archi-sotto */
  if (left->right->next == right->left)
  {
    if (debug) printf ("advance: crossing\n");
    /* attraverso un crossing */
    bgen = right->right;
    bgen = bgen->next;    /* e' il nuovo brr */
    right->right = bgen;
    brl = gettransborder (bgen);  /* e' il nuovo brl */
    right->left = brl;
    bgen = brl->next;             /* e' il nuovo blr */
    left->right = bgen;
    bll = gettransborder (bgen);  /* e' il nuovo bll */
    assert (bll->next == left->left);
    left->left = bll;
    outrowp (TYPE_CROSS, left, bll, brl);
    if (debug) printf ("Attraverso un nodo, questo e' un 'X'\n");
    if (debug) printf ("  nuovi archi: %d %d\n", left->left->info->tag,
                right->right->info->tag);
    enter_region (left, TYPE_CROSS);
  } else {
    if (debug) printf ("advance: attraverso arco\n");
    /* inserisco un attraversamento-arco */
    enter_arc_event (left->right->next, left);
  }
}

void
enter_region (struct morse *left, int type)
{
  struct morse *right, *destblock;
  struct border *binl, *binr;
  struct arc *arc;
  int arcbelow, arcbelow1, arcbelow2;

  if (debug) printf ("in enter_region: ");
  assert (left->inf);
  right = left->next;
  assert (right->inf);

  /* devo controllare se ho un "arco-sotto" */
  binl = left->right;
  binr = right->left;
  if (debug) printf ("%d\n", binl->border->region->tag);

  arcbelow = 0;
  if (binl == binr) /* stesso arco */
  {
    arc = binl->info;
    if (arc->endpoints == 0)   /* e' un S1 */
    {
      assert (type == TYPE_TOP);
      arcbelow = 1;
    } else {
      if (type == TYPE_CROSS) arcbelow = 1;
    }
  }

  if (arcbelow)
  {
    /* in base alle regole di consistenza questo
     * e' possibile solo se la regione prima e dopo
     * e' la regione esterna... ovvero posso chiudere
     * la descrizione di morse. Questa e' per forza
     * una regione nuova
     */
    if (debug) printf ("  e.r. arcbelow\n");
    //crbefore++;
    //crafter++;
    assert (binl->border == binl->border->region->border);
    nested_morse_islands (left, binl->border->next);
    assert (left->prev == left->inf);
    assert (right->next == right->inf);
    //crbefore--;
    //crafter--;
    morse_close (left, right);
    return;
  }

  /* altrimenti devo controllare se e' una regione nuova */

  destblock = segui_bordo (binl, left, type);
  if (destblock == right)
  {             /* e' una regione nuova */
    if (debug) printf ("  e.r. regione nuova\n");
    assert (binl->border == binl->border->region->border);
    nested_morse_islands (left, binl->border->next);
    return;
  }

  if (debug) printf ("  e.r. regione vecchia\n");
  assert (destblock != left);
  /* non e' una regione nuova... */
  //if (destblock == minf->next)
  //{
    /* e' la regione esterna */
    /* forse non va trattata in modo speciale */
  //}

  arcbelow1 = arcbelow2 = 0;
  if (destblock->left == left->right);
  if (type == TYPE_TOP)
  {
    assert (destblock->left != left->right);
    assert (destblock->prev->right != left->next->left);
  } else {
    if (destblock->left == left->right) arcbelow1 = 1;
    if (destblock->prev->right == left->next->left) arcbelow2 = 1;
  }
  if (isonleft (destblock, left))
  {
    nested_morse (destblock, left, arcbelow1);
    if (arcbelow2) morse_close (destblock->prev, left->next);
  } else {
    nested_morse (left->next, destblock->prev, arcbelow2);
    if (arcbelow1) morse_close (left, destblock);
  }
  return;
}

void
nested_morse_islands (struct morse *left, struct borderlist *bl)
{
  int cl, cr;

  cl = countleft (left) + 1;
  cr = countright (left->next) + 1;
  if (debug) printf ("in nested_morse_islands, cl = %d, cr = %d\n", cl, cr);
  crbefore += cl;
  crafter += cr;

  morse_islands (bl);

  crbefore -= cl;
  crafter -= cr;
  if (debug) printf ("exiting from nested_morse_islands\n");
}

void
nested_morse (struct morse *mb1, struct morse *mb2, int arcbelow)
{
  struct morse minf, *block;
  int cl, cr;

  if (debug) printf ("entering nested morse\n");

  cl = countleft (mb1);
  cr = countright (mb2);
  if (arcbelow)
  {
    if (mb1->next != mb2) nested_morse (mb1->next, mb2->next, 0);
    crbefore += cl;
    crafter += cr;
    morse_close (mb1, mb2);
    crbefore -= cl;
    crafter -= cr;
    return;
  }
  crbefore += cl;
  crafter += cr;
  /* devo creare una nuova linea di morse */
  minf.inf = 0;
  minf.next = mb1;
  minf.prev = mb2;

  /* ritaglio dalla vecchia linea di morse */
  mb1->prev->next = mb2->next;
  mb2->next->prev = mb1->prev;
  mb1->prev = &minf;
  mb2->next = &minf;
  for (block = mb1; block->inf; block = block->next)
    block->inf = &minf;

  do_morse (&minf); 
  crbefore -= cl;
  crafter -= cr;
}

int
isonleft (struct morse *mb1, struct morse *mb2)
{
  struct morse *mb;

  for (mb = mb1; mb->inf; mb = mb->next)
  {
    if (mb == mb2) return (1);
  }
  return (0);
}

/*
 * TYPE_TOP: non dobbiamo essere in una situazione di tipo "arco-sotto"
 * TYPE_CROSS: possiamo essere in una situazione di tipo "arco-sotto"
 */

struct morse *
segui_bordo (struct border *bin, struct morse *mbin, int type)
{
  struct morse *minf, *mblock;
  struct border *bp;

  assert (mbin->right == bin);
  minf = mbin->inf;
  bp = bin;
  if (type == TYPE_CROSS)
  {
    for (mblock = minf->next; mblock->inf; mblock = mblock->next)
    {
      if (mblock->left == bp) return (mblock);
    }
  }
  while (1)
  {
    bp = bp->next;
    //printf ("BP: %d\n", bp->info->tag);
    /* cerco bp sulla linea di morse */
    for (mblock = minf->next; mblock->inf; mblock = mblock->next)
    {
    //printf ("mblock->left = %d\n", mblock->left->info->tag);
      if (mblock->left == bp) return (mblock);
    }
    assert (bp != bin);
  }
  return (0);
}

void
morse_close (struct morse *left, struct morse *right)
{
  if (debug) printf ("closing a morse description (with an U).\n");
  outrow (TYPE_BOT, 0, 0);
  left->prev->next = right->next;
  right->next->prev = left->prev;
  free (left);
  free (right);
}

int
countleft (struct morse *m)
{
  struct morse *block;

  int count = 0;
  for (block = m->prev; block->inf; block = block->prev) count++;
  return (count);
}

int
countright (struct morse *m)
{
  struct morse *block;

  int count = 0;
  for (block = m->next; block->inf; block = block->next) count++;
  return (count);
}

void
outrowp (int type, struct morse *left, struct border *bll, struct border *brl)
{
  int cl, cr;

  cl = countleft (left);
  cr = countright (left->next);

  crbefore += cl;
  crafter += cr;
  outrow (type, bll, brl);
  crafter -= cr;
  crbefore -= cl;
}

void
outrow (int type, struct border *bll, struct border *brl)
{
  int i;

  if (debug) printf ("ROW: ");
  for (i = 0; i < crbefore; i++) printf (" | ");
  if (bll && bll->info->tag < 0) bll = 0;
  switch (type)
  {
    case TYPE_TOP:
      assert (brl == 0);
      printf (" ^");
      if (bll)
      {
        printf ("%c,", (bll->orientation > 0)?'r':'l');
        outinfodepths (bll->info);
        printf (" ");
      }
      break;

    case TYPE_BOT:
      assert (bll == 0 && brl == 0);
      printf (" U ");
      break;

    case TYPE_CROSS:
      printf (" X ");
      if (bll)
      {
        printf ("%c,", (bll->orientation > 0)?'u':'d');
        outinfodepths (bll->info);
      }
      if (brl && brl->info->tag < 0) brl = 0;
      if (brl)
      {
        if (bll) printf (" "); else printf ("[]");
        printf ("%c,", (brl->orientation > 0)?'u':'d');
        outinfodepths (brl->info);
      }
      if (brl || bll) printf (" ");
      break;

    default:
      printf (" ? ");
      break;
  }
  for (i = 0; i < crafter; i++) printf (" | ");
  printf (";\n");
}

void
outinfodepths (struct arc *arc)
{
  int d, i;

  assert (arc->tag >= 0);  /* we already wrote the information */
  arc->tag *= -1;
  d = arc->depths[0];
  printf ("%d", d);
  for (i = 1; i <= arc->cusps; i++)
  {
    if (abs (arc->depths[i] - d) != 1) fprintf (stderr, 
           "erroneous d values across cusp\n");
    if (arc->depths[i] > d) printf ("+");
    if (arc->depths[i] < d) printf ("-");
    d = arc->depths[i];
  }
}
