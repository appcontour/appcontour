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

/* prototipi */

void morse_islands (struct borderlist *bl);
void enter_arc_event (struct border *b, struct morse *before);
struct morse *newmblock (struct morse *before);
void do_morse (struct morse *minf);
void advance_morse_line (struct morse *minf);
void enter_region (struct morse *left, int type);
void nested_morse (struct morse *mb1, struct morse *mb2, int arcosotto);
int isonleft (struct morse *mb1, struct morse *mb2);
void morse_close (struct morse *left, struct morse *right);
struct morse *segui_bordo (struct border *bin, struct morse *mbin, int type);
int countleft (struct morse *m);
int countright (struct morse *m);
void outrowp (int type, struct morse *left);
void outrow (int type);

#define TYPE_TOP 1
#define TYPE_CROSS 2
#define TYPE_BOT 3

void
printmorse (struct sketch *sketch)
{
  if (debug) printf ("Entering in printmorse\n");
  printf ("morse {\n");
  morse_islands (sketch->regions->border->next);
  printf ("}\n");
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
    /* posso decidere di inserire tutto il bordo... */
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
  outrowp (TYPE_TOP, left);

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
  struct border *bgen;

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
    bgen = gettransborder (bgen);  /* e' il nuovo brl */
    right->left = bgen;
    bgen = bgen->next;             /* e' il nuovo blr */
    left->right = bgen;
    bgen = gettransborder (bgen);  /* e' il nuovo bll */
    assert (bgen->next == left->left);
    left->left = bgen;
    outrowp (TYPE_CROSS, left);
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
  int arcosotto, arcosotto1, arcosotto2;

  //printf ("in enter_region\n");
  assert (left->inf);
  right = left->next;
  assert (right->inf);

  /* devo controllare se ho un "arco-sotto" */
  binl = left->right;
  binr = right->left;

  arcosotto = 0;
  if (binl == binr) /* stesso arco */
  {
    arc = binl->info;
    if (arc->endpoints == 0)   /* e' un S1 */
    {
      assert (type == TYPE_TOP);
      arcosotto = 1;
    } else {
      if (type == TYPE_CROSS) arcosotto = 1;
    }
  }

  if (arcosotto)
  {
    /* in base alle regole di consistenza questo
     * e' possibile solo se la regione prima e dopo
     * e' la regione esterna... ovvero posso chiudere
     * la descrizione di morse. Questa e' per forza
     * una regione nuova
     */
    //printf ("  e.r. arcosotto\n");
    crbefore++;
    crafter++;
    assert (binl->border == binl->border->region->border);
    morse_islands (binl->border->next);
    assert (left->prev == left->inf);
    assert (right->next == right->inf);
    crbefore--;
    crafter--;
    morse_close (left, right);
    return;
  }

  /* altrimenti devo controllare se e' una regione nuova */

  destblock = segui_bordo (binl, left, type);
  if (destblock == right)
  {             /* e' una regione nuova */
    //printf ("  e.r. regione nuova\n");
    assert (binl->border == binl->border->region->border);
    morse_islands (binl->border->next);
    return;
  }

  //printf ("  e.r. regione vecchia\n");
  assert (destblock != left);
  /* non e' una regione nuova... */
  //if (destblock == minf->next)
  //{
    /* e' la regione esterna */
    /* forse non va trattata in modo speciale */
  //}

  arcosotto1 = 0;
  arcosotto2 = 0;
  if (destblock->left == left->right);
  if (type == TYPE_TOP)
  {
    assert (destblock->left != left->right);
    assert (destblock->prev->right != left->next->left);
  } else {
    if (destblock->left == left->right) arcosotto1 = 1;
    if (destblock->prev->right == left->next->left) arcosotto2 = 1;
  }
  if (isonleft (destblock, left))
  {
    nested_morse (destblock, left, arcosotto1);
    if (arcosotto2) morse_close (destblock->prev, left->next);
  } else {
    nested_morse (left->next, destblock->prev, arcosotto2);
    if (arcosotto1) morse_close (left, destblock);
  }
  return;
}

void
nested_morse (struct morse *mb1, struct morse *mb2, int arcosotto)
{
  struct morse minf, *block;
  int cl, cr;

  if (debug) printf ("entering nested morse\n");

  if (arcosotto)
  {
    if (mb1->next != mb2) nested_morse (mb1->next, mb2->next, 0);
    crbefore++;
    crafter++;
    morse_close (mb1, mb2);
    crbefore--;
    crafter--;
    return;
  }
  cl = countleft (mb1);
  cr = countright (mb2);
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
  outrow (TYPE_BOT);
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
outrowp (int type, struct morse *left)
{
  int cl, cr;

  cl = countleft (left);
  cr = countright (left->next);

  crbefore += cl;
  crafter += cr;
  outrow (type);
  crafter -= cr;
  crbefore -= cl;
}

void
outrow (int type)
{
  int i;

  if (debug) printf ("ROW: ");
  for (i = 0; i < crbefore; i++) printf (" | ");
  switch (type)
  {
    case TYPE_TOP:
      printf (" ^ ");
      break;

    case TYPE_BOT:
      printf (" U ");
      break;

    case TYPE_CROSS:
      printf (" X ");
      break;

    default:
      printf (" ? ");
      break;
  }
  for (i = 0; i < crafter; i++) printf (" | ");
  printf (";\n");
}
