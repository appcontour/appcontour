#include <assert.h>
#include "contour.h"

/*
 * scrive su stdout una descrizione "morse" di sketch
 */

#define MTYPE_MINF 1
#define MTYPE_PINF 2
#define MTYPE_TOP 3
#define MTYPE_BOT 4
#define MTYPE_CROSS 5
#define MTYPE_TRAN 6

struct morse {
  int type;
  struct border *left;
  struct border *right;
  struct morse *prev;
  struct morse *next;
  };

int advance_morseline (struct morse *mline);
struct morse *handle_mblock (struct morse *block);

extern int debug;
static struct morse minf;

void
printmorse (struct sketch *sketch)
{
  struct morse *block, *prevblock;
  struct borderlist *bl;
  struct border *b;
  struct region *r0;
  int res;

  fprintf (stderr, "NOT YET IMPLEMENTED\n");
  if (! debug) return;

  prevblock = &minf;
  minf.type = MTYPE_MINF;
  minf.prev = minf.next = 0;
  printf ("morse {\n");

  /* inizializzazione linea di morse */

  r0 = sketch->regions;
  bl = r0->border;
  assert (bl->sponda == 0);
  for (bl = bl->next; bl; bl = bl->next)
  {
    block = (struct morse *) malloc (sizeof (struct morse));
    block->type = MTYPE_TOP;
    b = bl->sponda;
    block->left = b;
    block->right = gettransborder (b);
    prevblock->next = block;
    block->prev = prevblock;
    prevblock = block;
  }

  while (minf.next)
  {
    res = advance_morseline (&minf);
  }
  printf ("LAVORI IN CORSO!\n");
  printf ("}\n");
}

/*
 * avanzamento linea di morse
 */

int
advance_morseline (struct morse *mline)
{
  struct morse *block;

  printf ("LINE: ");
  for (block = mline; block; block = block->next)
  {
    printf ("%d ", block->type);
    if (block->next) assert (block->next->prev == block);
  }
  printf ("\n");

  block = mline->next;
  while (1)
  {
    if (block == 0) break;
    block = handle_mblock (block);
  }
  printf (" ;\n");
  return (1);
}

struct morse *
handle_mblock (struct morse *block)
{
  struct morse *next, *inner;
  struct border *b;
  struct borderlist *bl, *blp;
  struct region *r;

  next = block->next;
  switch (block->type)
  {
    case MTYPE_TOP:
      printf (" ^ ");
      next = (struct morse *) malloc (sizeof (struct morse));
      next->type = MTYPE_TRAN;
      block->type = MTYPE_TRAN;
      next->next = block->next;
      if (block->next) block->next->prev = next;
      next->prev = block;
      block->next = next;
      next->left = block->right;
      next->right = block->left;
      b = block->right;
      bl = b->border;
      r = bl->region;
      assert (r->border == bl);
      for (blp = bl->next; blp; blp = blp->next)
      {
        inner = (struct morse *) malloc (sizeof (struct morse));
        inner->type = MTYPE_TOP;
        inner->left = b;
        inner->right = gettransborder (b);
        inner->prev = block;
        inner->next = block->next;
        block->next = inner;
        inner->next->prev = inner;
      }
      return (next->next);

    case MTYPE_BOT:
      printf (" U ");
      next = block->next;
      if (next) next->prev = block->prev;
      block->prev->next = next;
      free (block);
      return (next);

    case MTYPE_TRAN:
      printf (" | ");
      next = block->next;
      if (block->right == next->left)            /* circuito da chiudere */
      {
        block->type = MTYPE_BOT;
        block->next = next->next;
        if (next->next) next->next->prev = block;
        free (next);
        printf (" | ");
        return (block->next);
      }
      return (next);

    case MTYPE_CROSS:
      printf ("CROSS NON IMPLEMENTATO\n");
      exit (1);
  }
  return (next);
}
