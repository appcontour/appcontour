/*
 * trasforma una descrizione di un nodo (stile 'morse')
 * nella descrizione di superfici di tori sottili che
 * simulano il nodo
 *
 * si tratta di una semplice trasformazione dell'input:
 *
 * '^' ==>  '  ^  '
 *          ' / \ '
 *          '/ ^ \'
 *
 * '/' ==>  '  / /'
 *          ' / / '
 *          '/ /  '
 *
 * 'X' ==>  '\ X /'
 *          ' X X '
 *          '/ X \'
 *
 * 'U' ==>  '\ U /'
 *          ' \ / '
 *          '  U  '
 *
 * 'Y' ==>  '\ U /'
 *          ' \ / '
 *          ' | | '
 *
 * 'y' ==>  '\ U /'
 *          ' \ / '
 *          ' | | '
 *
 * 'h' ==>  ' | | '
 *          ' / \ '
 *          '/ ^ \'
 *
 * ',' ==>  '     '
 *          '  ^  '
 *          ' | | '
 *
 * ''' ==>  ' | | '
 *          '  U  '
 *          '     '
 *
 * there are two type of crossings: x and X
 * x has the overpass travelling in the direction NW-SE (backslash)
 * X has the overpass travelling in the direction NE-SW (slash)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "contour.h"
#include "parser.h"

#include "fundamental.h"

struct presentation *knot_update_presentation (struct presentation *p, int k_event, int i);

#define LINESIZE 200
#define MULTLINES 5
#define EATSPACES 4

char charlist[] = " ^XxU/\\()|Yyh,'";
char *subst[] = {
  "             ",
  "             ",
  "             ",
  "             ",
  "             ",

  "      ^l0    ",
  "    /   \\    ",
  "   /     \\   ",
  " /    ^r0  \\ ",
  "/   /   \\   \\",

  "\\     Xd0   /",
  " \\  /  \\u2/  ",
  "  Xd0   Xu0  ",
  " / \\d2 / \\u0 ",
  "/     Xu0d0\\ ",
  
  "\\     Xd2   /",
  " \\  /  \\u0/  ",
  "  Xd0   Xu2  ",
  " / \\d0 / \\u0 ",
  "/     Xu0d0\\ ",
  
  "\\   \\   /   /",
  " \\    U    / ",
  "   \\     /   ",
  "    \\   /    ",
  "      U      ",

  "        /   /",
  "      /   /  ",
  "    /   /    ",
  "  /   /      ",
  "/   /        ",

  "\\   \\        ",
  "  \\   \\      ",
  "    \\   \\    ",
  "      \\   \\  ",
  "        \\   \\",

  "        /   /",
  "       /   / ",
  "      |   |  ",
  "       \\   \\ ",
  "        \\   \\",

  "\\   \\        ",
  " \\   \\       ",
  "  |   |      ",
  " /   /       ",
  "/   /        ",

  "    |  |     ",
  "    |  |     ",
  "    |  |     ",
  "    |  |     ",
  "    |  |     ",

  "\\   \\   /   /",
  " \\    U    / ",
  "   \\     /   ",
  "    |   |    ",
  "    |   |    ",

  "\\   \\   /   /",
  " \\    U    / ",
  "   \\     /   ",
  "    |   |    ",
  "    |   |    ",

  "    |   |    ",
  "    |   |    ",
  "   /     \\   ",
  " /    ^r0  \\ ",
  "/   /   \\   \\",

  "             ",
  "             ",
  "      ^l0    ",
  "    /   \\    ",
  "    |   |    ",

  "    |   |    ",
  "    \\   /    ",
  "      U      ",
  "             ",
  "             ",

  0};

static int eatenspaces = 0;

int knot_getline (char *line, int linesize, FILE *file);
void outpatch (char *p);

int
knot2morse (FILE *file)
{
  int i, j, tok;
  char ch, *linept;
  char line[LINESIZE+2];

  printf ("morse {\n");
  tok = gettoken (file);
  if (tok != TOK_KNOT)
  {
    fprintf (stderr, "Keyword 'knot' expected\n");
    exit (1);
  }
  tok = gettoken (file);
  if (tok != TOK_LBRACE)
  {
    fprintf (stderr, "Left brace '{' expected\n");
    exit (1);
  }

  while (knot_getline (line, LINESIZE, file))
  {
    if (*line == '#' || *line == '\n') continue;
    for (j = 0; j < MULTLINES; j++)
    {
      linept = line;
      eatenspaces = 0;
      while ((ch = *linept++) && ch != ';')
      {
        if (ch == '%') ch = 'X';
        for (i = 0; charlist[i]; i++)
        {
          if (charlist[i] == ch) break;
        }
        if (charlist[i] != ch)
        {
          if (ch == '}') {printf ("}\n"); return (1);}
          fprintf (stderr, "Invalid char '%c'\n", ch);
          exit (1);
        }
        outpatch (subst[MULTLINES*i + j]);
      }
      outpatch ("       ;\n");
    }
  }
  printf ("}\n");
  return (1);
}

int
knot_getline (char *line, int linesize, FILE *file)
{
  int j = 0;
  char ch;

  while ((ch = fgetc (file)) != EOF)
  {
    if (ch == '\n') continue;
    if (ch == '#')
    {
      while ((ch = fgetc (file)) != '\n' && ch != EOF);
      continue;
    }
    assert (j < linesize - 1);
    line[j++] = ch;
    if (ch == ';')
    {
      line[j++] = 0;
      return (j);
    }
  }
  return (0);
}

void
outpatch (char *b)
{
  int i, lb;
  char bloc[80];

  for (i = 0; i < eatenspaces; i++)
  {
    if (*b == ' ') b++;
  }

  lb = strlen(b);
  assert (lb < 80);
  if (lb <= EATSPACES)
  {
    fprintf (stderr, "problems with '%s'\n", b);
    exit (1);
  }
  strcpy (bloc, b);
  for (i = 1; i <= EATSPACES; i++)
  {
    if (bloc[lb - i] == ' ') bloc[lb - i] = 0;
      else break;
  }
  eatenspaces = EATSPACES - lb + strlen(bloc);
  printf ("%s", bloc);
}

int
any2morse (FILE *file)
{
  int tok;
  char ch;
  struct sketch *sketch;

  tok = gettoken (file);
  switch (tok)
  {
    case TOK_KNOT:
    ungettoken (tok);
    return (knot2morse (file));

    case TOK_MORSE:
    /* pipe input file unmodified */
    printf ("morse {");
    if (gettoken (file) != TOK_LBRACE) fprintf (stderr, "Warning: `{' char expected\n");
    while ((ch = fgetc (file)) != EOF) printf ("%c", ch);
    return (1);

    case TOK_SKETCH:
    ungettoken (tok);
    if ((sketch = readcontour (file)) == 0) exit (14);
    printmorse (sketch);
    return (1);

    default:
    fprintf (stderr, "Invalid description type, token: %d\n", tok);
    return (0);
  }
  return (1);
}

/*
 * ==================================================================
 * build fundamental group from knot description
 * ==================================================================
 */

#define KNOT_TRAV 1
#define KNOT_TOP 2
#define KNOT_BOTTOM 3
#define KNOT_CROSSX 4
#define KNOT_CROSSx 5
#define KNOT_FORKTOP 6
#define KNOT_FORKBOTTOM 7
#define KNOT_COMMA 8
#define KNOT_QUOTE 9
#define KNOT_SPACE 100

//   char charlist[] = " ^XxU/\\()|Yyh,'";
// The following MUST correspond to the above as position
static int events[] = {KNOT_SPACE,
  KNOT_TOP,
  KNOT_CROSSX,
  KNOT_CROSSx,
  KNOT_BOTTOM,
  KNOT_TRAV,
  KNOT_TRAV,
  KNOT_TRAV,
  KNOT_TRAV,
  KNOT_TRAV,
  KNOT_FORKBOTTOM,
  KNOT_FORKBOTTOM,
  KNOT_FORKTOP,
  KNOT_COMMA,
  KNOT_QUOTE};


struct presentation *
knot2fg (FILE *file)
{
  char ch, *linept;
  char line[LINESIZE+2];
  int i, k, k_event, tok;
  struct presentation *p;

  tok = gettoken (file);
  if (tok != TOK_KNOT)
  {
    fprintf (stderr, "Keyword 'knot' expected\n");
    exit (1);
  }
  tok = gettoken (file);
  if (tok != TOK_LBRACE)
  {
    fprintf (stderr, "Left brace '{' expected\n");
    exit (1);
  }

  i = 0;

  p = (struct presentation *) malloc (sizeof (struct presentation));
  p->gennum = 0;
  p->rules = p->elements = 0;
  p->characteristic = 1000;  /* the spanning tree has undefined characteristic */

  while (knot_getline (line, LINESIZE, file))
  {
    if (*line == '#' || *line == '\n') continue;
    linept = line;
    while ((ch = *linept++) && ch != ';')
    {
      if (isspace (ch)) continue;
      if (ch == '%') ch = 'X';
      for (k = 0; charlist[k]; k++) if (charlist[k] == ch) {k_event = events[k]; break;}

      if (charlist[k] != ch)
      {
        if (ch == '}') return (p);
        fprintf (stderr, "Invalid char '%c'\n", ch);
        exit (1);
      }

      if (k_event == KNOT_TRAV)
      {
        i++;
        continue;
      }

      p = knot_update_presentation (p, k_event, i);

      switch (k_event)
      {
        case KNOT_FORKBOTTOM:
        case KNOT_COMMA:
        i++;
        break;

        case KNOT_TOP:
        case KNOT_CROSSX:
        case KNOT_CROSSx:
        case KNOT_FORKTOP:
        i += 2;
        break;
      }

    }
    if (ch == ';')
    {
      i = 0;
    } else {
      printf ("MI ASPETTAVO un ;\n");
    }
  }
  //printf ("END of BUILD\n");

  return (p);
}


struct presentation *
knot_update_presentation (struct presentation *p, int k_event, int pos_event)
{
  int num_head = 0;
  int num_tail = 0;
  int k, kk;
  struct presentationrule *word;
  struct presentationrule *tail;
  struct presentationrule *elem, *elemi, *elemii;

  if (pos_event == 0)
  {
    tail = p->elements;
    p->elements = 0;
    num_head = 0;
  } else {
    for (word = p->elements; word; word = word->next)
    {
      num_head++;
      if (num_head == pos_event)
      {
        tail = word->next;
        word->next = 0;
        break;
      }
    }
    assert (word);
  }
  for (word = tail; word; word = word->next) num_tail++;

  //printf ("In update presentation, event: %d, pos: %d, words before: %d, after: %d\n", k_event, pos_event, num_head, num_tail);

  switch (k_event)
  {
    case KNOT_TOP:
    p->gennum++;
    //printf ("Qui aggiungo due elementi: generatore %d e il suo inverso\n", p->gennum);
    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + 1*sizeof (int));
    elem->length = 1;
    elem->next = tail;
    elem->var[0] = - p->gennum;
    tail = elem;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + 1*sizeof (int));
    elem->length = 1;
    elem->next = tail;
    elem->var[0] = p->gennum;
    tail = elem;

    break;

    case KNOT_CROSSX:
    assert (tail);
    elemi = tail;
    tail = tail->next;
    assert (tail);
    elemii = tail;
    tail = tail->next;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + (2*elemii->length + elemi->length)*sizeof (int));
    elem->length = 2*elemii->length + elemi->length;
    kk = 0;
    for (k = elemii->length - 1; k >= 0; k--)
      elem->var[kk++] = - elemii->var[k];
    for (k = 0; k < elemi->length; k++)
      elem->var[kk++] = elemi->var[k];
    for (k = 0; k < elemii->length; k++)
      elem->var[kk++] = elemii->var[k];

    elem->next = tail;
    tail = elem;
    elemii->next = tail;
    tail = elemii;
    free (elemi);

    break;

    case KNOT_CROSSx:
    assert (tail);
    elemi = tail;
    tail = tail->next;
    assert (tail);
    elemii = tail;
    tail = tail->next;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + (2*elemi->length + elemii->length)*sizeof (int));
    elem->length = 2*elemi->length + elemii->length;
    kk = 0;
    for (k = 0; k < elemi->length; k++)
      elem->var[kk++] = elemi->var[k];
    for (k = 0; k < elemii->length; k++)
      elem->var[kk++] = elemii->var[k];
    for (k = elemi->length - 1; k >= 0; k--)
      elem->var[kk++] = - elemi->var[k];

    elemi->next = tail;
    tail = elemi;
    elem->next = tail;
    tail = elem;
    free (elemii);

    break;

    case KNOT_FORKTOP:
    p->gennum++;
    assert (tail);
    elemi = tail;
    tail = tail->next;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + (1 + elemi->length)*sizeof (int));
    elem->length = 1 + elemi->length;

    elem->var[0] = - p->gennum;
    for (k = 0; k < elemi->length; k++)
    {
      elem->var[k+1] = elemi->var[k];
    }
    elem->next = tail;
    tail = elem;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + (1)*sizeof (int));
    elem->length = 1;
    elem->var[0] = p->gennum;

    elem->next = tail;
    tail = elem;

    free (elemi);
    break;

    case KNOT_FORKBOTTOM:

    assert (tail);
    elemi = tail;
    tail = tail->next;
    assert (tail);
    elemii = tail;
    tail = tail->next;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + (elemi->length + elemii->length)*sizeof (int));
    elem->length = elemi->length + elemii->length;

    kk = 0;
    for (k = 0; k < elemi->length; k++)
      elem->var[kk++] = elemi->var[k];
    for (k = 0; k < elemii->length; k++)
      elem->var[kk++] = elemii->var[k];

    elem->next = tail;
    tail = elem;
    free (elemi);
    free (elemii);
    break;

    case KNOT_BOTTOM:

    assert (tail);
    elemi = tail;
    tail = tail->next;
    assert (tail);
    elemii = tail;
    tail = tail->next;

    elem = (struct presentationrule *) malloc (sizeof (struct presentationrule) + (elemi->length + elemii->length)*sizeof (int));
    elem->length = elemi->length + elemii->length;

    kk = 0;
    for (k = 0; k < elemi->length; k++)
      elem->var[kk++] = elemi->var[k];
    for (k = 0; k < elemii->length; k++)
      elem->var[kk++] = elemii->var[k];
    free (elemi);
    free (elemii);
    elem->next = p->rules;
    p->rules = elem;
    break;

    case KNOT_COMMA:
    case KNOT_QUOTE:
printf ("evento %d non gestito\n", k_event);
assert (false);
    break;

    default:
    printf ("Unhandled event %d\n", k_event);
    assert (false);
    break;
  }

  /* merge head with tail */

  if (p->elements == 0)
  {
    p->elements = tail;
  } else {
    for (word = p->elements; word; word = word->next)
    if (word->next == 0)
    {
      word->next = tail;
      break;
    }
  }

  //print_presentation (p);

  return (p);
}
