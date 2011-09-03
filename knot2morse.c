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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "contour.h"
#include "parser.h"

#define LINESIZE 200
#define MULTLINES 5
#define EATSPACES 4

char charlist[] = " ^XxU/\\()|Yh,'";
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

  while (fgets (line, LINESIZE, file))
  {
    if (*line == '#' || *line == '\n') continue;
    for (j = 0; j < MULTLINES; j++)
    {
      linept = line;
      eatenspaces = 0;
      while ((ch = *linept++) && ch != ';')
      {
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
