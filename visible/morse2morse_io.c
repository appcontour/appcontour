#include "morse2morse.h"

#define TOKENSTACKSIZE 100
#define MAXWORDSIZE 100

#define TOK_EOF 1
#define TOK_ERROR 2
#define TOK_COMMENT 3
#define TOK_SPACE 4
#define TOK_MORSE 5
#define TOK_LBRACE 6
#define TOK_RBRACE 7
#define TOK_SEMICOLON 8
#define TOK_MAX 9
#define TOK_MIN 10
#define TOK_VERT 11
#define TOK_MINUS 12
#define TOK_PLUS 13
#define TOK_LEFT 14
#define TOK_RIGHT 15
#define TOK_UP 16
#define TOK_DOWN 17
#define TOK_X 18
#define TOK_CUSPUP 19
#define TOK_CUSPDOWN 20
#define TOK_SLASH 22
#define TOK_BACKSLASH 23

int tok_stack[TOKENSTACKSIZE];
static int tok_stack_pointer = 0;
static int nowords = 0;

/* a couple of prototypes */
int token2number (int token);
int number2token (int val);
int isanumber (int token);
int getnumber (FILE *file);
void write_cuspsinfo (int dvalues, int *d);

struct mdesc *
read_contour (FILE *file)
{
  int token, tag;
  struct line *line, *prevline, *firstline;
  struct mdesc *contour;

  tok_stack_pointer = 0;
  nowords = 0;
  token = get_token (file);
  assert (token == TOK_MORSE);
  token = get_token (file);
  assert (token == TOK_LBRACE);

  prevline = 0;
  tag = 0;
  while (1)
  {
    nowords = 1;
    token = get_token (file);
    if (token == TOK_RBRACE) break;
    unget_token (token);
//printf ("leggo la riga %d\n", tag);
    line = get_line (file);
    assert (line);
    line->next = 0;
    line->prev = prevline;
    line->tag = tag++;
    if (prevline) prevline->next = line;
      else firstline = line;
    prevline = line;
  }
  contour = (struct mdesc *) malloc (sizeof (struct mdesc));
  contour->lines = firstline;
  contour->lastline = line;
  return (contour);
}

/* ------------------------------------------ */

struct line *
get_line (FILE *file)
{
  struct line *line;
  struct event *event, *prevevent, *firstevent;
  int token, tag;

  line = (struct line *) malloc (sizeof (struct line));
  prevevent = firstevent = 0;
  tag = 0;
  while ( (token = get_token (file)) != TOK_SEMICOLON )
  {
    assert (token != TOK_EOF);
    assert (token != TOK_ERROR);
    event = (struct event *) malloc (sizeof (struct event));
//printf ("event token: %d\n", token);
    if (prevevent) prevevent->next = event;
      else firstevent = event;
    event->next = 0;
    event->tag = tag++;
    event->line = line;
    event->orientation = event->orientation2 = 0;
    event->huffman = event->huffman2 = -1;
    event->cuspsign = 0;
    event->cusps = event->cusps2 = 0;  /* will hold cusps information */
    prevevent = event;
    switch (token)
    {
      case TOK_VERT:
      event->type = EVENT_VERT;
      break;
      case TOK_MAX:
      event->type = EVENT_MAX;
      break;
      case TOK_MIN:
      event->type = EVENT_MIN;
      break;
      case TOK_CUSPUP:
      event->type = EVENT_CUSP;
      event->orientation = ORIENT_UP;
      break;
      case TOK_CUSPDOWN:
      event->type = EVENT_CUSP;
      event->orientation = ORIENT_DOWN;
      break;
      case TOK_X:
      event->type = EVENT_CROSS;
      break;
      default:
      fprintf (stderr, "Invalid token for event: %d\n", token);
    }
    token = get_token (file);
    switch (token)
    {
      case TOK_LEFT:
      assert (event->orientation == 0);
      event->orientation = ORIENT_LEFT;
      break;
      case TOK_RIGHT:
      assert (event->orientation == 0);
      event->orientation = ORIENT_RIGHT;
      break;
      case TOK_UP:
      assert (event->orientation == 0);
      event->orientation = ORIENT_UP;
      break;
      case TOK_DOWN:
      assert (event->orientation == 0);
      event->orientation = ORIENT_DOWN;
      break;
      default:
      unget_token (token);
    }
    token = get_token (file);
    if (isanumber(token))
      event->huffman = token2number (token);
     else unget_token (token);
    if (event->type == EVENT_CUSP)
    {
      token = get_token (file);
      if (token == TOK_PLUS) event->cuspsign = 1;
      if (token == TOK_MINUS) event->cuspsign = -1;
      assert (event->cuspsign);
    }
    if (event->type == EVENT_CROSS)
    {
      token = get_token (file);
      if (token == TOK_UP)
      {
        event->orientation2 = ORIENT_UP;
        token = get_token (file);
        assert (token != TOK_DOWN);
      }
      if (token == TOK_DOWN)
      {
        event->orientation2 = ORIENT_DOWN;
        token = get_token (file);
      }
      if (isanumber (token))
      {
        event->huffman2 = token2number (token);
      } else unget_token (token);
    }

//printf ("orientation token: %d\n", token);
  }
  line->events = firstevent;
  return (line);
}

/* ------------------------------------------ */

int
get_tokenl (FILE *file)
{
  char ch, word[MAXWORDSIZE];
  int val;

  ch = fgetc (file);
  if (ch == EOF) return (TOK_EOF);
  if (ch == '#')
  {
    while ((ch = fgetc (file) != '\n' && ch != EOF));
    return (TOK_COMMENT);
  }
  if (isspace (ch))
  {
    while (isspace( (ch = fgetc(file)) ));
    ungetc (ch, file);
    return (TOK_SPACE);
  }
  if (isalpha (ch) && nowords == 0)
  {
    ungetc (ch, file);
    getword (file, word);
    if (strcmp (word, "morse") == 0) return (TOK_MORSE);
    fprintf (stderr, "unknown word: %s\n", word);
    return (TOK_ERROR);
  }
  if (isdigit (ch))
  {
    ungetc (ch, file);
    val = getnumber (file);
    assert (isanumber (number2token(val)));
    return (number2token(val));
  }
  switch (ch)
  {
    case '{':
      return (TOK_LBRACE);
      break;
    case '}':
      return (TOK_RBRACE);
      break;
    case ';':
      return (TOK_SEMICOLON);
      break;
    case '^':
      return (TOK_MAX);
      break;
    case 'U':
      return (TOK_MIN);
      break;
    case '|':
      return (TOK_VERT);
      break;
    case '/':
      return (TOK_SLASH);
      break;
    case '\\':
      return (TOK_BACKSLASH);
      break;
    case 'X':
      return (TOK_X);
      break;
    case '>':
      return (TOK_CUSPUP);
      break;
    case '<':
      return (TOK_CUSPDOWN);
      break;
    case '-':
      return (TOK_MINUS);
      break;
    case '+':
      return (TOK_PLUS);
      break;
    case 'l':
      return (TOK_LEFT);
      break;
    case 'r':
      return (TOK_RIGHT);
      break;
    case 'u':
      return (TOK_UP);
      break;
    case 'd':
      return (TOK_DOWN);
      break;
  }
  fprintf (stderr, "unknown token %c\n", ch);
  return (TOK_ERROR);
}

int
get_token (FILE *file)
{
  int token;

  if (tok_stack_pointer > 0)
  {
//printf ("token: %d (from stack)\n", tok_stack[tok_stack_pointer-1]);
    return (tok_stack[--tok_stack_pointer]);
  }

  while (1)
  {
    token = get_tokenl (file);
//printf ("tokenl: %d\n", token);
    if (token != TOK_COMMENT && token != TOK_SPACE) return (token);
  }
}

/* ------------------------------------------ */

char *
getword (FILE *file, char *word)
{
  char *chpt, ch;

  chpt = word;
  while ( isalnum(ch = fgetc (file)) ) *chpt++ = ch;
  *chpt = 0;
  ungetc (ch, file);
  return (word);
}

/* ------------------------------------------ */

int
getnumber (FILE *file)
{
  char ch;
  int val = 0;

  while ( isdigit(ch = fgetc (file)) ) val = 10*val + (ch - '0');
  ungetc (ch, file);
  return (val);
}

/* ------------------------------------------ */

void
unget_token (int token)
{
  assert (tok_stack_pointer < TOKENSTACKSIZE);
  tok_stack[tok_stack_pointer++] = token;
}

#define TOKNUMSHIFT 1000
#define TOKNUMTHRESHOLD (-100)

int
token2number (int token)
{
  return (-(TOKNUMSHIFT + token));
}

int
number2token (int val)
{
  return (-(TOKNUMSHIFT + val));
}

int
isanumber (int token)
{
  return (token < TOKNUMTHRESHOLD);
}

/* ------------------------------------------ */

void
write_contour (struct mdesc *contour)
{
  struct line *l;
  struct event *ev;
  char udori, udori2, lrori;
  struct patch *patch;

  printf ("#\n# generated by morse2morse\n#\n");
  printf ("morse {\n");
  for (l = contour->lines; l; l = l->next)
  {
    for (ev = l->events; ev; ev = ev->next)
    {
      udori = udori2 = lrori = '?';
      if (ev->orientation == ORIENT_UP) udori='u';
      if (ev->orientation == ORIENT_DOWN) udori='d';
      if (ev->orientation2 == ORIENT_UP) udori2='u';
      if (ev->orientation2 == ORIENT_DOWN) udori2='d';
      if (ev->orientation == ORIENT_LEFT) lrori='l';
      if (ev->orientation == ORIENT_RIGHT) lrori='r';
      switch (ev->type)
      {
        case EVENT_VERT:
        case EVENT_CUSP:
          assert (udori != '?');
          assert (ev->cusps == 0);
          //printf ("|%c ", udori);
          printf ("| ");
          break;

        case EVENT_CROSS:
          assert (udori != '?');
          assert (udori2 != '?');
          printf ("X ");
          if ((patch = ev->cusps) != 0)
          {
            printf ("%c,", udori);
            write_cuspsinfo (patch->dvalues, patch->d);
          }
          if ((patch = ev->cusps2) != 0)
          {
            if (ev->cusps) printf (" "); else printf ("[]");
            printf ("%c,", udori2);
            write_cuspsinfo (patch->dvalues, patch->d);
          }
          printf (" ");
          break;

        case EVENT_MAX:
          assert (lrori != '?');
          printf ("^");
          if ((patch = ev->cusps) != 0)
          {
            printf ("%c,", lrori);
            write_cuspsinfo (patch->dvalues, patch->d);
          }
          printf (" ");
          break;

        case EVENT_MIN:
          assert (lrori != '?');
          assert (ev->cusps == 0);
          //printf ("U %c ", lrori);
          printf ("U ");
          break;

        default:
          fprintf (stderr, "Unknown event of type %d\n", ev->type);
          exit (4);
      }
    }
    printf (";\n");
  }
  printf ("}\n");
}

void
write_cuspsinfo (int dvalues, int *d)
{
  int i, diff;

  printf ("%d", d[0]);
  for (i = 1; i < dvalues; i++)
  {
    diff = d[i] - d[i-1];
    assert (diff == 1 || diff == -1);
    printf ("%c", (diff == 1)?'+':'-');
  }
}
