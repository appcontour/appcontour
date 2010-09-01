#include "visible.h"

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
#define TOK_COMMA 12
#define TOK_APEX 13
#define TOK_LEFT 14
#define TOK_RIGHT 15
#define TOK_UP 16
#define TOK_DOWN 17
#define TOK_TJSE 18
#define TOK_TJSW 19
#define TOK_TJNE 20
#define TOK_TJNW 21
#define TOK_SLASH 22
#define TOK_BACKSLASH 23
#define TOK_DOT 24
#define TOK_EXT 25

int tok_stack[TOKENSTACKSIZE];
static int tok_stack_pointer = 0;
static int nowords = 0;

struct mdesc *
read_vis_contour (FILE *file)
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
    line->tag = tag++;
    if (prevline) prevline->next = line;
      else firstline = line;
    prevline = line;
  }
  contour = (struct mdesc *) malloc (sizeof (struct mdesc));
  contour->lines = firstline;
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
    event->external_on_right = event->external_on_left = 0;    /* unknown */
    event->external_inside = 0;
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
      case TOK_APEX:
      event->type = EVENT_EPBOT;
      break;
      case TOK_COMMA:
      event->type = EVENT_EPTOP;
      break;
      case TOK_TJNE:
      event->type = EVENT_TJNE;
      break;
      case TOK_TJNW:
      event->type = EVENT_TJNW;
      break;
      case TOK_TJSE:
      event->type = EVENT_TJSE;
      break;
      case TOK_TJSW:
      event->type = EVENT_TJSW;
      break;
      default:
      fprintf (stderr, "Invalid token for event: %d\n", token);
    }
    token = get_token (file);
    event->orientation = 0;
    switch (token)
    {
      case TOK_LEFT:
      event->orientation = ORIENT_LEFT;
      break;
      case TOK_RIGHT:
      event->orientation = ORIENT_RIGHT;
      break;
      case TOK_UP:
      event->orientation = ORIENT_UP;
      break;
      case TOK_DOWN:
      event->orientation = ORIENT_DOWN;
      break;
      default:
      unget_token (token);
    }
    token = get_token (file);
    switch (token)
    {
      case TOK_EXT:
      event->external_on_right = 1;
      break;
      default:
      unget_token (token);
    }
//printf ("orientation token: %d\n", token);
  }
  if (prevevent) prevevent->external_on_right = 1;
  if (firstevent) firstevent->external_on_left = 1;
  line->events = firstevent;
  return (line);
}

/* ------------------------------------------ */

int
get_tokenl (FILE *file)
{
  char ch, ch2, word[MAXWORDSIZE];

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
    case ',':
      return (TOK_COMMA);
      break;
    case '\'':
      return (TOK_APEX);
      break;
    case '`':
      ch2 = fgetc (file);
      if (ch2 == '/') return (TOK_TJNW);
      ungetc (ch2, file);
      return (TOK_APEX);
      break;
    case '.':
      ch2 = fgetc (file);
      if (ch2 == '\\') return (TOK_TJSW);
      ungetc (ch2, file);
      return (TOK_DOT);
      break;
    case '/':
      ch2 = fgetc (file);
      if (ch2 == '.') return (TOK_TJSE);
      ungetc (ch2, file);
      return (TOK_SLASH);
      break;
    case '\\':
      ch2 = fgetc (file);
      if (ch2 == '\'') return (TOK_TJNE);
      ungetc (ch2, file);
      return (TOK_BACKSLASH);
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
    case 'e':
      return (TOK_EXT);
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

void
unget_token (int token)
{
  assert (tok_stack_pointer < TOKENSTACKSIZE);
  tok_stack[tok_stack_pointer++] = token;
}
