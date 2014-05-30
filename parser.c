#include <stdio.h>
#include <stdlib.h>

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "parser.h"

extern int debug, heisemberg;
static int onecharword = 0;

char
mygetchar (FILE *file)
{
  skipblanks (file);
  return (fgetc (file));
}

int
getword (FILE *file, char *word, int wordmaxlen)
{
  int charcount = 0;
  char ch;

  while (isalnum(ch = fgetc (file)))
  {
    if (onecharword && ! isdigit(ch)) break;
    if (charcount++ >= wordmaxlen)
    {
      fprintf (stderr, "Error: word too long (%c)\n", ch);
      return (TOK_ERROR);
    }
    *word++ = ch;
  }

  if (charcount == 0)
  {
    *word++ = ch;
    *word = 0;
    if (ch == EOF) return (TOK_EOF);
    return (TOK_CHAR);
  }
  ungetc (ch, file);
  *word = 0;
  return (TOK_ID);
}

void
skipblanks (FILE *file)
{
  char ch;

  while ((ch = fgetc (file)) != EOF)
  {
    if (isspace(ch)) continue;
    if (ch == '\n') continue;
    if (ch == '#')
    {
      while ((ch = fgetc (file)) != '\n' && ch != EOF);
      continue;
    }
    break;
  }
  ungetc (ch, file);
}

static int ungetstackptr = 0;
static int ungetstack[100];

void
ungettoken (int tok)
{
  ungetstack[ungetstackptr++] = tok;
}

static char tokenword[80];

int
gettokens (FILE *file)
{
  int flagsaved, tok;

  flagsaved = onecharword;
  onecharword = 1;
  tok = gettoken (file);
  onecharword = flagsaved;
  return (tok);
}

int
gettoken (FILE *file)
{
  if (ungetstackptr)
  {
    return (ungetstack[--ungetstackptr]);
  }
  skipblanks (file);
  if (getword (file, tokenword, 80) == TOK_EOF) return (TOK_EOF);
//printf ("word = %s\n", tokenword);
  if (strcmp(tokenword,"{") == 0) return (TOK_LBRACE);
  if (strcmp(tokenword,"}") == 0) return (TOK_RBRACE);
  if (strcmp(tokenword,"(") == 0) return (TOK_LPAREN);
  if (strcmp(tokenword,")") == 0) return (TOK_RPAREN);
  if (strcmp(tokenword,"[") == 0) return (TOK_LBRACKET);
  if (strcmp(tokenword,"]") == 0) return (TOK_RBRACKET);
  if (strcmp(tokenword,"=") == 0) return (TOK_EQUAL);
  if (strcmp(tokenword,"+") == 0) return (TOK_PLUS);
  if (strcmp(tokenword,"-") == 0) return (TOK_MINUS);
  if (strcmp(tokenword,"A") == 0) return (KEY_A);
  if (strcmp(tokenword,"^") == 0) return (KEY_HAT);
  if (strcmp(tokenword,"V") == 0) return (KEY_V);
  if (strcmp(tokenword,"U") == 0) return (KEY_U);
  if (strcmp(tokenword,"I") == 0) return (KEY_I);
  if (strcmp(tokenword,"X") == 0) return (KEY_X);
  if (strcmp(tokenword,"O") == 0) return (KEY_O);
  if (strcmp(tokenword,"l") == 0) return (KEY_LEFT);
  if (strcmp(tokenword,"r") == 0) return (KEY_RIGHT);
  if (strcmp(tokenword,"u") == 0) return (KEY_UP);
  if (strcmp(tokenword,"d") == 0) return (KEY_DOWN);
  if (strcmp(tokenword,"f") == 0) return (KEY_F);
  if (strcmp(tokenword,"c") == 0) return (KEY_CUSP);
  if (strcmp(tokenword,"<") == 0) return (KEY_LT);
  if (strcmp(tokenword,">") == 0) return (KEY_GT);
//  if (strcmp(tokenword,"`") == 0) return (KEY_NWSE);   /* INCONSISTENT */
//  if (strcmp(tokenword,"'") == 0) return (KEY_NESW);   /* INCONSISTENT */
  if (strcmp(tokenword,"'") == 0) return (KEY_QUOTE);
  if (strcmp(tokenword,"ne") == 0) return (KEY_NE);
  if (strcmp(tokenword,"nw") == 0) return (KEY_NW);
  if (strcmp(tokenword,"se") == 0) return (KEY_SE);
  if (strcmp(tokenword,"sw") == 0) return (KEY_SW);
  if (strcmp(tokenword,";") == 0) return (TOK_SEMICOLON);
  if (strcmp(tokenword,":") == 0) return (TOK_COLON);
  if (strcmp(tokenword,",") == 0) return (TOK_COMMA);
  if (strcmp(tokenword,"/") == 0) return (KEY_SLASH);
  if (strcmp(tokenword,"\\") == 0) return (KEY_BSLASH);
  if (strcmp(tokenword,"`") == 0) return (KEY_BACKQUOTE);
  if (strcmp(tokenword,"|") == 0) return (KEY_PIPE);
  if (strcmp(tokenword,"_") == 0) return (KEY_UNDERSCORE);
  if (strcmp(tokenword,"morse") == 0) return (TOK_MORSE);
  if (strcmp(tokenword,"sketch") == 0) return (TOK_SKETCH);
  if (strcmp(tokenword,"fpgroup") == 0) return (TOK_FPGROUP);
  if (strcmp(tokenword,"alexander") == 0) return (TOK_ALEXANDER);
  if (strcmp(tokenword,"ideal") == 0) return (TOK_IDEAL);
  if (strcmp(tokenword,"knot") == 0) return (TOK_KNOT);
  if (strcmp(tokenword,"tag") == 0) return (TOK_TAG);
  if (strcasecmp(tokenword,"arc") == 0) return (TOK_ARC);
  if (strcasecmp(tokenword,"region") == 0) return (TOK_REGION);
  if (isdigit(tokenword[0])) return (ISNUMBER);
  fprintf (stderr, "Error: unrecognized keyword %s\n", tokenword);
exit (1);
  return (TOK_ERROR);
}

int
gettokennumber ()
{
  return (atoi(tokenword));
}

/*
 * read an expression for a Laurent polynomial in two indeterminates
 */

int
get_unsignednum (FILE *file)
{
  char ch;
  int val = 0;

  skipblanks (file);
  while (isdigit (ch = fgetc (file)))
  {
    val *= 10;
    val += ch - '0';
  }
  ungetc (ch, file);

  return (val);
}

int
get_exponent2 (FILE *file)
{
  char ch;
  int inparens = 0;
  int sign = 1;
  int val;

  ch = mygetchar (file);
  if (ch != '^')
  {
    ungetc (ch, file);
    return (1);
  }

  ch = mygetchar (file);
  if (ch == '(')
  {
    inparens = 1;
    ch = mygetchar (file);
    if (ch == '-') sign = -1;
    if (ch != '-' && ch != '+') ungetc (ch, file);
  } else ungetc (ch, file);

  val = get_unsignednum (file);

  if (inparens)
  {
    ch = mygetchar (file);
    assert (ch == ')');
  }
  return (sign*val);
}

int
get_factor2 (FILE *file, char indet_names[2], int *coefpt, int *exp1pt, int *exp2pt)
{
  char ch;

  *coefpt = 1;
  *exp1pt = *exp2pt = 0;

  ch = mygetchar (file);
  if (isdigit (ch))
  {
    *coefpt = get_unsignednum (file);
    return (1);
  }

  if (islower (ch))
  {
    if (ch == indet_names[0])
    {
      *exp1pt = get_exponent2 (file);
      return (1);
    }
    if (ch == indet_names[1])
    {
      *exp2pt = get_exponent2 (file);
      return (1);
    }
  }
  ungetc (ch, file);
  return (0);
}

int
get_unsignedmonomial2 (FILE *file, char indet_names[2], int *coefpt, int *exp1pt, int *exp2pt)
{
  int coef, exp1, exp2;

  *coefpt = 1;
  *exp1pt = *exp2pt = 0;

  while (get_factor2 (file, indet_names, &coef, &exp1, &exp2))
  {
    *coefpt *= coef;
    *exp1pt += exp1;
    *exp2pt += exp2;
  }
  return (1);
}
