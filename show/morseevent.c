#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "showcontour.h"
#include "morseevent.h"
#include "../parser.h"

/*
 * maxone = 1 implies that at most one nontrivial morse
 * event is allowed per line. If there is more than one
 * then a new line of morse events is simulated to give
 * to the caller what it wants.
 */

void
getmorseevent (struct morseevent *morseevent, int maxone)
{
  static int countdanglingup = 0;
  static int countdanglingdown = 0;
  static int prevlinedangling = 0;
  static int countspecial = 0;
  static struct morseevent mesaved;
  static int rettran1 = 0;
  static int rettran2 = 0;
  static int retnewrow = 0;
  static int retsaved = 0;

  morseevent->ori = morseevent->ori2 = 0;
  morseevent->arc = morseevent->arc2 = 0;

  if (rettran1 > 0)
  {
    rettran1--;
    morseevent->type = ME_TRAN;
    return;
  }

  if (retnewrow > 0)
  {
    retnewrow--;
    morseevent->type = ME_NEWROW;
    return;
  }

  if (rettran2 > 0)
  {
    rettran2--;
    morseevent->type = ME_TRAN;
    return;
  }

  if (retsaved)
  {
    retsaved--;
    morseevent->type = mesaved.type;
    morseevent->ori = mesaved.ori;
    morseevent->ori2 = mesaved.ori2;
    morseevent->arc = mesaved.arc;
    morseevent->arc2 = mesaved.arc2;
  } else {
    getmorseeventl (morseevent);
  }
  if (maxone == 0) return;

  switch (morseevent->type)
  {
    case ME_TOP:
    case ME_BOT:
    case ME_CROSS:

    if (countspecial > 0)
    {
      //printf ("Warning: more than one special!\n");
      /* must simulate another row! */
      /* 1. push this morse event
       * 2. return a number of ME_TRAN equal to
       *    prevlinedangling - countdanglingup
       * 3. return a ME_NEWROW
       * 4. return a number of ME_TRAN equal to countdanglingdown
       */
      mesaved.type = morseevent->type;
      mesaved.ori = morseevent->ori;
      mesaved.ori2 = morseevent->ori2;
      //mesaved.cusps = morseevent->cusps;
      //mesaved.cusps2 = morseevent->cusps2;
      mesaved.arc = morseevent->arc;
      mesaved.arc2 = morseevent->arc2;
      rettran1 = prevlinedangling - countdanglingup;
      retnewrow = 1;
      rettran2 = countdanglingdown;
      retsaved = 1;
      countspecial = 0;
      getmorseevent (morseevent, 1);
      return;
    }
  }
//  printf ("ME: %d\n", morseevent->type);
  switch (morseevent->type)
  {
    case ME_TRAN:
      countdanglingup++;
      countdanglingdown++;
      break;

    case ME_TOP:
      countdanglingdown += 2;
      countspecial++;
      break;

    case ME_BOT:
      countdanglingup += 2;
      countspecial++;
      break;

    case ME_CROSS:
      countdanglingdown += 2;
      countdanglingup += 2;
      countspecial++;
      break;

    case ME_NEWROW:
    case ME_LASTROW:
      assert (prevlinedangling == countdanglingup);
      prevlinedangling = countdanglingdown;
      countdanglingup = countdanglingdown = countspecial = 0;
      break;
  }
}

void
getmorseeventl (struct morseevent *morseevent)
{
  int tok;

  /* ho gia letto la graffa aperta */
  tok = gettokens (stdin);

  switch (tok)
  {
    case TOK_SEMICOLON:
      tok = gettokens (stdin);
      if (tok == TOK_RBRACE)
      {
        morseevent->type = ME_LASTROW;
      } else {
        morseevent->type = ME_NEWROW;
      }
      ungettoken (tok);
      break;

    case KEY_HAT:
    case KEY_A:
      morseevent->type = ME_TOP;
      getarcinfo (morseevent);
      break;

    case KEY_U:
    case KEY_V:
    case KEY_UNDERSCORE:
      morseevent->type = ME_BOT;
      getarcinfo (morseevent);
      break;

    case KEY_SLASH:
    case KEY_BSLASH:
    case KEY_BACKQUOTE:
    case TOK_LPAREN:
    case TOK_RPAREN:
    case KEY_PIPE:
    case KEY_I:
      morseevent->type = ME_TRAN;
      getarcinfo (morseevent);
      break;

    case KEY_X:
      morseevent->type = ME_CROSS;
      getarcinfo (morseevent);
      break;
  }
  return;
}

void
getarcinfo (struct morseevent *morseevent)
{
  getoricusps (&morseevent->ori, &morseevent->arc);
  if (abs(morseevent->ori) == 2) morseevent->ori /= 2;    /* significa speficicato left o right */
  if (morseevent->type == ME_CROSS)
  {
    getoricusps (&morseevent->ori2, &morseevent->arc2);
    if (abs(morseevent->ori) == 2) morseevent->ori /= -2;  /* significa specificato left o right */
  }
}

void
getoricusps (int *oript, struct arc **arcpt)
{
  struct arc *arc = 0;
  int i, tok, prevd;
  int require_rbr = 1;
  int depthind = 0;
  int dbuffer[200];

  assert (*arcpt == 0);
  assert (*oript == 0);
  tok = gettokens (stdin);
  if (tok == ISNUMBER || tok == KEY_LEFT ||
      tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN)
  {
    ungettoken (tok);
    tok = TOK_LBRACKET;
    require_rbr = 0;
  }
  if (tok != TOK_LBRACKET)
  {
    ungettoken (tok);
    return;
  }
  tok = gettokens (stdin);
  if (tok == TOK_RBRACKET) return;
  if (tok == KEY_LEFT || tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN)
  {
    if (tok == KEY_LEFT || tok == KEY_DOWN) *oript = 1;
    if (tok == KEY_RIGHT || tok == KEY_UP) *oript = -1;
    if (tok == KEY_LEFT || tok == KEY_RIGHT) *oript *= 2;    /* left/right restituisce +2 o -2 */
    tok = gettokens (stdin);
  }
  if (tok == TOK_COMMA || tok == ISNUMBER)
  {
    if (tok == ISNUMBER) ungettoken (tok);
    prevd = 0;
    while ((tok = gettokens (stdin)) == ISNUMBER ||
            tok == TOK_PLUS || tok == TOK_MINUS)
    {
      switch (tok)
      {
        case ISNUMBER:
        prevd = gettokennumber ();
        dbuffer[depthind++] = prevd;
        break;
        case TOK_PLUS:
        ++prevd;
        dbuffer[depthind++] = prevd;
        break;
        case TOK_MINUS:
        --prevd;
        dbuffer[depthind++] = prevd;
        break;
      }
      if (depthind >= 198)
      {
        fprintf (stderr, "too many cusps of a single arc\n");
        depthind--;
      }
    }
  }
  if (*oript)
  {
    arc = *arcpt = (struct arc *) malloc (sizeof (struct arc));
    arc->cusps = arc->cuspsinserted = 0;
    arc->first = arc->last = arc->loop = 0;
    arc->refcount = 0;
    arc->d = 0;
  }
  if (depthind > 0) 
  {
    assert (arc);
    arc->cusps = depthind - 1;
    arc->d = (int *) malloc (depthind * sizeof (int));
    for (i = 0; i < depthind; i++) arc->d[i] = dbuffer[i];
  }
  if (require_rbr == 0)
  {
    ungettoken (tok);
    tok = TOK_RBRACKET;
  }
  if (tok != TOK_RBRACKET)
  {
    fprintf (stderr, "Error: right paren expected: %d\n", tok);
    return;
  }
}

