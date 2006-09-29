#include <assert.h>
#include <string.h>
#include "contour.h"

extern int debug;

/*
 * possiamo idealmente pensare al tutto come contenuto in un
 * grande S1, che funge da bordo della regione esterna.
 * Quindi ogni riga immaginiamo sia del tipo "I xxx I"
 * con una riga in testa "A" e in fondo "V", che in
 * realta' non ci sono
 */

#define BUFSIZE 100

struct sketch *
readmorse (FILE *file)
{
  int tok;
  struct sketch *sketch;
  struct region *extregion, *region;
  struct border *actregions[BUFSIZE];
  int actregionsnum;
  struct arc infinity;
  struct border *infborder;
  struct borderlist *bl;

  infinity.tag = INFINITY_ARC;

  sketch = newsketch ();
  extregion = newregion(sketch);
  bl = extregion->border = newborderlist (extregion);
  sketch->regions = extregion;
  actregionsnum = 1;
  actregions[0] = infborder = newborder (bl);
  bl->sponda = actregions[0];
  //extregion->extborder = actregions[0];
  actregions[0]->info = &infinity;
  infinity.regionleft = actregions[0];
  actregions[0]->orientation = 1;

  if (debug) checkconsistency (sketch);
//printsketch (sketch);
  /* leggi una descrizione in formato 'morse theory' */
  tok = gettoken (file);
  if (tok != TOK_LBRACE)
  {
    fprintf (stderr, "Error: left brace expected\n");
    freesketch (sketch);
    return (0);
  }
  while ((tok = gettoken (file)) != TOK_RBRACE)
  {
    ungettoken (tok);
    if ((actregionsnum = readrow (file, sketch, actregions, 
                                  actregionsnum, BUFSIZE)) == 0)
    {
      freesketch(sketch);
      return(0);
    }
  }
  assert (actregionsnum == 1);
  assert (infborder->next == infborder);
  assert (infborder->info == &infinity);
  for (region = sketch->regions; region->next; region = region->next)
  {
    if (region->border->sponda == infborder)
    {
      region->border->sponda = 0;
    }
  }
  free (infborder);
  if (debug) printf ("postprocessing...\n");
  postprocesssketch (sketch);
  return (sketch);
}

int
readrow (FILE *file, struct sketch *sketch, 
         struct border *actregions[], int actregionsnum, int vecdim)
{
  int tok, actr = 0, i;
  int freespace = vecdim - actregionsnum;
  struct region *newreg;
  struct arc *narc, *arcleft, *arcright;;
  struct border *b1, *b2, *bleft, *bright, *actborder;
  struct border *r1, *r2, *r3, *rtemp;
  struct region *r3region, *r1region;
  struct borderlist *bl, *r3bl;

  if (debug) printf ("Nuova riga, regioni attive: %d\n", actregionsnum);
  if (debug) checkconsistency (sketch);
  while ((tok = gettoken(file)) != TOK_SEMICOLON)
  {
    switch (tok)
    {
      case KEY_HAT:
      tok = KEY_A;
      break;

      case KEY_U:
      case KEY_UNDERSCORE:
      tok = KEY_V;
      break;

      case KEY_SLASH:
      case KEY_BSLASH:
      case KEY_BACKQUOTE:
      case TOK_LPAREN:
      case TOK_RPAREN:
      case KEY_PIPE:
      tok = KEY_I;
      break;
    }

    if (tok == KEY_A)
    {
      if (debug) printf ("Trovata chiave A\n");
      /* ho bisogno di due ulteriori 'border' e una nuova regione */
      newreg = newregion (sketch);
      newreg->border = newborderlist (newreg);
      actborder = actregions[actr];
      b1 = newborder (actborder->border); /* nuovo bordo regione vecchia */
      b2 = newborder (newreg->border);    /* bordo regione nuova */
      narc = newarc (sketch);
      b1->next = actborder->next;
      actborder->next = b1;
      // b1->region = actborder->region;
      b1->info = narc;
      b1->orientation = -1;
      // b2->region = newreg;
      b2->info = narc;
      b2->orientation = 1;
      newreg->border->sponda = b2;

      freespace -= 2;
      actregionsnum += 2;
      actr++;
      assert (freespace >= 0);
      /* devo inserire una nuova regione attiva */
      for (i = actregionsnum; i > actr; i--)
      {
        actregions[i] = actregions[i-2];
      }
      actregions[actr++] = b2;
      actregions[actr] = b1;

      getarcinfo(tok, file, b1, b2);
      continue;
    }
    if (tok == KEY_V)
    {
      if (debug) printf ("Trovata chiave V\n");
      if (actr + 3 > actregionsnum)
      {
        fprintf (stderr, 
          "Error: too few active regions at V key, ignoring...\n");
        continue;
      }
      freespace += 2;
      actregionsnum -= 2;
      r1 = actregions[actr];
      r2 = actregions[actr+1];
      r3 = actregions[actr+2];
      /*
       * r2 corrisponde ad una regione che ora deve essere chiusa, 
       * mentre r1 e r3 devono essere unite insieme
       */
      for (i = actr+1; i < actregionsnum; i++)
      {
        actregions[i] = actregions[i+2];
      }
      arcleft = r2->info;
      arcright = r2->next->info;
      assert (arcleft == r1->next->info);
      assert (arcright == r3->info);
      bleft = r2->next;  /* serve dopo... */
      bright = r3;
      if (arcleft != arcright) /* devo identificare gli archi */
      {
        if (debug) printf ("identifico gli archi %d e %d\n",
                   arcleft->tag, arcright->tag);
        if (arcleft->depths && arcright->depths)
        {
          fprintf (stderr, "warning: arc info given twice\n");
          if (r2->orientation != r2->next->orientation)
          {
            fprintf (stderr, "fatal: incompatible orientation\n");
            exit (11);  /* comunque c'e' poi un brutto errore di memoria */
          }
          if (arcleft->depthsdim != arcright->depthsdim)
          {
            fprintf (stderr, "fatal: incompatible dup definition\n");
            exit (11);
          }
          for (i = 0; i < arcleft->depthsdim; i++)
          {
            if (arcleft->depths[i] != arcright->depths[i])
            {
              fprintf (stderr, "fatal: incompatible dup definition\n");
              exit (11);
            }
          }
          free (arcright->depths);   
          arcright->cusps = -1;
          arcright->depths = 0;
        }
        assert (r2->orientation != 0 && r2->next->orientation != 0);
        if (arcleft->depths)        /* rimuovo l'arco di destra */
        {
          if (debug) printf ("rimuovo l'arco di destra\n");
          r2->next->info = r3->info = arcleft;
          removearc (arcright, sketch);
          if (r2->orientation != r2->next->orientation)
          {
            if (debug) printf ("cambio l'orientazione a destra\n");
            r2->next->orientation *= -1;
            r3->orientation *= -1;
          }
        } else {                    /* rimuovo l'arco di sinistra */
          if (debug) printf ("rimuovo l'arco di sinistra\n");
          r1->next->info = r2->info = arcright;
          removearc (arcleft, sketch);
          if (r2->orientation != r2->next->orientation)
          {
            if (debug) printf ("cambio l'orientazione a sinistra\n");
            r1->next->orientation *= -1;
            r2->orientation *= -1;
          }
        }
      }
//printsketch(sketch);
      /* chiudo la regione 2 */
      if (r2 != r2->next)  /* altrimenti non devo fare nulla */
      {
        if (debug) printf ("chiudo la regione r2\n");
        rtemp = r2->next;
        ensurecanremoveborder (rtemp);
        r2->next = rtemp->next;
        for (i = 0; i < actregionsnum; i++)
          if (actregions[i] == rtemp) actregions[i] = r2;
        if (bleft == rtemp) bleft = r2;
        if (bright == rtemp) bright = r2;
        free (rtemp);
//printsketch(sketch);
      }
      /* riunisco le regioni r1 e r3 */
      assert (r3->info == r1->next->info);
      r3region = r3->border->region;  /* rischio di perdere r3->region... */
      if (r1->border->region == r3region)
      {
        if (debug) checkconsistency (sketch);
        if (debug) printf ("sono la stessa regione, si crea un buco\n");
        assert (r1->border == r3->border);
        assert (topo_change (r1, r3));
        if (r3 != r3->next)
        {
          rtemp = r3->next;
          ensurecanremoveborder (rtemp);
          r3->next = rtemp->next;
          for (i = 0; i < actregionsnum; i++) /* meglio controllarli tutti! */
            if (actregions[i] == rtemp) actregions[i] = r3;
          if (bleft == rtemp) bleft = r3;
          if (bright == rtemp) bright = r3;
          free (rtemp);
        }
        if (debug) checkconsistency (sketch);
      } else {
        if (debug) checkconsistency (sketch);
        if (debug) printf ("sono regioni diverse, le devo riunire\n");
        /* controllo che r3 non stia in un buco di r3region */
        assert (findinborder (r3, r3region->border->sponda));
        /* riunisco le regioni */
        r1region = regionunion (r1->border->region, r3region, sketch);
        assert (r1region == r1->border->region);
        /* ora devo effettuare il merge dei bordi esterni */
        r3bl = r3->border;
        assert (r3bl != r1->border);
        redefineborder (r3, r1->border);
        rtemp = r1->next;
        r1->next = r3->next;
        r3->next = rtemp;
        if (r3 != r3->next)
        {
          rtemp = r3->next;
          ensurecanremoveborder (rtemp);
          r3->next = rtemp->next;
          if (actregions[actr] == rtemp) actregions[actr] = r3;
          for (i = 0; i < actregionsnum; i++) /* meglio controllarli tutti! */
            if (actregions[i] == rtemp) actregions[i] = r3;
          free (rtemp);
        }
        /* cerco e rimuovo r3bl dalla lista delle c.c. */
        for (bl = r1->border->region->border; bl; bl = bl->next)
        {
          if (bl->next == r3bl)
          {
            bl->next = r3bl->next;
            r3bl->next = 0;
            r3bl->sponda = 0;
            freeborderlist (r3bl);
          }
        }
      }
      getarcinfo(tok, file, bleft, bright);
      if (debug) checkconsistency (sketch);

      continue;
    }
    if (tok == KEY_I)
    {
      if (debug) printf ("Trovata chiave I\n");
      bleft = actregions[actr]->next;

      if (++actr >= actregionsnum)
      {
        fprintf (stderr, 
          "Error: too few active regions at I key, ignoring...\n");
        actr--;
        bright = 0;
      } else bright = actregions[actr];

      getarcinfo(tok, file, bleft, bright);
      continue;
    }
    if (tok == KEY_X)
    {
      if (debug) printf ("Trovata chiave X\n");
      if (actr + 2 >= actregionsnum)
      {
        fprintf (stderr,
          "Error: too few active regions at X key, ignoring...\n");
        continue;
      } 
      newreg = newregion (sketch);
      newreg->border = newborderlist (newreg);
      r1 = actregions[actr++];
      // r2 = actregions[actr++];
      //actregions[actr] = newreg->???;
      r3 = actregions[actr + 1];

      b1 = newborder (r1->border); /* nuovo bordo regione 1 */
      b2 = newborder (newreg->border);
      narc = newarc (sketch);
      b1->next = r1->next;
      r1->next = b1;
      b1->info = b2->info = narc;
      b1->orientation = -1;
      b2->orientation = 1;
      newreg->border->sponda = b2;
      actregions[actr++] = b2;
//printborder(r1, r1->region); printf("\n");

      getarcinfo(tok, file, b1, b2);

      b1 = newborder (newreg->border);
      b2 = newborder (r3->border); /* nuovo bordo regione 1 */
      narc = newarc (sketch);
      actregions[actr] = b2;
      b2->next = r3->next;
      r3->next = b2;
      b1->info = b2->info = narc;
      b1->orientation = -1;
      b2->orientation = 1;
      b1->next = newreg->border->sponda;
      newreg->border->sponda->next = b1;
//printborder(r3, r3->region); printf("\n");

      getarcinfo(tok, file, b1, b2);
      continue;
    }
    fprintf (stderr, "Error: invalid token %d\n", tok);
  }
  if (actregionsnum != actr + 1)
  {
    fprintf (stderr, "Error: mismatch in regions number\n");
  }
  if (debug) printf ("fine riga, regioni attive: %d\n", actregionsnum);
  if (debug) printsketch (sketch);
  if (debug) checkconsistency (sketch);
  return (actregionsnum);
}

int
getarcinfo (int key, FILE *file, 
            struct border *bleft,
            struct border *bright)
{
  int tok, i, prevd;
  int orientation = ORIENT_EMPTY;
  int depths[50];
  int depthind = 0, require_rbr = 1;
  struct arc *arc;

  tok = gettoken (file);
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
    return (ORIENT_EMPTY);
  }
  tok = gettoken (file);
  if (tok == TOK_RBRACKET) return (ORIENT_EMPTY);
  if (tok == KEY_LEFT || tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN)
  {
    orientation = tok;
    tok = gettoken (file);
  }
  if (tok == TOK_COMMA || tok == ISNUMBER)
  {
    if (tok == ISNUMBER) ungettoken (tok);
    prevd = 0;
    while ((tok = gettoken (file)) == ISNUMBER || 
            tok == TOK_PLUS || tok == TOK_MINUS)
    {
      switch (tok)
      {
        case ISNUMBER:
        prevd = depths[depthind++] = gettokennumber ();
        break;
        case TOK_PLUS:
        depths[depthind++] = ++prevd;
        break;
        case TOK_MINUS:
        depths[depthind++] = --prevd;
        break;
      }
    }
  }
  if ((key == KEY_A || key == KEY_V) &&
      (orientation == KEY_UP || orientation == KEY_DOWN))
  {
    orientation = ORIENT_EMPTY;
    fprintf (stderr, 
    "Error: incompatible orientation for A or V key\n");
  }
  if ((key == KEY_I) &&
      (orientation == KEY_RIGHT || orientation == KEY_LEFT))
  {
    orientation = ORIENT_EMPTY;
    fprintf (stderr, 
    "Error: incompatible orientation for A or V key\n");
  }
  if (orientation == KEY_LEFT || orientation == KEY_DOWN)
  {
    orientation = ORIENT_LD;
  } else {
    orientation = ORIENT_RU;
  }
  if (require_rbr == 0)
  {
    ungettoken (tok);
    tok = TOK_RBRACKET;
  } 
  if (tok != TOK_RBRACKET)
  {
    fprintf (stderr, "Error: right paren expected: %d\n", tok);
    return (ORIENT_EMPTY);
  }
  /* setting arc information */
  if (debug)
  {
    printf ("Arc orientation: %d, depths:", orientation);
    for (i = 0; i < depthind; i++)
    {
      printf ("%d ", depths[i]);
    }
    printf ("\n");
  }
  assert (bleft);
  arc = bleft->info;
  assert (arc);
  if (arc->depths)
  {
    fprintf (stderr, 
             "duplicate information for arc %d\n",
             arc->tag);
    return (ORIENT_EMPTY);
  }
  arc->depths = (int *) malloc (depthind*sizeof(int));
  /* nota: in caso di loop senza nodi e' necessario   
   * indicare una profondita' in piu' (uguale al valore
   * iniziale) per segnalare il numero giusto di cuspidi
   */
  // arc->cusps = depthind - 1; (sistemato al postprocess)
  arc->depthsdim = depthind;
  for (i = 0; i < depthind; i++)
  {
    arc->depths[i] = depths[i];
  }
  if (bright == 0) orientation = 0;
  if (orientation)
  {
    assert (bright->info == arc);
    assert (bleft->orientation * bright->orientation < 0);
    if (bleft->orientation == -orientation)
    {
      bleft->orientation *= -1;
      bright->orientation *= -1;
    }
  }
  return (orientation);
}
