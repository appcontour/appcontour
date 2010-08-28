#include <assert.h>
#include <string.h>
#include "contour.h"
#include "parser.h"

extern int debug;

/*
 * possiamo idealmente pensare al tutto come contenuto in un
 * grande S1, che funge da bordo della regione esterna.
 * Quindi ogni riga immaginiamo sia del tipo "I xxx I"
 * con una riga in testa "A" e in fondo "V", che in
 * realta' non ci sono
 */

#define BUFSIZE 100

/* local prototypes */
void readregioninfo (FILE *file, struct border *actregion);
/* end local prototypes */

static int has_huffman_labelling = 0;
static int cusps_as_morse_events = 0;

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
  cusps_as_morse_events = 0;

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
  while ((tok = gettokens (file)) != TOK_RBRACE)
  {
    ungettoken (tok);
    if ((actregionsnum = readrow (file, sketch, actregions, 
                                  actregionsnum, BUFSIZE)) == 0)
    {
      freesketch(sketch);
      return(0);
    }
  }
  if (sketch->regions->next == 0)
  {
    fprintf(stderr, "Warning: empty contour\n");
    sketch->isempty = 1;
    has_huffman_labelling++;
    assert (sketch->regions->border->sponda == infborder);
    assert (sketch->regions->next == 0);
    assert (sketch->regions->border->next == 0);
  }
  assert (actregionsnum == 1);
  assert (infborder->next == infborder);
  assert (infborder->info == &infinity);
  for (region = sketch->regions; region; region = region->next)
  {
    if (region->border->sponda == infborder)
    {
      region->border->sponda = 0;
    }
  }
  free (infborder);
  if (has_huffman_labelling) sketch->huffman_labelling = 1;
  if (debug) printf ("postprocessing...\n");
  postprocesssketch (sketch);
  return (sketch);
}

int
readrow (FILE *file, struct sketch *sketch, 
         struct border *actregions[], int actregionsnum, int vecdim)
{
  int tok, actr = 0, i, oriprod, r1nori;
  int freespace = vecdim - actregionsnum;
  int dincrl, dincrr, orilup, orirup, havecrossinfo;
  struct region *newreg;
  struct arc *narc, *arcleft, *arcright, *ularc, *urarc;
  struct border *b1, *b2, *bleft, *bright, *actborder;
  struct border *r1, *r2, *r3, *rtemp;
  struct region *r3region, *r1region;
  struct borderlist *bl, *r3bl;
  int *newdepths;

  if (debug) printf ("Nuova riga, regioni attive: %d\n", actregionsnum);
  if (debug) checkconsistency (sketch);
  while ((tok = gettokens(file)) != TOK_SEMICOLON)
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

    if (tok == TOK_LBRACE)
    {
      readregioninfo (file, actregions[actr]);
      continue;
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
      b1->orientation = -1*2;
      // b2->region = newreg;
      b2->info = narc;
      b2->orientation = 1*2;
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
        /* prima mi occupo delle orientazioni */
        oriprod = r2->orientation * r2->next->orientation;
	assert (oriprod != 0);
	if (oriprod == -1)
        {
          fprintf (stderr, "fatal: incompatible orientation\n");
          exit (11);  /* comunque c'e' poi un brutto errore di memoria */
        }
        if (oriprod < 0)
        {
          /* must reverse orientation of one arc where
           * the orientation is \pm 2 (not given by user)
           */
          if (abs(r2->orientation) == 2)
          {
            /* changing orientation on the left */
            if (debug) printf ("cambio l'orientazione a sinistra\n");
            r1->next->orientation *= -1;
            r2->orientation *= -1;
          } else {
            /* changing orientation on the right */
            assert (abs(r2->next->orientation) == 2);
            if (debug) printf ("cambio l'orientazione a destra\n");
            r2->next->orientation *= -1;
            r3->orientation *= -1;
          }
          oriprod *= -1;
        }
        if (oriprod == 2)
        {
          /* this means that on exactly one arc the orientation
           * is given by user; extend on the other arc
           */
          if (debug) printf ("propago l'orientazione\n");
          if (abs(r2->orientation) == 2)
          {
            r1->next->orientation /= 2;
            r2->orientation /= 2;
          } else {
            assert (abs(r2->next->orientation) == 2);
            r2->next->orientation /= 2;
            r3->orientation /= 2;
          }
        }
        if (arcleft->depths && arcright->depths)
        {
          if (cusps_as_morse_events)
          { /* in this case cusps are morse events, so we must concatenate the
             * two depth vectors
             */
            assert (r2->orientation * r2->next->orientation == 1);
            if (r2->orientation > 0) newdepths = concatenate_depths (arcleft, arcright);
              else newdepths = concatenate_depths (arcright, arcleft);
            arcleft->cusps += arcright->cusps;
            free (arcright->depths);
            free (arcleft->depths);
            arcright->depths = 0;
            arcleft->depths = newdepths;
            arcleft->depthsdim += arcright->cusps;
            arcleft->dvalues += arcright->cusps;
            arcright->cusps = -1;
          } else {
            // fprintf (stderr, "warning: arc info given twice\n");
            if (arcleft->cusps != arcright->cusps)
            {
              fprintf (stderr, "fatal: incompatible dup definition\n");
              exit (11);
            }
            for (i = 0; i < arcleft->cusps + 1; i++)
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
        }
//        assert (r2->orientation != 0 && r2->next->orientation != 0);
        if (arcleft->depths)        /* rimuovo l'arco di destra */
        {
          if (debug) printf ("rimuovo l'arco di destra\n");
          r2->next->info = r3->info = arcleft;
          removearc (arcright, sketch);
//          if (r2->orientation != r2->next->orientation)
//          {
//            if (debug) printf ("cambio l'orientazione a destra\n");
//            r2->next->orientation *= -1;
//            r3->orientation *= -1;
//          }
        } else {                    /* rimuovo l'arco di sinistra */
          if (debug) printf ("rimuovo l'arco di sinistra\n");
          r1->next->info = r2->info = arcright;
          removearc (arcleft, sketch);
//          if (r2->orientation != r2->next->orientation)
//          {
//            if (debug) printf ("cambio l'orientazione a sinistra\n");
//            r1->next->orientation *= -1;
//            r2->orientation *= -1;
//          }
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
    if (tok == KEY_LT || tok == KEY_GT)
    {
      if (debug) printf ("Trovata chiave %c\n", (tok == KEY_LT)?'<':'>');
      bleft = actregions[actr]->next;

      if (++actr >= actregionsnum)
      {
        fprintf (stderr, 
          "Error: too few active regions at cusp key, ignoring...\n");
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
      r2 = actregions[actr];
      r3 = actregions[actr + 1];
      r1nori = r1->next->orientation;
      ularc = r1->next->info;
      urarc = r3->info;
      orilup = orirup = 0;
      if (r3->orientation < 0) orilup = 1;  /* the  left arc oriented up */
      if (r1nori > 0) orirup = 1;           /* the right arc oriented up */

      dincrl = dincrr = 0;
      havecrossinfo = 1;
      switch (getcrossinfo(file))
      {
        case KEY_NWSE:
          dincrl = 2;
          if (r1nori < 0) dincrl = -2;
          break;
        case KEY_NESW:
          dincrr = 2;
          if (r3->orientation < 0) dincrr = -2;
          break;
        default:
          havecrossinfo = 0;
          break;
      }

      b1 = newborder (r1->border); /* nuovo bordo regione 1 */
      b2 = newborder (newreg->border);
      narc = newarc (sketch);
      b1->next = r1->next;
      r1->next = b1;
      b1->info = b2->info = narc;
      /* inherit orientation from upper-right arc */
      b1->orientation = r2->next->orientation;
      b2->orientation = r3->orientation;
      assert (b1->orientation*b2->orientation < 0);
      // b1->orientation = -1*2;
      // b2->orientation = 1*2;
      newreg->border->sponda = b2;
      actregions[actr++] = b2;
//printborder(r1, r1->region); printf("\n");

      getarcinfo(tok, file, b1, b2);
      if (havecrossinfo) adjustarcinfo (narc, urarc, dincrl, orilup);

      b1 = newborder (newreg->border);
      b2 = newborder (r3->border); /* nuovo bordo regione 1 */
      narc = newarc (sketch);
      actregions[actr] = b2;
      b2->next = r3->next;
      r3->next = b2;
      b1->info = b2->info = narc;
      /* inherit orientation from upper-left arc */
      b1->orientation = r1nori;
      b2->orientation = r2->orientation;
      assert (b1->orientation*b2->orientation < 0);
      // b1->orientation = -1*2;
      // b2->orientation = 1*2;
      b1->next = newreg->border->sponda;
      newreg->border->sponda->next = b1;
//printborder(r3, r3->region); printf("\n");

      getarcinfo(tok, file, b1, b2);
      if (havecrossinfo) adjustarcinfo (narc, ularc, dincrr, orirup);
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

int *
concatenate_depths (struct arc *arcbefore, struct arc *arcafter)
{
  int *newdepths, newdepthsdim, i;

  assert (arcbefore->depths && arcafter->depths);
  newdepthsdim = arcbefore->cusps + arcafter->cusps + 1;
  newdepths = (int *) malloc (newdepthsdim*sizeof(int));
  assert (arcbefore->depths[arcbefore->cusps] == arcafter->depths[0]);
  for (i = 0; i < arcbefore->cusps; i++)
    newdepths[i] = arcbefore->depths[i];
  for (i = 0; i <= arcafter->cusps; i++)
    newdepths[arcbefore->cusps + i] = arcafter->depths[i];

  return (newdepths);
}

int
getcrossinfo (FILE *file)
{
  int tok;

  tok = gettokens (file);
  if (tok == KEY_NWSE || tok == KEY_NESW)
    return (tok);

  ungettoken (tok);
  return (0);
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
  int cusps_no_d = 0, cusp_as_morse = 0;
  int *newdepths;

  if (key == KEY_LT || key == KEY_GT)
  {
    if (debug) printf ("getarcinfo: cusp indicated as morse event\n");
    if (key == KEY_LT) ungettoken (KEY_DOWN);
      else ungettoken (KEY_UP);
    cusps_as_morse_events = 1;
    key = KEY_I;
    cusp_as_morse = 1;
    assert (bleft);
    arc = bleft->info;
    assert (arc);
    if (debug && arc->cusps > 0)
    {
      printf ("Ci sono gia' cuspidi: %d, tag:%d\n", arc->cusps, arc->tag);
    }
  }
  tok = gettokens (file);
  if (tok == ISNUMBER || tok == KEY_CUSP || tok == KEY_LEFT || 
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
  tok = gettokens (file);
  if (tok == TOK_RBRACKET) return (ORIENT_EMPTY);
  if (tok == KEY_LEFT || tok == KEY_RIGHT || tok == KEY_UP || tok == KEY_DOWN)
  {
    orientation = tok;
    tok = gettokens (file);
  }
  if (tok == TOK_COMMA || tok == ISNUMBER || tok == KEY_CUSP || 
      tok == TOK_PLUS || tok == TOK_MINUS)
  {
    if (tok == TOK_PLUS || tok == TOK_MINUS)
    {
      /* must simulate a 0 as a starting number */
      assert (cusps_as_morse_events == 0);  /* incompatible! */
      prevd = depths[depthind++] = 0;
      ungettoken (tok);
    }
    if (tok == ISNUMBER || tok == KEY_CUSP) ungettoken (tok);
    prevd = 0;
    while ((tok = gettokens (file)) == ISNUMBER || 
            tok == TOK_PLUS || tok == TOK_MINUS || tok == KEY_CUSP)
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
        case KEY_CUSP:
        cusps_no_d++;
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
  }
  if (orientation == KEY_RIGHT || orientation == KEY_UP)
  {
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
    printf ("Arc orientation: %d -> %d, depths:", bleft->orientation, orientation);
    for (i = 0; i < depthind; i++)
    {
      printf ("%d ", depths[i]);
    }
    printf ("\n");
  }
  assert (bleft);
  arc = bleft->info;
  assert (arc);
  if (abs(orientation*bleft->orientation) == 1)
  {
    fprintf (stderr, 
             "%s orientation for arc %d\n",
             (orientation*bleft->orientation < 0)?"INCOMPATIBLE":"duplicate",
             arc->tag);
    orientation = bleft->orientation;
    // return (ORIENT_EMPTY);
  }
  if (depthind > 0 && cusps_no_d > 0)
  {
    fprintf (stderr,
             "cannot give cusps information in both ways for arc %d\n",
             arc->tag);
    return (ORIENT_EMPTY);
  }
  if (cusps_as_morse_events && key != KEY_A && key != KEY_X)
  /*
   * note that case KEY_A and KEY_X will allways start new arcs, so
   * we are safe to use the old syntax
   */
  {
    assert (depthind <= 2);
    if (depthind == 1)
    {
      if (key != KEY_V)
      {
        if (orientation > 0)
        {
          assert (depths[0] == arc->depths[0]);
        } else {
          assert (depths[0] == arc->depths[arc->cusps]);
        }
      }
    } else {
      newdepths = (int *) malloc (arc->cusps + 2);
      if (orientation > 0)
      {
        newdepths[0] = depths[0];
        assert (depths[1] == arc->depths[0]);
        for (i = 0; i <= arc->cusps; i++)
        {
          newdepths[i+1] = arc->depths[i];
        }
      } else {
        for (i = 0; i <= arc->cusps; i++)
        {
          newdepths[i] = arc->depths[i];
        }
        assert (arc->depths[arc->cusps] == depths[0]);
        newdepths[arc->cusps + 1] = depths[1];
      }
      free (arc->depths);
      arc->depths = newdepths;
      arc->cusps++;
      arc->dvalues++;
      arc->depthsdim = arc->cusps + 1;
    }
  } else {
    // if (arc->depths) /* <--- it seems there is a memory leak here */
    arc->depths = (int *) malloc (depthind*sizeof(int));
    arc->depthsdim = depthind;
    arc->cusps = depthind - 1;  /* questo e' il valore che fa testo */
    if (arc->cusps < 0) arc->cusps = 0;
    if (depthind > 0)
    {
      /* nota: in caso di loop senza nodi e' necessario   
       * indicare una profondita' in piu' (uguale al valore
       * iniziale) per segnalare il numero giusto di cuspidi
       */
      has_huffman_labelling++;
      for (i = 0; i < depthind; i++)
      {
        arc->depths[i] = depths[i];
      }
    }
  }
  if (cusps_no_d) arc->cusps = cusps_no_d;
  if (bright == 0) orientation = 0;
  if (orientation)
  {
    assert (bright->info == arc);
    assert (bleft->orientation * bright->orientation < 0);
    if (bleft->orientation == -orientation)
    {
      fprintf (stderr, "incompatible orientation\n");
    }
//    if (bleft->orientation == -orientation)
    assert (abs(orientation) == 1);
    bleft->orientation = orientation;
    bright->orientation = -orientation;
//    if (bleft->orientation*orientation < 0)
//    {
//      bleft->orientation *= -1;
//      bright->orientation *= -1;
//    }
  }
  return (orientation);
}

/*
 * adjust the d values on arc in such a way to have a given
 * difference with respect to the previous/following arc
 * oriup = 1 means that we have info between the last d
 * value of arc and the first of prevarc (notation comes
 * from the adjusting after a crossing).
 */

int
adjustarcinfo (struct arc *arc, struct arc *prevarc, int dincr, int oriup)
{
  int refd, arefd, diff, k;

  if (prevarc->depths == 0 || prevarc->depthsdim <= 0)
  {
    fprintf (stderr, 
      "cannot adjust d value across node, info needed of adjacent arc\n");
    return (0);
  }
  if (oriup)
    refd = prevarc->depths[0];
  else
    refd = prevarc->depths[prevarc->cusps];

  refd += dincr;
  if (arc->depths == 0)
  {
    //fprintf (stderr, "arc with no cusps across node...\n");
    arc->depths = (int *) malloc (sizeof (int));
    arc->depthsdim = arc->dvalues = 1;
    arc->cusps = 0;
    arc->depths[0] = 0;
  }
  if (arc->depthsdim <= 0)
  {
    fprintf (stderr, "d value required in adjusting cross info\n");
    return (0);
  }
  if (oriup)
    arefd = arc->depths[arc->cusps];
  else
    arefd = arc->depths[0];

  diff = refd - arefd;
  for (k = 0; k <= arc->cusps; k++) arc->depths[k] += diff;
  //printf ("adjusting: refd = %d, arefd = %d, dincr = %d, oriup = %d\n",
  //       refd, arefd, dincr, oriup);

  return (1);
}

void
readregioninfo (FILE *file, struct border *actregion)
{
  int tok, tag;
  extern int doretagregions;

  while ((tok = gettoken (file)))
  {
    if (tok == TOK_RBRACE) break;
    if (tok == TOK_TAG)
    {
      if ((tok = gettoken (file)) != TOK_COLON)
      {
        ungettoken (tok);
        continue;
      }
      if ((tok = gettoken (file)) != ISNUMBER)
      {
        ungettoken (tok);
        continue;
      }
      tag = gettokennumber ();
      actregion->border->region->tag = tag;
      doretagregions = 0;
    }
  }
  return;
}

