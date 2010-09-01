#include "visible.h"

int
build_contour (struct mdesc *contour)
{
  struct dangarc *darc, *darcn, *darcnn, *darcnnn, darclist;
  struct dangarc *newdarc1, *newdarc2, *idarc, *idarcn, *idarcnn, *idarcnnn;
  struct line *l;
  struct event *ev;
  int fval, targetval, decrease, i, ori;
  int d1, d2;

  printf ("morse {\n");
  /* initial arc list is empty */
  darclist.next = 0;
  for (l = contour->lines; l; l = l->next)
  {
    //print_dang_list (darclist.next, 0);
    darc = &darclist;
    for (ev = l->events; ev; ev = ev->next)
    {
      /* skip invisible dangling arcs */
      darc = skip_invis (darc);
      switch (ev->type)
      {
        case EVENT_VERT:
        darc = darc->next;
        assert (darc && darc->d == 0);
        assert (ev->orientation == darc->orientation);
        break;

        case EVENT_MAX:
        if (ev->orientation == ORIENT_LEFT)
        {
          print_event (EVENT_MAX, ORIENT_LEFT, /* d */ 0, -9999, darclist.next, darc->next);
          newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc1->next = newdarc2;
          newdarc2->next = darc->next;
          darc->next = newdarc1;
          newdarc1->d = newdarc2->d = 0;
          newdarc1->orientation = ORIENT_DOWN;
          newdarc2->orientation = ORIENT_UP;
          darc = newdarc2;
        } else {
          assert (ev->orientation == ORIENT_RIGHT);
          /* figura 6 del paper BBP09b */
          darcn = darc->next;  /* evento successivo */
          fval = compute_f_value (darclist.next, darc);
          targetval = fval - 2;
          if (targetval < 2) targetval = 2;
          if (ev->external_inside) targetval = 0;
          decrease = (fval - targetval)/2 - 1;
          //printf ("fval = %d, targetval = %d, decrease = %d\n", fval, targetval, decrease);
          ori = 1;
          if (decrease < 0) {ori = -1; decrease = -decrease;}
          for (i = 0; i < decrease; i++)
          {
            print_event (EVENT_MAX, ori*ORIENT_RIGHT, 1, -9999, darclist.next, darcn);
            newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
            newdarc1->d = 1;
            newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
            newdarc2->d = 1;
            newdarc1->orientation = ori*ORIENT_UP;
            newdarc2->orientation = ori*ORIENT_DOWN;
            darc->next = newdarc1;
            newdarc1->next = newdarc2;
            newdarc2->next = darcn;
            darc = newdarc1;
            darcn = newdarc2;
          }
          print_event (EVENT_MAX, ORIENT_RIGHT, 0, -9999, darclist.next, darcn);
          newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc1->d = 0;
          newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc2->d = 0;
          newdarc1->orientation = ORIENT_UP;
          newdarc2->orientation = ORIENT_DOWN;
          darc->next = newdarc1;
          newdarc1->next = newdarc2;
          newdarc2->next = darcn;
          darc = newdarc2;
        }
        break;

        case EVENT_EPTOP:
        print_event (EVENT_MAX, ORIENT_LEFT, /* d */ 1, -9999, darclist.next, darc->next);
        newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
        newdarc1->d = 1;
        newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
        newdarc2->d = 0;

        if (ev->orientation == ORIENT_UP)
        {
          newdarc1->next = darc->next;  /* temporarily, to get print_event work properly */
          newdarc2->next = darc->next;
          newdarc1->orientation = ORIENT_DOWN;
          newdarc2->orientation = ORIENT_UP;
        
          darc->next = newdarc1;
          print_event (EVENT_CUSP, ORIENT_UP, /* d */ 0, /* d2 */ 1, darclist.next, newdarc1->next);

          newdarc1->next = newdarc2;
          darc = newdarc2;
        } else {
          darcn = darc->next;
          newdarc1->next = darcn;
          newdarc2->orientation = ORIENT_DOWN;
          newdarc1->orientation = ORIENT_UP;

          darc->next = newdarc1;
          print_event (EVENT_CUSP, ORIENT_DOWN, /* d */ 1, /* d2 */ 0, darclist.next, newdarc1);

          darc->next = newdarc2;
          newdarc2->next = newdarc1;
          darc = newdarc1;
        }
        break;

        case EVENT_EPBOT:
        d1 = 0;
        d2 = 1;
        if (ev->orientation == ORIENT_UP)
        {
          d1 = 1;
          d2 = 0;
        }
        darcn = darc->next;
        darcnn = darcn->next;
        darc->next = darcnn;
        print_event (EVENT_CUSP, ev->orientation, d1, d2, darclist.next, darcnn);
        darc->next = darcn;
        darcn->d = 1;
        darc = darcn;
        break;

        case EVENT_TJSE:
        if (ev->orientation == ORIENT_UP)
        {
          /* caso figura 10 (simmetrico) */
          print_event (EVENT_MAX, ORIENT_LEFT, 2, -9999, darclist.next, darc->next);
          newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc1->d = 2;
          newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc2->d = 0;
          newdarc1->orientation = ORIENT_DOWN;
          newdarc2->orientation = ORIENT_UP;
          darcn = darc->next;
          darcnn = darcn->next;
          newdarc1->next = darcnn;
          newdarc2->next = darcn;
          darc->next = newdarc1;
          print_event (EVENT_XSE, ORIENT_UP, 0, 0, darclist.next, darcnn);
          newdarc1->next = darcn;
          darcn->next = newdarc2;
          newdarc2->next = darcnn;
          darc = newdarc2;
          //print_dang_list (darclist.next, darc);
        } else {
          /* caso figura 11 (simmetrico) */
          darcn = darc->next;    /* ramo NE della tj */
          darcnn = darcn->next;
          fval = compute_f_value (darclist.next, darc);
          targetval = fval - 4;
          if (targetval < 2) targetval = 2;
          if (ev->external_inside) targetval = 0;
          decrease = (fval - targetval)/2 - 2;
          ori = 1;
          if (decrease < 0) {ori = -1; decrease = -decrease;}
          for (i = 0; i < decrease; i++)
          {
            print_event (EVENT_MAX, ori*ORIENT_RIGHT, 3, -9999, darclist.next, darcn);
            newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
            newdarc1->d = 3;
            newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
            newdarc2->d = 3;
            newdarc1->orientation = ori*ORIENT_UP;
            newdarc2->orientation = ori*ORIENT_DOWN;
            darc->next = newdarc1;
            newdarc1->next = darcnn;
            print_event (EVENT_XSE, ori*ORIENT_DOWN, 0, 1, darclist.next, darcnn);
            newdarc2->d = 1;
            darc->next = newdarc1;
            newdarc1->next = darcn;
            darcn->next = newdarc2;
            newdarc2->next = darcnn;
            darc = newdarc1;
            assert (darcn == darc->next);
            darcnn = darcn->next;
          }
          print_event (EVENT_MAX, ORIENT_RIGHT, 2, -9999, darclist.next, darcn);
          newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc1->d = 2;
          newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc2->d = 2;
          newdarc1->orientation = ORIENT_UP;
          newdarc2->orientation = ORIENT_DOWN;
          darc->next = newdarc1;
          newdarc1->next = darcnn;
          print_event (EVENT_XSE, ORIENT_DOWN, 0, 0, darclist.next, darcnn);
          newdarc2->d = 0;
          darc->next = newdarc1;
          newdarc1->next = darcn;
          darcn->next = newdarc2;
          newdarc2->next = darcnn;
          darc = newdarc1;
          assert (darcn == darc->next);
          darcnn = darcn->next;
          darc = newdarc2;
        }
        break;

        case EVENT_TJSW:
        if (ev->orientation == ORIENT_DOWN)
        {
          /* caso figura 10 */
          darcn = darc->next;
          darcnn = darcn->next;
          print_event (EVENT_MAX, ORIENT_LEFT, 2, -9999, darclist.next, darcnn);
          newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc1->d = 2;
          newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc2->d = 2;
          newdarc1->orientation = ORIENT_DOWN;
          newdarc2->orientation = ORIENT_UP;
          darc->next = newdarc2;
          newdarc2->next = darcnn;
          print_event (EVENT_XSW, ORIENT_DOWN, 0, 0, darclist.next, newdarc2);
          newdarc1->d = 0;
          darc->next = newdarc1;
          newdarc1->next = darcn;
          darcn->next = newdarc2;
          darc = newdarc1;
          assert (darcn == darc->next);
          darc = darcn;
        } else {
          /* caso figura 11 */
          darcn = darc->next;    /* ramo NW della tj */
          darcnn = darcn->next;
          fval = compute_f_value (darclist.next, darcn);
          targetval = fval - 4;
          if (targetval < 2) targetval = 2;
          if (ev->external_inside) targetval = 0;
          decrease = (fval - targetval)/2 - 2;
          ori = 1;
          if (decrease < 0) {ori = -1; decrease = -decrease;}
          for (i = 0; i < decrease; i++)
          {
            print_event (EVENT_MAX, ori*ORIENT_RIGHT, 3, -9999, darclist.next, darcnn);
            newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
            newdarc1->d = 3;
            newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
            newdarc2->d = 3;
            newdarc1->orientation = ori*ORIENT_UP;
            newdarc2->orientation = ori*ORIENT_DOWN;
            darc->next = newdarc2;
            newdarc2->next = darcnn;
            print_event (EVENT_XSW, ori*ORIENT_UP, 0, 1, darclist.next, newdarc2);
            newdarc1->d = 1;
            darc->next = newdarc1;
            newdarc1->next = darcn;
            darcn->next = newdarc2;
            darc = newdarc1;
            assert (darcn == darc->next);
            darcnn = darcn->next;
          }
          print_event (EVENT_MAX, ORIENT_RIGHT, 2, -9999, darclist.next, darcnn);
          newdarc1 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc1->d = 2;
          newdarc2 = (struct dangarc *) malloc (sizeof (struct dangarc));
          newdarc2->d = 2;
          newdarc1->orientation = ORIENT_UP;
          newdarc2->orientation = ORIENT_DOWN;
          darc->next = newdarc2;
          newdarc2->next = darcnn;
          print_event (EVENT_XSW, ORIENT_UP, 0, 0, darclist.next, newdarc2);
          newdarc1->d = 0;
          darc->next = newdarc1;
          newdarc1->next = darcn;
          darcn->next = newdarc2;
          darc = newdarc1;
          assert (darcn == darc->next);
          darcnn = darcn->next;
          darc = darcn;
        }
        break;

        case EVENT_TJNE:
	darcn = darc->next;    /* arco NW */
        darcnn = darcn->next;  /* primo dangling invisibile dentro la regione */
        while (darcnn->d > 0)
        {
          darc->next = darcnn->next;
          darcnnn = darcnn->next;
          print_event (EVENT_XNE, darcnn->orientation, 0, darcnn->d, darclist.next, darcnnn);
          darcnn->d += 2;
          darc->next = darcnn;
          darcnn->next = darcn;
          darcn->next = darcnnn;
          darc = darcnn;
          darcnn = darcnnn;
        }
        assert (darcnn->d == 0);
        darc->next = darcnn->next;
        darcnnn = darcnn->next;
        print_event (EVENT_XNE, darcnn->orientation, 0, 0, darclist.next, darcnnn);
        darcnn->d += 2;
        darc->next = darcnn;
        darcnn->next = darcn;
        darcn->next = darcnnn;
        darc = darcn;
        //print_dang_list (darclist.next, darc);
        break;

        case EVENT_TJNW:
	darcn = darc->next;    /* arco NW */
        darcnn = darcn->next;  /* primo dangling invisibile dentro la regione */
        while (darcnn->d > 0)
        {
          /* devo cercare l'ultimo dangling invisibile, 
           * poiche' e' il primo che devo far uscire
           * (figura 9 speculare)
           */
          for (idarc = darcn, idarcn = darcn->next; idarcn->next->d > 0;
                                            idarc = idarc->next, idarcn = idarc->next);
          idarcnn = idarcn->next;
          idarcnnn = idarcnn->next;
          idarc->next = idarcnnn;
          print_event (EVENT_XNW, idarcn->orientation, 0, idarcn->d, darclist.next, idarcnnn);
          idarcn->d += 2;
          idarc->next = idarcnn;
          idarcnn->next = idarcn;
          idarcn->next = idarcnnn;
          assert (darcn == darc->next);
          darcnn = darcn->next;
        }
        assert (darcnn->d == 0);
        darcnnn = darcnn->next;
        darc->next = darcnnn;
        print_event (EVENT_XNW, darcn->orientation, 0, 0, darclist.next, darcnnn);
        darcn->d += 2;
        darc->next = darcnn;
        darcnn->next = darcn;
        darcn->next = darcnnn;
        darc = darcn;
        //print_dang_list (darclist.next, darc);
        break;

        case EVENT_MIN:
        if (ev->orientation == ORIENT_LEFT)
        {
          darcn = darc->next;   /* ramo di sx */
          darcnn = darcn->next; /* primo dangling invisibile dentro la regione */
          while (darcnn->d > 0)
          {
            darc->next = darcnn->next;
            darcnnn = darcnn->next;
            print_event (EVENT_XNE, darcnn->orientation, 0, darcnn->d, darclist.next, darcnnn);
            darcnn->d += 2;
            darc->next = darcnn;
            darcnn->next = darcn;
            darcn->next = darcnnn;
            darc = darcnn;
            darcnn = darcnnn;
          }
          assert (darcnn->d == 0);
          darcnnn = darcnn->next;
          assert (darcn == darc->next);
          darc->next = darcnnn;
          print_event (EVENT_MIN, ORIENT_LEFT, 0, -9999, darclist.next, darcnnn);
          free (darcn);
          free (darcnn);
        } else if (ev->external_on_right)
        {
          /* figure 14; e.g. global minimum */
          darcn = darc->next;   /* descending visible arc */
          fval = compute_f_value (darclist.next, darcn);
          assert (fval == 2);
          while (darcn->next->d > 0)
          {
            /* look for the first upward invisible arc */
            for (idarc = darcn; idarc; idarc = idarc->next)
            {
              assert (idarc->next->d > 0 && idarc->next->next->d > 0);
              if (idarc->next->next->orientation == ORIENT_UP) break;
            }
            assert (idarc);
            /* merge untill there are no more invisible dangling arcs */
            idarcn = idarc->next;
            idarcnn = idarcn->next;
            assert (idarcn->orientation == ORIENT_DOWN);
            assert (idarcnn->orientation == ORIENT_UP);
            idarc->next = idarcnn;
            while (idarcn->d < idarcnn->d)
            {
              print_event (EVENT_CUSP, ORIENT_DOWN, idarcn->d, idarcn->d+1, darclist.next, idarcnn);
              idarcn->d++;
            }
            idarc->next = idarcn;
            idarcn->next = idarcnn->next;
            while (idarcn->d > idarcnn->d)
            {
              print_event (EVENT_CUSP, ORIENT_UP, idarcnn->d+1, idarcnn->d, darclist.next, idarcn->next);
              idarcnn->d++;
            }
            idarcn->next = idarcnn;
            assert (idarcn->d == idarcnn->d);
            /* now we can join the two arcs */
            idarc->next = idarcnn->next;
            print_event (EVENT_MIN, ORIENT_RIGHT, idarcn->d, -9999, darclist.next, idarc->next);
            free (idarcn);
            free (idarcnn);
          }
          darcnn = darcn->next;
          assert (darcnn->d == 0);
          assert (darcnn->orientation == ORIENT_UP);
          darc->next = darcnn->next;
          print_event (EVENT_MIN, ORIENT_RIGHT, 0, -9999, darclist.next, darcnn->next);
          free (darcn);
          free (darcnn);
        } else {
          /* figura 13 */
          darcn = darc->next;   /* ramo di sx */
          darcnn = darcn->next; /* primo dangling invisibile dentro la regione */
          for (idarc = darcn; idarc->next->d > 0; idarc = idarc->next)
          { /* first: raise all d values at least to 2 */
            idarcn = idarc->next;
            if (idarcn->d == 1)
            {
              idarcnn = idarcn ->next;
              idarc->next = idarcnn;
              ori = idarcn->orientation;
              print_event (EVENT_CUSP, ori, (ori==ORIENT_DOWN)?1:2, (ori==ORIENT_DOWN)?2:1,
                 darclist.next, idarcnn);
              idarcn->d = 2;
              idarc->next = idarcn;
            }
          }
          while (darcnn->d > 0)
          {
            /* devo controllare il valore di f a dx di darcnn */
            darcnnn = darcnn->next;
            fval = compute_f_value (darclist.next, darcnn);
            assert (fval >= 2);
            if (fval == 2)
            {
              idarc = darcn;
              idarcn = idarc->next;   /* left arc of the two to join */
              idarcnn = idarcn->next; /* right arc of the two to join */
              assert (idarcn->orientation == ORIENT_UP);
              assert (idarcnn->orientation == ORIENT_DOWN);
              assert (idarcn->d == 2 && idarcnn->d == 2);
              idarcnnn = idarcnn->next;
              idarc->next = idarcnnn;
              print_event (EVENT_MIN, ORIENT_LEFT, 2, -9999, darclist.next, idarcnnn);
              free (idarcn);
              free (idarcnn);
              darcnn = darcn->next;
            } else {
              assert (darcnn->d >= 2);
              if (darcnn->d == 2);
              {
                darcnnn = darcnn ->next;
                darcn->next = darcnnn;
                ori = darcnn->orientation;
                print_event (EVENT_CUSP, ori, (ori==ORIENT_DOWN)?2:3, (ori==ORIENT_DOWN)?3:2,
                   darclist.next, darcnnn);
                darcnn->d = 3;
                darcn->next = darcnn;
              }
              /* ora devo incrociare darcnn con darcn */
              darcnnn = darcnn->next;
              darc->next = darcnnn;
              darcnn->d -= 2;
              print_event (EVENT_XSW, darcnn->orientation, 0, darcnn->d, darclist.next, darcnnn);
              darc->next = darcnn;
              darcnn->next = darcn;
              darcn->next = darcnnn;
              darc = darcnn;
              darcnn = darcnnn;
            }
          }
          //print_dang_list (darclist.next, darc);
          /* non ci sono piu' dangling arcs */
          darcnn = darcn->next;
          assert (darcnn->d == 0);
          darcnnn = darcnn->next;
          darc->next = darcnnn;
          print_event (EVENT_MIN, ORIENT_RIGHT, 0, -9999, darclist.next, darcnnn);
          free (darcn);
          free (darcnn);
        }
        break;

        default:
        fprintf (stderr, "Unrecognized event %d!\n", ev->type);
        exit (10);
      }
    }
    assert (darc);
    // print_dang_list (darclist.next, darc);
    assert (darc->next == 0);
    fval = compute_f_value (darclist.next, 0);
    assert (fval == 0);
  }
  printf ("}\n");
  return (1);
}

/* ------------------------------------------------- */

struct dangarc *
skip_invis (struct dangarc *darc)
{
  while (darc && darc->next && darc->next->d > 0) darc = darc->next;
  return (darc);
}

/* ------------------------------------------------- */
/*
 * scrivi l'evento dopo dangnext
 */

void
print_event (int type, int orientation, int d, int d2,
             struct dangarc *danglist, struct dangarc *dangnext)
{
  struct dangarc *darc;

  for (darc = danglist; darc != dangnext; darc = darc->next)
  {
    printf ("|%c%d ", (darc->orientation == ORIENT_DOWN)?'d':'u', darc->d);
  }
  switch (type)
  {
    case EVENT_MAX:
    printf ("^%c%d ", (orientation == ORIENT_LEFT)?'l':'r', d);
    break;

    case EVENT_MIN:
    printf ("U%c%d ", (orientation == ORIENT_LEFT)?'l':'r', d);
    break;

    case EVENT_CUSP:
    assert (d2 == d + 1 || d2 == d - 1);
    printf ("%c%d%c ", (orientation == ORIENT_DOWN)?'<':'>', d,
      (d2>d)?'+':'-');
    break;

    /*
     * the syntax for crossing for "contour" app is the following:
     *  " X <o1><n1> <o2><n2> "
     * where o1 is the orientation (u/d) of the sw [and ne] arc
     *       n1 is its huffman labelling
     *       o2 is the orientation (u/d) of the se [and nw] arc
     *       n2 is its huffman labelling
     */
    case EVENT_XSE:
    //printf ("X.%c[%d %d] ", (orientation == ORIENT_DOWN)?'d':'u', d, d2);
    printf ("Xu%d%c%d ", d, (orientation == ORIENT_DOWN)?'d':'u', d2);
    break;

    case EVENT_XSW:
    //printf (".X%c[%d %d] ", (orientation == ORIENT_DOWN)?'d':'u', d, d2);
    printf ("X%c%dd%d ", (orientation == ORIENT_DOWN)?'d':'u', d2, d);
    break;

    case EVENT_XNE:
    //printf ("X'%c[%d %d] ", (orientation == ORIENT_DOWN)?'d':'u', d, d2);
    printf ("X%c%du%d ", (orientation == ORIENT_DOWN)?'d':'u', d2+2, d);
    break;

    case EVENT_XNW:
    //printf ("`X%c[%d %d] ", (orientation == ORIENT_DOWN)?'d':'u', d, d2);
    printf ("Xd%d%c%d ", d, (orientation == ORIENT_DOWN)?'d':'u', d2+2);
    break;

    default:
    fprintf (stderr, "Unknown event in print_event: %d\n", type);
    printf ("[??? %d(%d,%d) ???] ", type, d, d2);
    break;
  }
  for (; darc; darc = darc->next)
  {
    printf ("|%c%d ", (darc->orientation == ORIENT_DOWN)?'d':'u', darc->d);
  }
  printf (";\n");
}

/* ------------------------------------------------- */

void
print_dang_list (struct dangarc *list, struct dangarc *target)
{
  struct dangarc *darc;

  printf ("Dangling list: ");
  for (darc = list; darc; darc = darc->next)
  {
    if (darc == target) printf ("*");
    printf ("%c%d ", (darc->orientation==ORIENT_DOWN)?'d':'u', darc->d);
  }
  printf ("\n");
}

/* ------------------------------------------------- */

int
compute_f_value (struct dangarc *list, struct dangarc *target)
{
  struct dangarc *darc;
  int fval;

  fval = 0;
  for (darc = list; darc; darc = darc->next)
  {
    if (darc->orientation == ORIENT_DOWN) fval += 2;
      else fval -= 2;
    if (darc == target) break;
  }
  return (fval);
}
