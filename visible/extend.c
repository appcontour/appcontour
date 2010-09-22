#include "visible.h"

#ifdef WITHCONNECTIONS
void
connect_events_to_each_other (struct mdesc *contour)
{
  int count, i;
  struct line *line, *nl;
  struct event *ev1, *ev2;
  struct event **evdown, **evup;

  for (line = contour->lines; line && line->next; line = line->next)
  {
    nl=line->next;
    count = 0;
    for (ev1 = line->events; ev1; ev1 = ev1->next)
    {
      switch (ev1->type)
      {
        case EVENT_MIN:
        case EVENT_EPBOT:
        break;

        case EVENT_VERT:
        case EVENT_EPTOP:
        case EVENT_TJNE:
        case EVENT_TJNW:
        count++;
        break;

        case EVENT_MAX:
        case EVENT_TJSE:
        case EVENT_TJSW:
        count += 2;
        break;

        default: assert (0);
      }

    }
    evdown = (struct event **) malloc (count*sizeof (struct event *));
    evup = (struct event **) malloc (count*sizeof (struct event *));

    /* fill up evdown */
    for (ev1 = line->events, i=0; ev1; ev1 = ev1->next)
    {
      switch (ev1->type)
      {
        case EVENT_MIN:
        case EVENT_EPBOT:
        break;

        case EVENT_VERT:
        case EVENT_EPTOP:
        case EVENT_TJNE:
        case EVENT_TJNW:
        evdown[i++] = ev1;
        break;

        case EVENT_MAX:
        case EVENT_TJSE:
        case EVENT_TJSW:
        evdown[i++] = ev1;
        evdown[i++] = ev1;
        break;
      }
    }
    assert (i == count);

    /* fill up evup */
    i = 0;
    for (ev2 = nl->events; ev2; ev2 = ev2->next)
    {
      switch (ev2->type)
      {
        case EVENT_MAX:
        case EVENT_EPTOP:
        break;

        case EVENT_VERT:
        case EVENT_EPBOT:
        case EVENT_TJSE:
        case EVENT_TJSW:
        evup[i++] = ev2;
        break;

        case EVENT_MIN:
        case EVENT_TJNE:
        case EVENT_TJNW:
        evup[i++] = ev2;
        evup[i++] = ev2;
        break;
      }
    }
    assert (i == count);

    /* create links down */
    i = 0;
    for (ev1 = line->events; ev1; ev1 = ev1->next)
    {
      switch (ev1->type)
      {
        case EVENT_MIN:
        case EVENT_EPBOT:
        break;

        case EVENT_VERT:
        case EVENT_EPTOP:
        case EVENT_TJNE:
        case EVENT_TJNW:
        ev1->a2 = evup[i++];
        break;

        case EVENT_MAX:
        ev1->a1 = evup[i++];
        ev1->a2 = evup[i++];
        break;

        case EVENT_TJSE:
        ev1->a2 = evup[i++];
        ev1->a3 = evup[i++];
        break;

        case EVENT_TJSW:
        ev1->a3 = evup[i++];
        ev1->a2 = evup[i++];
        break;
      }
    }

    /* create links up */
    i = 0;
    for (ev2 = nl->events; ev2; ev2 = ev2->next)
    {
      switch (ev2->type)
      {
        case EVENT_MAX:
        case EVENT_EPTOP:
        break;

        case EVENT_VERT:
        case EVENT_EPBOT:
        case EVENT_TJSE:
        case EVENT_TJSW:
        ev2->a1 = evdown[i++];
        break;

        case EVENT_MIN:
        ev2->a1 = evdown[i++];
        ev2->a2 = evdown[i++];
        break;

        case EVENT_TJNE:
        ev2->a1 = evdown[i++];
        ev2->a3 = evdown[i++];
        break;

        case EVENT_TJNW:
        ev2->a3 = evdown[i++];
        ev2->a1 = evdown[i++];
        break;
      }
    }
    free (evup);
    free (evdown);
  }
}
#endif

/* --------------------------------------------- */

void
extend_orientations (struct mdesc *contour)
{
  int goon;
  struct line *line;

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (line = contour->lines; line; line = line->next)
    {
      if (line->next) goon += inherit_orientation_lines (line);
    }
  }
}

/* --------------------------------------------- */

int
check_if_all_oriented (struct mdesc *contour)
{
  struct line *line;
  struct event *ev;

  /* check that everything is oriented */
  for (line = contour->lines; line; line = line->next)
  {
    for (ev = line->events; ev; ev = ev->next)
        if (ev->orientation == 0)
        {
          fprintf (stderr, "missing orientation for event %d, line %d\n",
                   ev->tag, line->tag);
          return (0);
        }
  }
  return (1);
}

/* --------------------------------------------- */

void
orientation_from_external_region (struct mdesc *contour)
{
  struct line *line;
  struct event *ev;

  for (line = contour->lines; line; line = line->next)
  {
    for (ev = line->events; ev; ev = ev->next)
    {
      if (ev->orientation) continue;   /* line already oriented */
      switch (ev->type)
      {
        case EVENT_VERT:
        if (ev->external_on_right) ev->orientation = ORIENT_UP;
        if (ev->external_on_left) ev->orientation = ORIENT_DOWN;
        break;

        case EVENT_MAX:
        if (ev->external_on_right) ev->orientation = ORIENT_LEFT;
        if (ev->external_on_left) ev->orientation = ORIENT_LEFT;
        if (ev->external_inside) ev->orientation = ORIENT_RIGHT;
        break;

        case EVENT_MIN:
        if (ev->external_on_right) ev->orientation = ORIENT_RIGHT;
        if (ev->external_on_left) ev->orientation = ORIENT_RIGHT;
        if (ev->external_inside) ev->orientation = ORIENT_LEFT;
        break;

        case EVENT_EPTOP:
        case EVENT_EPBOT:
        /* there is no way to orient endpoints from external region! */
        break;

        case EVENT_TJNE:
        case EVENT_TJSE:
        if (ev->external_on_right) ev->orientation = ORIENT_UP;
        if (ev->external_inside) ev->orientation = ORIENT_DOWN;
        break;

        case EVENT_TJNW:
        case EVENT_TJSW:
        if (ev->external_on_left) ev->orientation = ORIENT_DOWN;
        if (ev->external_inside) ev->orientation = ORIENT_UP;
        break;

        default: assert (0);
      }
    }
  }
}

/* --------------------------------------------- */

/*
 * inherit orientations between two adjacent lines
 */

int
inherit_orientation_lines (struct line *l)
{
  struct line *nl;
  int i, count, dangdown, dangup, *ordown, *orup;
  int or1, or2;
  int prod, sum;
  struct event *ev1, *ev2;

  nl = l->next;
  assert (nl);

  dangdown = count_dangling_down (l);
  dangup = count_dangling_up (nl);

  if (dangdown != dangup)
  {
    fprintf (stderr, "Error: nonmatching events between lines %d and %d\n", l->tag, nl->tag);
    fprintf (stderr, "dangling down: %d, dangling up %d\n", dangdown, dangup);
  }
  assert (dangdown == dangup);

  ordown = (int *) malloc (dangdown * sizeof(int));
  orup = (int *) malloc (dangup * sizeof(int));

  i = 0;
  for (ev1 = l->events; ev1; ev1 = ev1->next)
  {
    switch (ev1->type)
    {
      case EVENT_MIN:
      case EVENT_EPBOT:
      break;

      case EVENT_VERT:
      case EVENT_EPTOP:
      ordown[i++] = ev1->orientation;
      break;

      case EVENT_TJNE:
      ordown[i++] = ORIENT_UP;
      break;

      case EVENT_TJNW:
      ordown[i++] = ORIENT_DOWN;
      break;

      case EVENT_MAX:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      ordown[i++] = -ev1->orientation;
      ordown[i++] = ev1->orientation;
      break;

      case EVENT_TJSE:
      ordown[i++] = ORIENT_UP;
      ordown[i++] = ev1->orientation;
      break;

      case EVENT_TJSW:
      ordown[i++] = ev1->orientation;
      ordown[i++] = ORIENT_DOWN;
      break;

      default: assert (0);
    }
  }

  i = 0;
  for (ev2 = nl->events; ev2; ev2 = ev2->next)
  {
//printf ("ev2->type = %d\n", ev2->type);
    switch (ev2->type)
    {
      case EVENT_MAX:
      case EVENT_EPTOP:
      break;

      case EVENT_VERT:
      case EVENT_EPBOT:
      orup[i++] = ev2->orientation;
      break;

      case EVENT_TJSE:
      orup[i++] = ORIENT_UP;
      break;

      case EVENT_TJSW:
      orup[i++] = ORIENT_DOWN;
      break;

      case EVENT_MIN:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      orup[i++] = ev2->orientation;
      orup[i++] = -ev2->orientation;
      break;

      case EVENT_TJNE:
      orup[i++] = ORIENT_UP;
      orup[i++] = ev2->orientation;
      break;

      case EVENT_TJNW:
      orup[i++] = ev2->orientation;
      orup[i++] = ORIENT_DOWN;
      break;

      default: assert (0);
    }
  }

  /* now inherit between the two vectors */
  count = 0;
  for (i = 0; i < dangdown; i++)
  {
    prod = ordown[i]*orup[i];
    sum = ordown[i] + orup[i];
    if (prod < 0)
    {
      fprintf (stderr, "Inconsistent orientation\n");
      exit (7);
    }
    if (prod == 0 && sum != 0)
    {
      count++;
      ordown[i] = orup[i] = sum;
    }
  }

  if (count == 0) {free(orup); free(ordown); return (0);}

  i = 0;
  for (ev1 = l->events; ev1; ev1 = ev1->next)
  {
    switch (ev1->type)
    {
      case EVENT_MIN:
      case EVENT_EPBOT:
      break;

      case EVENT_VERT:
      case EVENT_EPTOP:
      ev1->orientation = ordown[i++];
      break;

      case EVENT_TJNE:
      case EVENT_TJNW:
      i++;
      break;

      case EVENT_MAX:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      or1 = -ordown[i++];
      or2 = ordown[i++];
      assert (or1*or2 >= 0);
      if (or1 != or2) {or2 += or1; count++;}
      ev1->orientation = or2;
      break;

      case EVENT_TJSE:
      i++;
      ev1->orientation = ordown[i++];
      break;

      case EVENT_TJSW:
      ev1->orientation = ordown[i++];
      i++;
      break;

      default: assert (0);
    }
  }

  i = 0;
  for (ev2 = nl->events; ev2; ev2 = ev2->next)
  {
    switch (ev2->type)
    {
      case EVENT_MAX:
      case EVENT_EPTOP:
      break;

      case EVENT_VERT:
      case EVENT_EPBOT:
      ev2->orientation = orup[i++];
      break;

      case EVENT_TJSE:
      case EVENT_TJSW:
      i++;
      break;

      case EVENT_MIN:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      or1 = orup[i++];
      or2 = -orup[i++];
      assert (or1*or2 >= 0);
      if (or1 != or2) {or2 += or1; count++;}
      ev2->orientation = or2;
      break;

      case EVENT_TJNE:
      i++;
      ev2->orientation = orup[i++];
      break;

      case EVENT_TJNW:
      ev2->orientation = orup[i++];
      i++;
      break;

      default: assert (0);
    }
  }
  free (orup);
  free (ordown);
  return (count);
}

/* --------------------------------------------- */

#define NEWALG 1 

void
extend_external_region (struct mdesc *contour)
{
#ifdef NEWALG
  int goon;
  struct line *l;

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (l = contour->lines; l; l = l->next)
    {
      goon += extend_ext_region_line (l);
      if (l->next) goon += extend_ext_region_lines (l);
    }
  }
#else
  int goon;
  struct line *l;
  struct event *ev;

  goon = 1;
  while (goon)
  {
    goon = 0;
    for (l = contour->lines; l; l = l->next)
    {
dump (contour);
      for (ev = l->events; ev; ev = ev->next)
      {
        if (ev->next)
        {
          if (ev->external_on_right + ev->next->external_on_left == 1)
          {goon++; ev->external_on_right = ev->next->external_on_left = 1;}
        }
        switch (ev->type)
        {
          case EVENT_VERT:
          if (ev->external_on_right)
          {
            if (ev->a1->external_on_right == 0)
            {goon++; ev->a1->external_on_right = 1;}
assert(0); /* sbagliato nel caso ev->a2 e' un minimo */
            if (ev->a2->external_on_right == 0)
            {goon++; ev->a2->external_on_right = 1;}
          }
          if (ev->external_on_left)
          {
            if (ev->a1->external_on_left == 0)
            {goon++; ev->a1->external_on_left = 1;}
            if (ev->a2->external_on_left == 0)
            {goon++; ev->a2->external_on_left = 1;}
          }
          break;

          case EVENT_MAX:
          case EVENT_MIN:
          if (ev->external_on_left + ev->external_on_right == 1)
          {goon++; ev->external_on_left = ev->external_on_right = 1;}
          if (ev->external_on_left && ev->a1->external_on_left == 0)
            {goon++; ev->a1->external_on_left = 1;}
          if (ev->external_on_right && ev->a2->external_on_right == 0)
            {goon++; ev->a2->external_on_right = 1;}
          if (ev->external_inside)
          {
            if (ev->a1->external_on_right == 0)
            {goon++; ev->a1->external_on_right = 1;}
            if (ev->a2->external_on_left == 0)
            {goon++; ev->a2->external_on_left = 1;}
          }
          break;

          case EVENT_EPTOP:
          if (ev->external_on_left + ev->external_on_right == 1)
          {goon++; ev->external_on_left = ev->external_on_right = 1;}
          if (ev->external_on_right)
          {
            if (ev->a2->external_on_right == 0 || ev->a2->external_on_left == 0)
            {goon++; ev->a2->external_on_right = ev->a2->external_on_left = 1;}
          }
          break;

          case EVENT_EPBOT:
          if (ev->external_on_left + ev->external_on_right == 1)
          {goon++; ev->external_on_left = ev->external_on_right = 1;}
          if (ev->external_on_right)
          {
            if (ev->a1->external_on_right == 0 || ev->a1->external_on_left == 0)
            {goon++; ev->a1->external_on_right = ev->a1->external_on_left = 1;}
          }
          break;

          case EVENT_TJNE:
          if (ev->external_on_right)
          {
            if (ev->a3->external_on_right * ev->a2->external_on_right == 0)
            {goon++; ev->a3->external_on_right = ev->a2->external_on_right = 1;}
          }
          if (ev->external_on_left)
          {
            if (ev->a1->external_on_left * ev->a2->external_on_left == 0)
            {goon++; ev->a1->external_on_left = ev->a2->external_on_left = 1;}
          }
          if (ev->external_inside)
          {
            if (ev->a1->external_on_right * ev->a3->external_on_left == 0)
            {goon++; ev->a1->external_on_right = ev->a3->external_on_left = 1;}
          }
          break;

          case EVENT_TJNW:
          if (ev->external_on_left)
          {
            if (ev->a3->external_on_left * ev->a2->external_on_left == 0)
            {goon++; ev->a3->external_on_left = ev->a2->external_on_left = 1;}
          }
          if (ev->external_on_right)
          {
            if (ev->a1->external_on_right * ev->a2->external_on_right == 0)
            {goon++; ev->a1->external_on_right = ev->a2->external_on_right = 1;}
          }
          if (ev->external_inside)
          {
            if (ev->a1->external_on_left * ev->a3->external_on_right == 0)
            {goon++; ev->a1->external_on_left = ev->a3->external_on_right = 1;}
          }
          break;

          case EVENT_TJSE:
          if (ev->external_on_right)
          {
            if (ev->a3->external_on_right * ev->a1->external_on_right == 0)
            {goon++; ev->a3->external_on_right = ev->a1->external_on_right = 1;}
          }
          if (ev->external_on_left)
          {
            if (ev->a1->external_on_left * ev->a2->external_on_left == 0)
            {goon++; ev->a1->external_on_left = ev->a2->external_on_left = 1;}
          }
          if (ev->external_inside)
          {
            if (ev->a2->external_on_right * ev->a3->external_on_left == 0)
            {goon++; ev->a2->external_on_right = ev->a3->external_on_left = 1;}
          }
          break;

          case EVENT_TJSW:
          if (ev->external_on_left)
          {
            if (ev->a3->external_on_left * ev->a1->external_on_left == 0)
            {goon++; ev->a3->external_on_left = ev->a1->external_on_left = 1;}
          }
          if (ev->external_on_right)
          {
            if (ev->a1->external_on_right * ev->a2->external_on_right == 0)
            {goon++; ev->a1->external_on_right = ev->a2->external_on_right = 1;}
          }
          if (ev->external_inside)
          {
            if (ev->a2->external_on_left * ev->a3->external_on_right == 0)
            {goon++; ev->a2->external_on_left = ev->a3->external_on_right = 1;}
          }
          break;
        }
      }
    }
  }
#endif
}

/* --------------------------------------------- */

int
extend_ext_region_line (struct line *l)
{
  struct event *ev, *nev;
  int count = 0;
  int sum;

  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MAX:
      case EVENT_MIN:
      case EVENT_EPTOP:
      case EVENT_EPBOT:
      sum = ev->external_on_right + ev->external_on_left;
      if (sum == 1)
      {
        count++;
        ev->external_on_right = ev->external_on_left = 1;
      }
      break;
    }
  }
  for (ev = l->events; ev && ev->next; ev = ev->next)
  {
    nev = ev->next;
    sum = ev->external_on_right + nev->external_on_left;
    if (sum == 1)
    {
      count++;
      ev->external_on_right = nev->external_on_left = 1;
    }
  }
  return (count);
}

/* --------------------------------------------- */

int
extend_ext_region_lines (struct line *l)
{
  struct line *nl;
  struct event *ev;
  int *regvec;
  int dangdown, dangup, i;
  int count = 0;

  nl = l->next;
  assert (nl);
  dangdown = count_dangling_down (l);
  dangup = count_dangling_up (nl);
  assert (dangdown == dangup);

  regvec = (int *) malloc ( (dangdown+1) * sizeof (int) );
  regvec[0] = 1;
  i = 1;
  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MIN:
      case EVENT_EPBOT:
      break;

      case EVENT_VERT:
      case EVENT_EPTOP:
      case EVENT_TJNE:
      case EVENT_TJNW:
      regvec[i++] = ev->external_on_right;
      break;

      case EVENT_MAX:
      case EVENT_TJSE:
      case EVENT_TJSW:
      regvec[i++] = ev->external_inside;
      regvec[i++] = ev->external_on_right;
      break;

      default: assert (0);
    }
  }

  /* now scan the following line and inherit */
  if (nl->events) assert (nl->events->external_on_left);
  i = 1;
  for (ev = nl->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MAX:
      case EVENT_EPTOP:
      break;

      case EVENT_VERT:
      case EVENT_EPBOT:
      case EVENT_TJSE:
      case EVENT_TJSW:
      if (regvec[i] > ev->external_on_right)
      {
        count++;
        ev->external_on_right = 1;
      }
      regvec[i++] = ev->external_on_right;
      break;

      case EVENT_MIN:
      case EVENT_TJNE:
      case EVENT_TJNW:
      if (regvec[i] > ev->external_inside)
      {
        count++;
        ev->external_inside = 1;
      }
      regvec[i++] = ev->external_inside;
      if (regvec[i] > ev->external_on_right)
      {
        count++;
        ev->external_on_right = 1;
      }
      regvec[i++] = ev->external_on_right;
      break;

      default: assert (0);
    }
  }

  /* now scan again the first line and inherit */
  i = 1;
  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MIN:
      case EVENT_EPBOT:
      break;

      case EVENT_VERT:
      case EVENT_EPTOP:
      case EVENT_TJNE:
      case EVENT_TJNW:
      if (regvec[i++] > ev->external_on_right)
      {
        count++;
        ev->external_on_right = 1;
      }
      break;

      case EVENT_MAX:
      case EVENT_TJSE:
      case EVENT_TJSW:
      if (regvec[i++] > ev->external_inside)
      {
        count++;
        ev->external_inside = 1;
      }
      if (regvec[i++] > ev->external_on_right)
      {
        count++;
        ev->external_on_right = 1;
      }
      break;

      default: assert (0);
    }
  }

  return (count);
}

/* --------------------------------------------- */

int
count_dangling_down (struct line *l)
{
  int dangdown = 0;
  struct event *ev;

  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MIN:
      case EVENT_EPBOT:
      break;

      case EVENT_VERT:
      case EVENT_EPTOP:
      case EVENT_TJNE:
      case EVENT_TJNW:
      dangdown++;
      break;

      case EVENT_MAX:
      case EVENT_TJSE:
      case EVENT_TJSW:
      dangdown += 2;
      break;

      default: assert (0);
    }
  }
  return (dangdown);
}

/* --------------------------------------------- */

int
count_dangling_up (struct line *l)
{
  int dangup = 0;
  struct event *ev;

  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MAX:
      case EVENT_EPTOP:
      break;

      case EVENT_VERT:
      case EVENT_EPBOT:
      case EVENT_TJSE:
      case EVENT_TJSW:
      dangup++;
      break;

      case EVENT_MIN:
      case EVENT_TJNE:
      case EVENT_TJNW:
      dangup += 2;
      break;

      default: assert (0);
    }
  }
  return (dangup);
}

/* --------------------------------------------- */

/*
 * check if the conditions of the completion theorem are satisfied
 */

int
check_conditions (struct mdesc *contour)
{
  struct line *l;
  struct event *ev;
  int wrongorientation = 0;
  int externalendpoint = 0;
  int errline, errevent;

  for (l = contour->lines; l; l = l->next)
  {
    for (ev = l->events; ev; ev = ev->next)
    {
      /* check orientation compatibility with the external region
       * and that endpoints do not touch the external region
       */
      switch (ev->type)
      {
        case EVENT_EPTOP:
        case EVENT_EPBOT:
          if (ev->external_on_right)
          {
            fprintf (stderr, "External endpoint found: line %d, event %d\n",
              l->tag, ev->tag);
            errline = l->tag;
            errevent = ev->tag;
            externalendpoint++;
          }
          break;

        case EVENT_VERT:
          if ((ev->orientation == ORIENT_UP && ev->external_on_left) ||
              (ev->orientation == ORIENT_DOWN && ev->external_on_right))
            wrongorientation++;
          break;

        case EVENT_MAX:
          if ((ev->orientation == ORIENT_RIGHT && ev->external_on_right) ||
              (ev->orientation == ORIENT_LEFT && ev->external_inside))
            wrongorientation++;
          break;

        case EVENT_MIN:
          if ((ev->orientation == ORIENT_LEFT && ev->external_on_right) ||
              (ev->orientation == ORIENT_RIGHT && ev->external_inside))
            wrongorientation++;
          break;

        case EVENT_TJNE:
        case EVENT_TJSE:
          if (ev->external_on_left) wrongorientation++;
          if ((ev->orientation == ORIENT_UP && ev->external_inside) ||
              (ev->orientation == ORIENT_DOWN && ev->external_on_right))
            wrongorientation++;
          break;

        case EVENT_TJNW:
        case EVENT_TJSW:
          if (ev->external_on_right) wrongorientation++;
          if ((ev->orientation == ORIENT_DOWN && ev->external_inside) ||
              (ev->orientation == ORIENT_UP && ev->external_on_left))
            wrongorientation++;
          break;
      }
    }
  }
  if (externalendpoint)
    fprintf (stderr, "External endpoint found: line %d, event %d\n",
              errline, errevent);
  if (wrongorientation)
    fprintf (stderr, "Wrong orientation for external boundary\n");

  if (externalendpoint || wrongorientation) return (0);
  return (1);
}
