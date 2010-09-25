#include "morse2morse.h"

void
extend_orientations (struct mdesc *contour)
{
  int goon;
  struct line *line;

  goon = 1;
  while (goon)
  {
//printf ("extend_orientation: new iteration\n");
    goon = 0;
    for (line = contour->lines; line; line = line->next)
    {
//printf ("extend_orientation: line: %d\n", line->tag);
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
    {
      if (ev->orientation == 0)
      {
        fprintf (stderr, "missing orientation for event %d, line %d\n",
                 ev->tag, line->tag);
        return (0);
      }
      if (ev->type == EVENT_CROSS && ev->orientation2 == 0)
      {
        fprintf (stderr, "missing orientation2 for CROSS event %d, line %d\n",
                 ev->tag, line->tag);
        return (0);
      }
      if (ev->huffman < 0)
      {
        fprintf (stderr, "missing labelling for event %d, line %d\n",
                 ev->tag, line->tag);
        return (0);
      }
      if (ev->type == EVENT_CROSS && ev->huffman2 < 0)
      {
        fprintf (stderr, "missing huffman2 for CROSS event %d, line %d\n",
                 ev->tag, line->tag);
        return (0);
      }
    }
  }
  return (1);
}

/* --------------------------------------------- */

/*
 * inherit orientations between two adjacent lines
 * also inherit huffman labelling...
 */

int
inherit_orientation_lines (struct line *l)
{
  struct line *nl;
  int i, count, dangdown, dangup, *ordown, *orup;
  int *hufdown, *hufup;
  int or1, or2, huf1, huf2;
  int prod, sum, hufmax, hufmin;
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
  hufdown = (int *) malloc (dangdown * sizeof(int));
  hufup = (int *) malloc (dangup * sizeof(int));

  i = 0;
  for (ev1 = l->events; ev1; ev1 = ev1->next)
  {
//printf ("  event %d: %d\n", ev1->tag, ev1->type);
    switch (ev1->type)
    {
      case EVENT_MIN:
      break;

      case EVENT_VERT:
      ordown[i] = ev1->orientation;
      hufdown[i++] = ev1->huffman;
      break;

      case EVENT_CUSP:
      ordown[i] = ev1->orientation;
      hufdown[i] = ev1->huffman;
      if (ev1->orientation == ORIENT_DOWN)
        hufdown[i] += ev1->cuspsign;
      i++;
      break;

/*
 * information orientation and huffman refers to the
 * sw arc
 * information orientation2 and huffman2 refers to the
 * se arc
 */      
      case EVENT_CROSS:
      ordown[i] = ev1->orientation;
      hufdown[i++] = ev1->huffman;
      ordown[i] = ev1->orientation2;
      hufdown[i++] = ev1->huffman2;
      break;

      case EVENT_MAX:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      ordown[i] = -ev1->orientation;
      hufdown[i++] = ev1->huffman;
      ordown[i] = ev1->orientation;
      hufdown[i++] = ev1->huffman;
      break;

      default: assert (0);
    }
  }

  i = 0;
  for (ev2 = nl->events; ev2; ev2 = ev2->next)
  {
//printf ("  eventnl %d: %d\n", ev2->tag, ev2->type);
    switch (ev2->type)
    {
      case EVENT_MAX:
      break;

      case EVENT_VERT:
      orup[i] = ev2->orientation;
      hufup[i++] = ev2->huffman;
      break;

      case EVENT_CUSP:
      orup[i] = ev2->orientation;
      hufup[i] = ev2->huffman;
      if (ev2->orientation == ORIENT_UP)
        hufup[i] += ev2->cuspsign;
      i++;
      break;
      
      case EVENT_CROSS:
      /*
       * we do not propagate huffman information upwords!
       * on the contrary, we trivially inherit what is
       * coming from up in order to prevent increasing the
       * inheritance count which would force another loop!
       */
      orup[i] = ev2->orientation2;
      hufup[i] = hufdown[i];
      i++;
      orup[i] = ev2->orientation;
      hufup[i] = hufdown[i];
      i++;
      break;

      case EVENT_MIN:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      orup[i] = ev2->orientation;
      hufup[i++] = ev2->huffman;
      orup[i] = -ev2->orientation;
      hufup[i++] = ev2->huffman;
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
    assert (prod >= 0);
    if (prod == 0 && sum != 0)
    {
      count++;
      ordown[i] = orup[i] = sum;
    }
    hufmax = hufmin = hufdown[i];
    if (hufmax < hufup[i]) hufmax = hufup[i];
    if (hufmin > hufup[i]) hufmin = hufup[i];
    if (hufmin >= 0) assert (hufmin == hufmax);
    if (hufmax >= 0 && hufmin < 0)
    {
      count++;
      hufdown[i] = hufup[i] = hufmax;
    }
  }

  if (count == 0) {free(orup); free(ordown); free(hufup); free(hufdown); return (0);}

//printf ("information back...\n");
  i = 0;
  for (ev1 = l->events; ev1; ev1 = ev1->next)
  {
    switch (ev1->type)
    {
      case EVENT_MIN:
      break;

      case EVENT_VERT:
      ev1->orientation = ordown[i];
      ev1->huffman = hufdown[i++];
      break;

      case EVENT_CUSP:
      ev1->orientation = ordown[i];
      ev1->huffman = hufdown[i++];
      if (ev1->orientation == ORIENT_DOWN) ev1->huffman -= ev1->cuspsign;
      break;

      case EVENT_MAX:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      or1 = -ordown[i];
      huf1 = hufdown[i++];
      or2 = ordown[i];
      huf2 = hufdown[i++];
      assert (or1*or2 >= 0);
      if (or1 != or2) {or2 += or1; count++;}
      ev1->orientation = or2;
      if (ev1->huffman < 0)
      {
        if (huf1 >= 0)
        {
          ev1->huffman = huf1;
          count++;
          if (huf2 >= 0) assert (huf2 == huf1);
        } else if (huf2 >= 0)
        {
          ev1->huffman = huf2;
          count++;
        }
      }
      break;

      case EVENT_CROSS:
      ev1->orientation = ordown[i];
      ev1->huffman = hufdown[i++];
      ev1->orientation2 = ordown[i];
      ev1->huffman2 = hufdown[i++];
      break;

      default: assert (0);
    }
  }

//printf ("information back (nl)...\n");
  i = 0;
  for (ev2 = nl->events; ev2; ev2 = ev2->next)
  {
    switch (ev2->type)
    {
      case EVENT_MAX:
      break;

      case EVENT_VERT:
      ev2->orientation = orup[i];
      ev2->huffman = hufup[i++];
      break;

      case EVENT_CUSP:
      ev2->orientation = orup[i];
      ev2->huffman = hufup[i++];
      if (ev2->orientation == ORIENT_UP) ev2->huffman -= ev2->cuspsign;
      break;

      case EVENT_MIN:  /* qui e' importante che ORIENT_LEFT == ORIENT_UP */
      or1 = orup[i];
      huf1 = hufup[i++];
      or2 = -orup[i];
      huf2 = hufup[i++];
      assert (or1*or2 >= 0);
      if (or1 != or2) {or2 += or1; count++;}
      ev2->orientation = or2;
      if (ev2->huffman < 0)
      {
        if (huf1 >= 0)
        {
          ev2->huffman = huf1;
          count++;
          if (huf2 >= 0) assert (huf2 == huf1);
        } else if (huf2 >= 0)
        {
          ev2->huffman = huf2;
          count++;
        }
      }
      break;

      case EVENT_CROSS:
      ev2->orientation2 = orup[i++];
      ev2->orientation = orup[i++];
      break;

      default: assert (0);
    }
  }
  free (orup);
  free (ordown);
  free (hufup);
  free (hufdown);
  return (count);
}

/* --------------------------------------------- */

int
count_dangling_down (struct line *l)
{
  int dangdown = 0;
  struct event *ev;

  if (l == 0) return (0);
  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MIN:
      break;

      case EVENT_VERT:
      case EVENT_CUSP:
      dangdown++;
      break;

      case EVENT_MAX:
      case EVENT_CROSS:
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

  if (l == 0) return (0);
  for (ev = l->events; ev; ev = ev->next)
  {
    switch (ev->type)
    {
      case EVENT_MAX:
      break;

      case EVENT_VERT:
      case EVENT_CUSP:
      dangup++;
      break;

      case EVENT_MIN:
      case EVENT_CROSS:
      dangup += 2;
      break;

      default: assert (0);
    }
  }
  return (dangup);
}

