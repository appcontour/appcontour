#include "morse2morse.h"

/*
 * convert "simple-arcs" morse description to "extended-arcs" morse
 * description
 */

int
main (int argc, char *argv[])
{
  int i;
  FILE *file;
  char *optp;
  char *descfile = 0;
  struct mdesc *morsedesc;

  for (i = 1; i < argc; i++)
  {
    if (*argv[i] == '-')
    {
      optp = argv[i];
      for (optp = argv[i] + 1; *optp; optp++)
      {
        switch (*optp)
        {
          default:
            fprintf (stderr, "Invalid option %c\n", *optp);
            usage (argc, argv);
            exit (1);
        }
      }
    } else {
      if (descfile)
      {
        fprintf (stderr, "Only one argument allowed\n");
        usage (argc, argv);
        exit (1);
      }
      descfile = argv[i];
    }
  }

  if (descfile == 0)
     file = stdin;
    else file = fopen (descfile, "r");

  if (file == 0)
  {
    perror ("Reading description file");
    exit (2);
  }
  morsedesc = read_contour (file);
  extend_orientations (morsedesc);

  if (check_if_all_oriented (morsedesc) == 0)
  {
    fprintf (stderr, "Insufficient orienting information\n");
    exit (5);
  }

  collect_simple_arcs (morsedesc);
  write_contour (morsedesc);
  exit (0);
}

void
usage (int argc, char *argv[])
{
  printf ("Usage: %s [<file>]\n", argv[0]);
}

void
ll_pa (struct event *ev)
{
  printf ("[%d:%d]", ev->line->tag, ev->tag);
}

void
dump (struct mdesc *contour)
{
  struct line *l;
  struct event *ev;
  int udori, udori2, lrori;

  for (l = contour->lines; l; l = l->next)
  {
    printf ("Line %d:\n", l->tag);
    for (ev = l->events; ev; ev = ev->next)
    {
      printf (" Event %d (%p): ", ev->tag, ev);
      udori = lrori = udori2 = '?';
      if (ev->orientation == ORIENT_UP) udori = 'u';
      if (ev->orientation == ORIENT_DOWN) udori = 'd';
      if (ev->orientation == ORIENT_LEFT) lrori = 'l';
      if (ev->orientation == ORIENT_RIGHT) lrori = 'r';
      if (ev->orientation2 == ORIENT_UP) udori2 = 'u';
      if (ev->orientation2 == ORIENT_DOWN) udori2 = 'd';
      switch (ev->type)
      {
        case EVENT_VERT:
        printf ("VERT %c%d", udori, ev->huffman);
        break;

        case EVENT_MAX:
        printf ("MAX  %c%d", lrori, ev->huffman);
        break;

        case EVENT_MIN:
        printf ("MIN  %c%d", lrori, ev->huffman);
        break;

        case EVENT_CUSP:
        printf ("CUSP %c%d%c", udori,
           ev->huffman,
           (ev->cuspsign==1)?'+':'-');
        break;

        case EVENT_CROSS:
        printf ("CROSS %c%d", udori, ev->huffman);
        printf ("%c%d", udori2, ev->huffman2);
        break;

        default: assert (0);
      }
      printf ("\n");
    }
  }
}
