#include "visible.h"

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
        fprintf (stderr, "Only one argument required\n");
        usage (argc, argv);
        exit (1);
      }
      descfile = argv[i];
      //printf ("arg: %s\n", descfile);
    }
  }

  if (descfile == 0)
  {
    usage (argc, argv);
    exit (1);
  }
  file = fopen (descfile, "r");
  if (file == 0)
  {
    perror ("Reading description file");
    exit (2);
  }
  morsedesc = read_vis_contour (file);
  //printf ("Extending external region:\n");
  extend_external_region (morsedesc);
  //printf ("Extending orientations:\n");
  extend_orientations (morsedesc);
  orientation_from_external_region (morsedesc);
  extend_orientations (morsedesc);

  //  //printf ("Connecting events:\n");
  //  // this should no-longer be necessary!
  //  connect_events_to_each_other (morsedesc);
  if (check_if_all_oriented (morsedesc) == 0)
  {
    fprintf (stderr, "Insufficient orienting information\n");
    exit (5);
  }
  if (check_conditions (morsedesc) == 0)
  {
    fprintf (stderr, "Visible contour does not meet all necessary conditions\n");
    exit (3);
  }
  build_contour (morsedesc);
  //dump (morsedesc);
  exit (0);
}

void
usage (int argc, char *argv[])
{
  printf ("Usage: %s <file>\n", argv[0]);
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

  for (l = contour->lines; l; l = l->next)
  {
    printf ("Line %d:\n", l->tag);
    for (ev = l->events; ev; ev = ev->next)
    {
      printf (" Event %d (%p): ", ev->tag, ev);
      switch (ev->type)
      {
        case EVENT_VERT:
        printf ("VERT %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c - ", ev->external_on_left?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
#endif
        break;

        case EVENT_MAX:
        printf ("MAX  %c - ", (ev->orientation==ORIENT_LEFT)?'l':'r');
        printf ("%c|%c|%c - ", ev->external_on_left?'e':'i',
          ev->external_inside?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
#endif
        break;

        case EVENT_MIN:
        printf ("MIN  %c - ", (ev->orientation==ORIENT_LEFT)?'l':'r');
        printf ("%c|%c|%c - ", ev->external_on_left?'e':'i',
          ev->external_inside?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
#endif
        break;

        case EVENT_EPTOP:
        printf ("EPTOP %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c - ", ev->external_on_left?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a2);
#endif
        break;

        case EVENT_EPBOT:
        printf ("EPBOT %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c - ", ev->external_on_left?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
#endif
        break;

        case EVENT_TJNE:
        printf ("TJNE %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c|%c - ", ev->external_on_left?'e':'i',
          ev->external_inside?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
        ll_pa (ev->a3);
#endif
        break;

        case EVENT_TJNW:
        printf ("TJNW %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c|%c - ", ev->external_on_left?'e':'i',
          ev->external_inside?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
        ll_pa (ev->a3);
#endif
        break;

        case EVENT_TJSE:
        printf ("TJSE %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c|%c - ", ev->external_on_left?'e':'i',
          ev->external_inside?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
        ll_pa (ev->a3);
#endif
        break;

        case EVENT_TJSW:
        printf ("TJSW %c - ", (ev->orientation==ORIENT_UP)?'u':'d');
        printf ("%c|%c|%c - ", ev->external_on_left?'e':'i',
          ev->external_inside?'e':'i', ev->external_on_right?'e':'i');
        printf ("links: ");
#ifdef WITHCONNECTIONS
        ll_pa (ev->a1);
        ll_pa (ev->a2);
        ll_pa (ev->a3);
#endif
        break;

        default: assert (0);
      }
      printf ("\n");
    }
  }
}
