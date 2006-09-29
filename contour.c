#include <assert.h>
#include <string.h>
#include "config.h"
#include "contour.h"

#define ACTION_NONE 0
#define ACTION_PRINTSKETCH 1
#define ACTION_COMPARE 2
#define ACTION_ISAPPCON 3
#define ACTION_CANONIFY 4
#define ACTION_APPLYRULE 5
#define ACTION_COUNTCC 6
#define ACTION_EXTRACTCC 7
#define ACTION_TESTALLRULES 8
#define ACTION_KNOT2MORSE 9
#define ACTION_PRINTMORSE 10
#define ACTION_CHARACTERISTIC 11
#define ACTION_REMOVECC 12

int debug = 0;
int quiet = 0;

int
main (int argc, char *argv[])
{
  struct sketch *sketch, *s1, *s2;
  int i, count, action=ACTION_NONE, res, ccid = 0;
  char rule[10];
  FILE *infile = 0;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"--version") == 0)
    {
      printf ("%s version %s\n", argv[0], PACKAGE_VERSION);
      return (0);
    }
    if (strcmp(argv[i],"-q") == 0)
    {
      quiet = 1;
      continue;
    }
    if (strcmp(argv[i],"--debug") == 0)
    {
      debug = 1;
      continue;
    }
    if (strcmp(argv[i],"--help") == 0)
    {
      printf ("usage: %s [print|isappcon|compare|canonify|applyrule|\n",
                                          argv[0]);
      printf ("       testallrules|countcc|extractcc|characteristic|knot2morse]\n");
      exit (0);
    }
    if (infile)
    {
      printf ("Too many arguments: %s\n", argv[i]);
    }
    if (action != ACTION_NONE)
    {
      infile = fopen (argv[i], "r");
      if (infile == 0)
      {
        perror ("Cannot open input file");
        exit (10);
      }
    }
    if (strcmp(argv[i],"applyrule") == 0)
    {
      action = ACTION_APPLYRULE;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a rule\n"); exit (11);}
      fprintf (stderr, "applying rule %s\n", argv[i]);
      strncpy (rule, argv[i], 9);
    }
    if (strcmp(argv[i],"extractcc") == 0)
    {
      action = ACTION_EXTRACTCC;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      ccid = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"removecc") == 0)
    {
      action = ACTION_REMOVECC;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      ccid = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"testallrules") == 0) action = ACTION_TESTALLRULES;
    if (strcmp(argv[i],"countcc") == 0) action = ACTION_COUNTCC;
    if (strcmp(argv[i],"print") == 0) action = ACTION_PRINTSKETCH;
    if (strcmp(argv[i],"isappcon") == 0) action = ACTION_ISAPPCON;
    if (strcmp(argv[i],"compare") == 0) action = ACTION_COMPARE;
    if (strcmp(argv[i],"canonify") == 0) action = ACTION_CANONIFY;
    if (strcmp(argv[i],"knot2morse") == 0) action = ACTION_KNOT2MORSE;
    if (strcmp(argv[i],"printmorse") == 0) action = ACTION_PRINTMORSE;
    if (strcmp(argv[i],"characteristic") == 0) action = ACTION_CHARACTERISTIC;
    if (action == ACTION_NONE)
    {
      fprintf (stderr, "invalid arg[%d] = %s\n", i, argv[i]);
      exit (10);
    }
  }
  if (action == ACTION_NONE) action = ACTION_PRINTSKETCH;
  if (infile == 0) infile = stdin;
  switch (action)
  {
    case ACTION_PRINTSKETCH:
    sketch = readcontour (infile);
    // printsketch (sketch);
    canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_APPLYRULE:
    sketch = readcontour (infile);
    canonify (sketch);
    res = apply_rule (rule, sketch);
    printsketch (sketch);
    if (res == 0)
    {
      fprintf (stderr, "Non trovo alcun match!\n");
      exit(14);
    }
    break;

    case ACTION_TESTALLRULES:
    sketch = readcontour (infile);
    canonify (sketch);
    if (testallrules (sketch) == 0) exit(15);
    break;

    case ACTION_EXTRACTCC:
    sketch = readcontour (infile);
    canonify (sketch);
    res = extract_connected_component (ccid, sketch);
    printsketch (sketch);
    if (res == 0) exit (15);
    break;

    case ACTION_REMOVECC:
    sketch = readcontour (infile);
    canonify (sketch);
    count = count_connected_components (sketch);
    res = remove_connected_component (ccid, sketch);
    printsketch (sketch);
    if (res == 0) exit (15);
    break;

    case ACTION_COUNTCC:
    sketch = readcontour (infile);
    canonify (sketch);
    count = count_connected_components (sketch);
    if (quiet) printf ("%d\n", count);
      else printf ("Connected components: %d\n", count);
    break;

    case ACTION_ISAPPCON:
    sketch = readcontour (infile);
    canonify (sketch);
    appcontourcheck (sketch, 1);
    break;

    case ACTION_CANONIFY:
    sketch = readcontour (infile);
    canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_COMPARE:
    s1 = readcontour (infile);
    canonify (s1);
    s2 = readcontour (infile);
    canonify (s2);
    res = sketchcmp (s1, s2);
    switch (res)
    {
      case 0: printf ("s1 = s2\n");
      break;
      case 1: printf ("s1 > s2\n");
      break;
      case -1: printf ("s1 < s2\n");
      res = 2;
      break;
    }
    exit (res);
    break;

    case ACTION_KNOT2MORSE:
    knot2morse (infile);
    break;

    case ACTION_PRINTMORSE:
    sketch = readcontour (infile);
    //canonify (sketch);
    printmorse (sketch);
    break;

    case ACTION_CHARACTERISTIC:
    sketch = readcontour (infile);
    if (quiet) printf ("%d\n", euler_characteristic (sketch));
      else printf ("Euler characteristic: %d\n", euler_characteristic (sketch));
    break;

    default:
    printf ("invalid action %d\n", action);
    exit (1);
    break;
  }
  if (infile) fclose (infile);
  exit (0);
}

struct sketch *
readcontour (FILE *file)
{
  int tok;
  tok = gettoken (file);

  if (tok == TOK_MORSE) return (readmorse (file));
  if (tok == TOK_SKETCH) return (readsketch (file));
  fprintf (stderr, "Only 'morse'/'sketch' formats implemented\n");
  return (0);
}

