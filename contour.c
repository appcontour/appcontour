#include <assert.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "contour.h"
#include "parser.h"

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
#define ACTION_INFO 13
#define ACTION_FRONTBACK 14
#define ACTION_LEFTRIGHT 15
#define ACTION_EVERT 16

int debug = 0;
int quiet = 0;
int heisemberg = -1;
int docanonify = 1;
int dorecomputef = 1;
int finfinity = 0;

int
main (int argc, char *argv[])
{
  struct sketch *sketch, *s1, *s2;
  int i, count, action=ACTION_NONE, res, ccid = 0;
  int newextregion = 0;
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
    if (strcmp(argv[i],"--nocanonify") == 0)
    {
      docanonify = 0;
      continue;
    }
    if (strcmp(argv[i],"--heisemberg") == 0 ||
        strcmp(argv[i],"--ti") == 0 ||
        strcmp(argv[i],"--transfer_islands") == 0)
    {
      heisemberg = atoi (argv[++i]);
      continue;
    }
    if (strcmp(argv[i],"--darkmatter") == 0 ||
        strcmp(argv[i],"--finfinity") == 0 ||
        strcmp(argv[i],"--fi") == 0)
    {
      finfinity = atoi (argv[++i]);
      continue;
    }
    if (strcmp(argv[i],"--help") == 0)
    {
      printf ("usage: %s [options] command [file]\n", argv[0]);
      printf ("  possible commands are:\n");
      printf ("  print, printmorse, isappcon, canonify, countcc\n");
      printf ("  testallrules, info, characteristic, knot2morse\n");
      printf ("  compare: lessicographic comparison between two contours, in this\n");
      printf ("    case the stdin (or file) must contain two descriptions\n");
      printf ("  applyrule <rule>: apply indicated rule to contour\n");
      printf ("  extractcc <int>: extract 3D connected component\n");
      printf ("  removecc <int>: remove 3D connected component from contour\n");
      printf ("\n  possible options are:\n");
      printf ("  --help: this help\n");
      printf ("  --transfer_islands <int_coded_flags>: information on island\n");
      printf ("      location in case of ambiguity (e.g. rule C2)\n");
      printf ("\n  if file is not given, description is taken from standard input\n");
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
    if (strcmp(argv[i],"info") == 0) action = ACTION_INFO;
    if (strcmp(argv[i],"frontback") == 0) action = ACTION_FRONTBACK;
    if (strcmp(argv[i],"leftright") == 0) action = ACTION_LEFTRIGHT;
    if (strcmp(argv[i],"evert") == 0)
    {
      action = ACTION_EVERT;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a region tag\n"); exit (11);}
      newextregion = atoi (argv[i]);
    }
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
    if (docanonify) canonify (sketch);
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
    res = appcontourcheck (sketch, 1);
    if (res == 0) exit (15);
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
    if (docanonify) canonify (sketch);
    printmorse (sketch);
    break;

    case ACTION_CHARACTERISTIC:
    sketch = readcontour (infile);
    if (quiet) printf ("%d\n", euler_characteristic (sketch));
      else printf ("Euler characteristic: %d\n", euler_characteristic (sketch));
    break;

    case ACTION_INFO:
    sketch = readcontour (infile);
    showinfo (sketch);
    break;

    case ACTION_FRONTBACK:
    sketch = readcontour (infile);
    frontback (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_LEFTRIGHT:
    sketch = readcontour (infile);
    leftright (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_EVERT:
    sketch = readcontour (infile);
    changeextregion (sketch, newextregion);
    dorecomputef = 0;
    if (docanonify) postprocesssketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
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

