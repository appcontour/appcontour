#include <assert.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#ifdef HAVE_UNISTD_H
  #include <unistd.h>
  #include <sys/types.h>
  #include <sys/wait.h>
#endif

#include "contour.h"
#include "parser.h"
#include "mendes.h"

#define ACTION_NONE 0
#define ACTION_PRINTSKETCH 1
#define ACTION_COMPARE 2
#define ACTION_ISCONTOUR 3
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
#define ACTION_MENDES 17
#define ACTION_ISHUFFMAN 18

int debug = 0;
int quiet = 0;
int verbose = 0;
int heisemberg = -1;
int docanonify = 1;
int dorecomputef = 1;
int doretagregions = 1;
int finfinity = 0;
int mendesge = HGE_TEXT;

int
main (int argc, char *argv[])
{
  struct sketch *sketch, *s1, *s2;
  int i, count, action=ACTION_NONE, res, ccid = 0;
  unsigned int rndseed = 0;
  int newextregion = 0;
  char rule[10];
  FILE *infile = 0;
  struct mendesgraph *mendes;

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
    if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0)
    {
      verbose = 1;
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
    if (strcmp(argv[i],"--seed") == 0)
    {
      rndseed = atoi (argv[++i]);
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
    if (strcmp(argv[i],"--mendes_ge") == 0 ||
        strcmp(argv[i],"--mge") == 0)
    {
      i++;
      if (strcmp(argv[i],"text") == 0) mendesge = HGE_TEXT;
      else if (strcmp(argv[i],"pykig") == 0) mendesge = HGE_PYKIG;
      else if (strcmp(argv[i],"kig") == 0) mendesge = HGE_KIG;
      else
      {
        printf ("Invalid mendes graphic engine selection: %s\n",
                    argv[i]);
        printf ("Valid choices are: text(default), pykig\n");
        exit (111);
      }
      continue;
    }
    if (strcmp(argv[i],"--help") == 0)
    {
      printf ("usage: %s [options] command [file]\n", argv[0]);
      printf ("  possible commands are:\n");
      printf ("  print, printmorse, iscontour, ishuffman, canonify, countcc\n");
      printf ("  testallrules, info, characteristic, knot2morse\n");
      printf ("  compare: lessicographic comparison between two contours, in this\n");
      printf ("    case the stdin (or file) must contain two descriptions\n");
      printf ("  applyrule <rule>: apply indicated rule to contour\n");
      printf ("  extractcc <int>: extract 3D connected component\n");
      printf ("  removecc <int>: remove 3D connected component from contour\n");
      printf ("  leftright: left-right reflection\n");
      printf ("  frontback: front-back reflection\n");
      printf ("  mendes: compute Mendes graph (see Hacon-Mendes-Romero Fuster)\n");
      printf ("  evert <int>: make region <int> become the unbounded region\n");
      printf ("\n  possible options are:\n");
      printf ("  --help: this help\n");
      printf ("  --version: print program version\n");
      printf ("  -q: be quiet\n");
      printf ("  -v|--verbose: be more verbose\n");
      printf ("  --nocanonify: do not canonify region description before printing\n");
      printf ("  --transfer_islands|--ti <int_coded_flags>: information on island\n");
      printf ("      location in case of ambiguity (e.g. rule C2)\n");
      printf ("  --finfinity|--fi <int>: value of f at infinity (default 0)\n");
      printf ("  --seed <int>: initialize random number generator\n");
      printf ("      e.g. for Hacon graph graphic presentation\n");
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
    if (strcmp(argv[i],"isappcon") == 0) action = ACTION_ISHUFFMAN;
    if (strcmp(argv[i],"iscontour") == 0) action = ACTION_ISCONTOUR;
    if (strcmp(argv[i],"ishuffman") == 0) action = ACTION_ISHUFFMAN;
    if (strcmp(argv[i],"compare") == 0) action = ACTION_COMPARE;
    if (strcmp(argv[i],"canonify") == 0) action = ACTION_CANONIFY;
    if (strcmp(argv[i],"knot2morse") == 0) action = ACTION_KNOT2MORSE;
    if (strcmp(argv[i],"printmorse") == 0) action = ACTION_PRINTMORSE;
    if (strcmp(argv[i],"characteristic") == 0) action = ACTION_CHARACTERISTIC;
    if (strcmp(argv[i],"info") == 0) action = ACTION_INFO;
    if (strcmp(argv[i],"frontback") == 0) action = ACTION_FRONTBACK;
    if (strcmp(argv[i],"leftright") == 0) action = ACTION_LEFTRIGHT;
    if (strcmp(argv[i],"mendes") == 0) action = ACTION_MENDES;
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
  srandom (rndseed);
  switch (action)
  {
    case ACTION_PRINTSKETCH:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    // printsketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_APPLYRULE:
    if ((sketch = readcontour (infile)) == 0) exit (14);
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
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (testallrules (sketch) == 0) exit(15);
    break;

    case ACTION_EXTRACTCC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    res = extract_connected_component (ccid, sketch);
    printsketch (sketch);
    if (res == 0) exit (15);
    break;

    case ACTION_REMOVECC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    count = count_connected_components (sketch);
    res = remove_connected_component (ccid, sketch);
    printsketch (sketch);
    if (res == 0) exit (15);
    break;

    case ACTION_COUNTCC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    count = count_connected_components (sketch);
    if (quiet) printf ("%d\n", count);
      else printf ("Connected components: %d\n", count);
    break;

    case ACTION_ISCONTOUR:
    case ACTION_ISHUFFMAN:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    res = appcontourcheck (sketch, 
          (action == ACTION_ISHUFFMAN)?1:0, 
          (quiet)?0:1);
    if (res == 0) exit (15);
    break;

    case ACTION_CANONIFY:
    if ((sketch = readcontour (infile)) == 0) exit (14);
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
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    printmorse (sketch);
    break;

    case ACTION_CHARACTERISTIC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (quiet) printf ("%d\n", euler_characteristic (sketch));
      else printf ("Euler characteristic: %d\n", euler_characteristic (sketch));
    break;

    case ACTION_INFO:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    showinfo (sketch);
    break;

    case ACTION_FRONTBACK:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    frontback (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_LEFTRIGHT:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    leftright (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_EVERT:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    changeextregion (sketch, newextregion);
    dorecomputef = 0;
    if (docanonify) postprocesssketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_MENDES:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    mendes = compute_mendes (sketch);
    print_mendes (mendes);
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
#ifdef HAVE_UNISTD_H
  int retcode, status, cpid;
  int pipedes[2];
  struct sketch *s;
  FILE *tomorsefile;
#endif

  tok = gettoken (file);

  if (tok == TOK_MORSE) return (readmorse (file));
  if (tok == TOK_SKETCH) return (readsketch (file));
#ifdef HAVE_UNISTD_H
  if (tok == TOK_KNOT)
  {
    fprintf (stderr, "reading knot description\n");
    retcode = pipe (pipedes);
    cpid = fork ();
    if (cpid < 0) exit (123);
    if (cpid == 0)  /* this is the child */
    {
      /* the child only writes on pipe */
      close (pipedes[0]);
      /* function knot2morse writes on stdout, changing descriptor */
      retcode = dup2 (pipedes[1], 1);
      if (retcode < 0) {perror ("error in dup2"); exit (222);}
      ungettoken (tok);
      retcode = knot2morse (file);
      close (pipedes[1]);
      if (retcode) exit (0);
      exit (333);
    } else {        /* this is the parent */
      /* the parent reads from the pipe */
      close (pipedes[1]);
      tomorsefile = fdopen (pipedes[0],"r");
      tok = gettoken (tomorsefile);
      if (tok != TOK_MORSE)
      {
        fprintf (stderr, "invalid token %d from knot2morse\n", tok);
        fclose (tomorsefile);
        return (0);
      }
      s = readmorse (tomorsefile);
      fclose (tomorsefile);
      waitpid (cpid, &status, 0);
      retcode = WEXITSTATUS (status);
      return (s);
    }
  }
#endif
  fprintf (stderr, "Only 'morse'/'sketch' formats implemented\n");
  return (0);
}

