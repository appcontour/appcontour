/*
 * BIG credits to the authors of knotscape for their knotTables
 *
        a        n
------------------
03      1        0
04      1        0
05      2        0
06      3        0
07      7        0
08     18        3
09     41        8
10    123       42
11    367      185
12   1288      888
13   4878     5110
14  19536    27436
15  85263   168030
16 379799  1008906

 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

static int quiet = 0;

void printcode (int codelen, int dtcode[], int dtpositive[]);

int
main (int argc, char *argv[])
{
  int start = 0;
  int end = INT_MAX;
  int nodeid = 0;
  int i, crossings, codelen, alternate;
  char path[]="knotTable/";
  char filename[80];
  char *argpt;
  FILE *pakfile;
  int *dtcode, *dtpositive;
  int ch, tailstart, high, low;

  if (argc >= 2 && strcmp (argv[1], "-q") == 0)
  {
    quiet = 1;
    argc--;
    argv++;
  }

  if (argc < 2)
  {
    printf ("Usage: %s xx[a|n]_nnnnn[-nnnnn]\n", argv[0]);
    exit (1);
  }
  crossings = strtol (argv[1], &argpt, 10);
  assert (crossings >= 3 && crossings <= 16);
  dtcode = (int *) malloc (crossings *sizeof (int));
  dtpositive = (int *) malloc (crossings *sizeof (int));
  switch (*argpt++)
  {
    case 'a':
    alternate = 1;
    break;

    case 'n':
    alternate = 0;
    break;

    default:
    printf ("Invalid character after number of crossings: %c (valid values: 'a' or 'n')\n", *argpt);
    exit (2);
    break;
  }
  if (*argpt == '_')
  {
    argpt++;
    assert (isdigit(*argpt));
    start = strtol (argpt, &argpt, 10);
    if (*argpt == '-')
    {
      argpt++;
      assert (isdigit (*argpt));
      end = strtol (argpt, &argpt, 10);
      assert (*argpt == 0);
    } else {
      end = start;
      assert (*argpt == 0);
    }
  } else {
    assert (*argpt == 0);
  }
  sprintf (filename, "%d%c.pak", crossings, (alternate)?'a':'n');
  if (access (filename, R_OK))
  {
    sprintf (filename, "%s%d%c.pak", path, crossings, (alternate)?'a':'n');
    if (access (filename, R_OK))
    {
      printf ("Cannot open file %s\n", filename);
      exit (3);
    }
  }

  //printf ("filename = %s\n", filename);

  pakfile = fopen (filename, "r");
  codelen = crossings - 1;
  for (i = 0; i < crossings; i++) {dtcode[i] = 0; dtpositive[i] = 1;}
  while ((ch = fgetc (pakfile)) != EOF)
  {
    nodeid++;
    high = ch/16;
    low = ch - 16*high;
    assert (high == low);
    tailstart = codelen - high;
    for (i = tailstart; i < codelen;)
    {
      ch = fgetc (pakfile);
      high = ch/16;
      low = ch - 16*high;
      dtcode[i++] = high;
      if (i >= codelen) break;
      dtcode[i++] = low;
    }
    if (alternate == 0)
    {
      ch = fgetc (pakfile);
      for (i = 7; i >= 0; i--)
      {
        if (i < crossings) dtpositive[i] = ch & 1;
        ch = ch >> 1;
      }
      ch = fgetc (pakfile);
      for (i = 15; i >= 8; i--)
      {
        if (i < crossings) dtpositive[i] = ch & 1;
        ch = ch >> 1;
      }
    }
    if (nodeid >= start && nodeid <= end) printcode (codelen, dtcode, dtpositive);
  }
}

void
printcode (int codelen, int dtcode[], int dtpositive[])
{
  int i, sign;
  int sum1 = 0;
  int sum2 = 0;

  if (!quiet) printf ("dtcode {[");
  for (i = 0; i < codelen; i++)
  {
    sign = 2*dtpositive[i]-1;
    printf ("%d ", sign*(2*dtcode[i] + 2));
    sum1 += i;
    sum2 += dtcode[i];
  }
  sum1 += codelen;
  sign = 2*dtpositive[codelen]-1;
  printf ("%d%s\n", sign*(2*(sum1 - sum2) + 2), (quiet)?"":"]}");
}
