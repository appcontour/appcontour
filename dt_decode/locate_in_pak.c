#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

void printcode (int codelen, int dtcode[], int dtpositive[]);

int
main (int argc, char *argv[])
{
  int start = 0;
  int end = INT_MAX;
  int nodeid = 0;
  int i, crossings, codelen, alternate;
  char path[]="./";
  char ch, filename[80];
  FILE *pakfile;
  int *dtcode, *dtpositive;
  int tailstart, high, low;

  if (argc < 3)
  {
    printf ("Usage: %s crossings a|n [startnode [endnode]]\n", argv[0]);
    exit (1);
  }
  crossings = atoi (argv[1]);
  assert (crossings >= 3 && crossings <= 16);
  dtcode = (int *) malloc (crossings *sizeof (int));
  dtpositive = (int *) malloc (crossings *sizeof (int));
  alternate = -1;
  if (strcmp (argv[2], "a") == 0) alternate = 1;
  if (strcmp (argv[2], "n") == 0) alternate = 0;
  if (alternate < 0)
  {
    printf ("invalid alternate value: %s, valid values: a or n\n", argv[2]);
    exit (2);
  }
  if (argc > 3)
  {
    start = atoi (argv[3]);
    if (argc > 4)
    {
      end = atoi (argv[4]);
    } else {
      end = start;
    }
  }
  sprintf (filename, "%s%d%c.pak", path, crossings, (alternate)?'a':'n');

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

  printf ("[");
  for (i = 0; i < codelen; i++)
  {
    sign = 2*dtpositive[i]-1;
    printf ("%d ", sign*(2*dtcode[i] + 2));
    sum1 += i;
    sum2 += dtcode[i];
  }
  sum1 += codelen;
  sign = 2*dtpositive[codelen]-1;
  printf ("%d]\n", sign*(2*(sum1 - sum2) + 2));
}
