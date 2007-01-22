#include <string.h>
#include <stdlib.h>
#include "grcommon.h"

int motion = 1;
int steps = 10000;

struct grflags grflags;

void
grparser (int *argcpt, char *argv[])
{
  int goon = 1;
  int i, j;

  grflags.xfigspecial = 0;
  grflags.onlyvisible = 0;
  grflags.visiblewidth = 5;
  grflags.dashlength = 8.0;
  grflags.dotspacing = 6.0;

  while (goon)
  {
    goon = 0;
    for (i = 1; i < *argcpt; i++)
    {
      if (strcmp (argv[i], "--pause") == 0)
      {
        motion = 0;
        goon = 1;
        (*argcpt)--;
        for (j = i; j < *argcpt; j++)
        {
          argv[j] = argv[j+1];
        }
      }
      if (strcmp (argv[i], "--steps") == 0)
      {
        steps = atoi (argv[i+1]);
        goon = 1;
        (*argcpt) -= 2;
        for (j = i; j < *argcpt; j++)
        {
          argv[j] = argv[j+2];
        }
      }
      if (strcmp (argv[i], "--xfigspecial") == 0)
      {
        grflags.xfigspecial = 1;
        goon = 1;
        (*argcpt)--;
        for (j = i; j < *argcpt; j++)
        {
          argv[j] = argv[j+1];
        }
      }
      if (strcmp (argv[i], "--onlyvisible") == 0)
      {
        grflags.xfigspecial = 1;
        goon = 1;
        (*argcpt)--;
        for (j = i; j < *argcpt; j++)
        {
          argv[j] = argv[j+1];
        }
      }
    }
  }
}
