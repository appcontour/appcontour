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
    }
  }
}
