#include <assert.h>
#include "contour.h"
#include "readdtcode.h"
#include "fundamental.h"
#include "wirtinger.h"

struct presentation *
wirtingerfromloiv (struct vecofintlist *loiv)
{
  int i, selflink, numnodes, numlabels;
  int over, label, signature;
  extern int verbose;
  int *dtcode, *dt_involution, *dt_realization, *gregionsign;
  int *overgen, *ingen;
  int underpasslabel, overpasslabel, nextnode, nextlabel, nextunder;

  printf ("Work in progress...\n");
  if (verbose) printloiv (loiv);
  assert (loiv->next == 0);
  assert (loiv->handedness);
  dtcode = loiv->vec;
  numnodes = loiv->len;
  numlabels = 2*numnodes;
  gregionsign = loiv->handedness;

  dt_involution = (int *) malloc (numlabels * sizeof (int));
  dt_realization = (int *) malloc (numlabels * sizeof (int));
  overgen = (int *) malloc (numnodes * sizeof (int));
  ingen = (int *) malloc (numnodes * sizeof (int));

  selflink = 0;
  for (i = 0; i < numnodes; i++)
  {
    label = abs(dtcode[i]) - 1;
    dt_involution[2*i] = label;
    dt_involution[label] = 2*i;
    dt_realization[2*i] = gregionsign[i];
    dt_realization[dt_involution[2*i]] = -gregionsign[i];
  }

  for (i = 0; i < numnodes; i++)
  {
    underpasslabel = dt_involution[2*i];
    if (dtcode[i] < 0) underpasslabel = 2*i;
    /* generator of exiting arc has the same index as the node */
    nextlabel = underpasslabel;
    while (1)
    {
      nextlabel++;
      if (nextlabel >= numlabels) nextlabel = 0;
      nextnode = nextlabel/2;
      if ( (nextlabel % 2) == 1 ) nextnode = dt_involution[nextlabel]/2;
      nextunder = dt_involution[2*nextnode];
      //nextover = 2*nextnode;
      if (dtcode[nextnode] < 0)
      {
        nextunder = 2*nextnode;
        //nextover = dt_involution[2*nextnode];
      }
      if (nextunder == nextlabel)
      {
        ingen[nextnode] = i;
        break;
      }
      overgen[nextnode] = i;
    }
  }

  for (i = 0; i < numnodes; i++)
  {
    underpasslabel = dt_involution[2*i];
    if (dtcode[i] < 0) underpasslabel = 2*i;
    overpasslabel = dt_involution[underpasslabel];
    if (verbose) printf ("Node: %d involving labels %d %d\n", i + 1, 2*i + 1, dt_involution[2*i] + 1);
    over = 1;
    if (dtcode[i] < 0) over = -1;
    if (verbose) printf (" overpass label: %d\n", overpasslabel + 1);
    if (verbose) printf (" overpass generator: %d\n", overgen[i] + 1);
    if (verbose) printf (" incoming generator: %d\n", ingen[i] + 1);
    if (verbose) printf (" outgoing generator: %d\n", i+1);
    if (dt_realization[underpasslabel] > 0)
    {
      printf (" relator: %c = ", i + 'a');
      printf ("%c", overgen[i] + 'A');
      printf ("%c", ingen[i] + 'a');
      printf ("%c", overgen[i] + 'a');
      printf ("\n");
    } else {
      printf (" relator: %c = ", i + 'a');
      printf ("%c", overgen[i] + 'a');
      printf ("%c", ingen[i] + 'a');
      printf ("%c", overgen[i] + 'A');
      printf ("\n");
    }
    signature = over*gregionsign[i];
    assert (signature);
    selflink += signature;
    if (verbose) printf (" selflinking: %c\n", (signature > 0)?'+':'-');
  }
  if (verbose) printf ("Selflinking number: %d\n", selflink);

  free (dt_involution);
  free (dt_realization);
  free (overgen);
  free (ingen);
  return (0);
}
