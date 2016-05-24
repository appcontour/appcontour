#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*
 * the code vector is twice as the number of nodes, we allow
 * for each odd/even label
 */

static int numnodes;
static int numlabels;
static int numregions;
static int *abscode;
static int *dtsign;
static int *regionsign;

void reconstruct_sign (void);
int nextlabel (int label);
int prevlabel (int label);
int inherit (int label);
int isinpath (int i, int start, int end);
int propagate (int sign, int label);
void display_arcs_from_arcs (void);
void display_arcs_from_nodes (void);
void display_regions (void);
void display_regions_from_arcs (void);
void display_regions_from_nodes (void);

int
main (int argc, char *argv[])
{
  int i, curoddnode, curevennode, rcode;

  /*
   * read dt code
   */

  numnodes = argc - 1;
  numlabels = 2*numnodes;
  assert (numnodes >= 2);

  abscode = (int *) malloc ((numlabels+1)*(sizeof (int)));
  dtsign = (int *) malloc ((numlabels+1)*(sizeof (int)));
  regionsign = (int *) malloc ((numlabels+1)*(sizeof (int)));
  for (i = 1; i <= numlabels; i++) regionsign[i] = 0;

  curoddnode = 1;
  for (i = 0; i < numnodes; i++)
  {
    rcode = atoi(argv[i+1]);
    assert ((rcode/2)*2 == rcode);
    abscode[curoddnode] = abs(rcode);
    dtsign[curoddnode] = 1;
    if (rcode < 0) dtsign[curoddnode] = -1;
    abscode[abs(rcode)] = curoddnode;
    curoddnode += 2;
  }
  curevennode = 2;
  for (i = 0; i < numnodes; i++)
  {
    dtsign[curevennode] = -dtsign[abscode[curevennode]];
    curevennode += 2;
  }

  /* decide arbitrarily the region sign of a node */

  regionsign[1] = 1;
  regionsign[abscode[1]] = -regionsign[1];

  reconstruct_sign ();
  numregions = 2 + numlabels - numnodes;
  // printf ("numregions: %d\n", numregions);
  printf ("#\n# tubular knot with Dowker-Thistletwait notation\n");
  printf ("# [");
  for (i = 1; i <= numlabels; i += 2)
  {
    if (i > 1) printf (", ");
    printf ("%d", dtsign[i]*abscode[i]);
  }
  printf ("]\n#\nsketch {\n");
  display_arcs_from_arcs ();
  display_arcs_from_nodes ();
  display_regions ();
  display_regions_from_arcs ();
  display_regions_from_nodes ();
  printf ("}\n");
}

/*
 * display_arcs_from_arcs
 * node arcs are numbered from 1 to numlabels, each produces two
 * apparent contour arcs, one on the right, same orientation, i -> 2*i - 1
 * one on the left, with opposite orientation, i -> 2*i
 */

void
display_arcs_from_arcs (void)
{
  int i;

  for (i = 1; i <= numlabels; i++)
  {
    printf ("Arc %d: [0];\n", 2*i - 1);
    printf ("Arc %d: [0];\n", 2*i);
  }
}

/*
 * display_arcs_from_nodes
 */

void
display_arcs_from_nodes (void)
{
  int i, overpass;
  int offset = 2*numlabels;

  printf ("# start of node arcs\n");
  for (i = 1; i <= numnodes; i++)
  {
    overpass = 0;
    if (dtsign[2*i-1] < 0) overpass = 2;
    printf ("Arc %d: [%d];\n", offset + 4*i - 3, 2 - overpass);
    printf ("Arc %d: [%d];\n", offset + 4*i - 2, overpass);
    printf ("Arc %d: [%d];\n", offset + 4*i - 1, 2 - overpass);
    printf ("Arc %d: [%d];\n", offset + 4*i - 0, overpass);
  }
}

/*
 * display_regions
 */

void
display_regions (void)
{
  int *redtagged, *bluetagged;
  int redregionnum = 0;
  int blueregionnum = 0;
  int i, arc;

  redtagged = (int *) malloc ( (numlabels+1) * sizeof(int) );
  bluetagged = (int *) malloc ( (numlabels+1) * sizeof(int) );
  for (i = 1; i <= numlabels; i++) redtagged[i] = bluetagged[i] = 0;

  //printf ("Displaying red regions...\n");

  for (i = 1; i <= numlabels; i++)
  {
    if (redtagged[i]) continue;
    redregionnum++;  /* found new region */
    //printf ("New red region: %d:\n", redregionnum);
    printf ("Region %d: ", redregionnum - 1);
    if (redregionnum == 1) printf ("() ");
    printf ("(");
    for (arc = i;;)
    {
      redtagged[arc] = redregionnum;
      if ( (arc % 2) == 1 )
      { /* odd arc */
        if (regionsign[arc] > 0)
        {
          arc = nextlabel(abscode[arc]);
          printf ("-a%d ", 2*arc);
          //printf ("+a%d ", arc);
        } else {
          arc = abscode[arc];
          printf ("-a%d ", 2*arc - 1);
          //printf ("-a%d ", arc);
        }
      } else {  /* even arc */
        if (regionsign[prevlabel(arc)] > 0)
        {
          arc = abscode[prevlabel(arc)];
          printf ("-a%d ", 2*arc - 1);
          //printf ("-a%d ", arc);
        } else {
          arc = nextlabel(abscode[prevlabel(arc)]);
          printf ("-a%d ", 2*arc);
          //printf ("+a%d ", arc);
        }
      }
      if (arc == i)
      {
        printf (");\n");
        break;
      }
    }
  }

  //printf ("Displaying blue regions (counterclockwise)...\n");

  for (i = 1; i <= numlabels; i++)
  {
    if (bluetagged[i]) continue;
    blueregionnum++;  /* found new region */
    //printf ("New blue region: %d:\n", blueregionnum);
    printf ("Region %d: ", redregionnum + blueregionnum - 1);
    printf ("(");
    for (arc = i;;)
    {
      bluetagged[arc] = blueregionnum;
      if ( (arc % 2) == 1 )
      { /* odd arc */
        if (regionsign[prevlabel(arc)] > 0)
        {
          arc = abscode[prevlabel(arc)];
          printf ("-a%d ", 2*arc - 1);
          //printf ("-a%d ", arc);
        } else {
          arc = nextlabel(abscode[prevlabel(arc)]);
          printf ("-a%d ", 2*arc);
          //printf ("+a%d ", arc);
        }
      } else {  /* even arc */
        if (regionsign[arc] > 0)
        {
          arc = nextlabel(abscode[arc]);
          printf ("-a%d ", 2*arc);
          //printf ("+a%d ", arc);
        } else {
          arc = abscode[arc];
          printf ("-a%d ", 2*arc - 1);
          //printf ("-a%d ", arc);
        }
      }
      if (arc == i)
      {
        printf (");\n");
        break;
      }
    }
  }

  assert (redregionnum + blueregionnum == numregions);
}

/*
 * display_regions_from_arcs
 * nodes are indexed by using the odd label i: (i-1)/2
 * then we have 4 arcs for each node, so there is a multiplying factor 4
 * for a final 2*(i-1) + h
 * h = 0,1,2,3
 */

void
display_regions_from_arcs ()
{
  int i;
  int offset = numregions - 1;
  int arcsoffset = 2*numlabels + 1;
  int arcahead, arcbehind;

  printf ("# elongated regions corresponding to arcs\n");
  for (i = 1; i <= numlabels; i++)
  {
    if ( (i % 2) == 1)
    {
      arcahead = 2*(i-1);
      arcbehind = 2*(abscode[prevlabel(i)] - 1) + 3;
    } else {
      arcahead = 2*(abscode[i]-1) + 1;
      arcbehind = 2*(prevlabel(i) - 1) + 2;
    }
    printf ("Region %d: ", offset + i);
    printf ("(+a%d -a%d +a%d -a%d);\n", 2*i - 1, arcsoffset + arcahead, 2*i, arcsoffset + arcbehind);
  }
}

/*
 * display_regions_from_nodes
 */

void
display_regions_from_nodes ()
{
  int i, oddarc, ori;
  int offset = numregions + numlabels;
  int arcsoffset = 2*numlabels - 1;

  printf ("# small nodal regions\n");
  for (i = 0; i < numnodes; i++)
  {
    oddarc = 2*i + 1;
    ori = regionsign[oddarc];
    oddarc = 2*oddarc + arcsoffset;
    printf ("Region %d: (+a%d +a%d +a%d +a%d);\n", offset + i, oddarc,
                                                               oddarc + 2 - ori,
                                                               oddarc + 2,
                                                               oddarc + 2 + ori);
  }
}

/*
 * reconstruct regionsign of nodes
 */

void
reconstruct_sign (void)
{
  int *tagged;
  int i, j, node, expansions, totexpansions = 0;
  int insist;

  tagged = (int *) malloc ( (2*numnodes+1) * sizeof(int) );

  insist = 1;
  while (insist)
  {
    insist = 0;
    for (i = 1; i <= numlabels; i++)
    {
      for (j = 1; j <= numlabels; j++) tagged[j] = 0;
      /* follow cicle */
      node = i;
      while (1)
      {
        if (tagged[node])
        {
          // printf ("found cycle starting at node: %d\n", abscode[node]);
          expansions = inherit (abscode[node]);
          if (expansions) insist = 1;
          totexpansions += expansions;
          // printf ("Expanded by %d nodes\n", expansions);
          break;
        }
        assert (tagged[abscode[node]] == 0);
        tagged[abscode[node]] = 1;
        node = nextlabel (node);
      }
    }
  }
  free (tagged);
  assert (totexpansions == numnodes - 1);
  // printf ("totexpansions: %d\n", totexpansions);

  //for (i = 1; i <= numlabels; i++)
  //{
  //  printf ("regionsign[%d] = %d\n", i, regionsign[i]);
  //}
}

int
inherit (int startlabel)
{
  int isinside;
  int endlabel, label;
  int expansions = 0;

  // printf ("searching inheritance for cycle starting at %d\n", startlabel);

  /* there is a cicle starting ad startlabel and ending
   * at abscode[startlabel]
   */

  endlabel = abscode[startlabel];
  assert (nextlabel(startlabel) != abscode[startlabel]);  /* no tight loops allowed */
  assert (nextlabel(endlabel) != abscode[endlabel]);  /* no tight loops allowed */

  /* now find where the path continuing from abscode[startlabel]
   * first intersects the loop
   */

  isinside = 0;
  for (label = nextlabel(endlabel); label != startlabel; label = nextlabel(label))
  {
    if (isinpath (abscode[label], startlabel, endlabel))
    {
      if (isinside)
      {
        // printf ("Sign Inheritance between %d and %d (EXITING)\n", startlabel, label);
        expansions += propagate (-regionsign[startlabel], label);
        expansions += propagate (-regionsign[label], startlabel);
      } else {
        // printf ("Sign Inheritance between %d and %d (ENTERING)\n", startlabel, label);
        expansions += propagate (regionsign[startlabel], label);
        expansions += propagate (regionsign[label], startlabel);
      }
      isinside = 1-isinside;
    }
  }
  return (expansions);
}

int
nextlabel (int label)
{
  label++;
  if (label > numlabels) label -= numlabels;
  return (label);
}

int
prevlabel (int label)
{
  label--;
  if (label < 1) label += numlabels;
  return (label);
}

int
isinpath (int i, int start, int end)
{
  if (end >= start) return (i >= start && i <= end);

  if (i >= start) return (1);
  if (i <= end) return (1);
  return (0);
}

int
propagate (int sign, int label)
{
  if (sign == 0) return (0);
  if (regionsign[label] != 0) return (0);

  regionsign[label] = sign;
  regionsign[abscode[label]] = -sign;
  return (1);
}
