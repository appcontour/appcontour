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
static int *tagged;

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
void first_completion (int stackdim);
int next_completion (void);
int isconsistent (void);
int tour_of_region (int label, int velocity);
void walk_left (int *labelpt, int *velocitypt);

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
  tagged = (int *) malloc ( (2*numnodes+1) * sizeof(int) );
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

  numregions = 2 + numlabels - numnodes;
  reconstruct_sign ();
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
  free (tagged);
  free (regionsign);
  free (dtsign);
  free (abscode);
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
  int i, j, node, expansions, totexpansions = 0;
  int insist;

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
  if (totexpansions < numnodes - 1)
  {
    //printf ("totexpansions: %d instead of %d\n", totexpansions, numnodes - 1);
    first_completion (numnodes - 1 - totexpansions);
    while (! isconsistent () )
    {
      //printf ("completion is not consistent...\n");
      if (! next_completion () )
      {
        printf ("Cannot find any consistent knot diagram!\n");
        exit (11);
      }
    }
  }
  //assert (totexpansions == numnodes - 1);
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

/*
 *
 */

static int *stack;
static int *stacknode;
static int stackpt;
static int stackdim;

void
first_completion (int stdim)
{
  int i, label;

  stackdim = stdim;
  stack = (int *) malloc (stackdim * sizeof (int));
  stacknode = (int *) malloc (stackdim * sizeof (int));
  for (i = 0; i < numnodes; i++)
  {
    label = 2*i+1;
    if (regionsign[label]) continue;
    stack[stackpt] = 1;
    stacknode[stackpt++] = i;
    assert (stackpt <= stackdim);
    regionsign[label] = 1;
    regionsign[abscode[label]] = -1;
  }
  assert (stackpt == stackdim);
}

/*
 *
 */

int
next_completion (void)
{
  int node, label;

  while (stackpt > 0 && stack[--stackpt] < 0);
  if (stack[stackpt] < 0)
  {
    /* no more configurations possible */
    free (stack);
    free (stacknode);
    return (0);
  }
  while (stackpt < stackdim)
  {
    stack[stackpt] = -stack[stackpt];
    node = stacknode[stackpt++];
    label = 2*node + 1;
    regionsign[label] *= -1;  
    regionsign[abscode[label]] *= -1;  
  }
  return (1);
}

/*
 *
 */

int
isconsistent (void)
{
  int i;
  int countregions = 0;

  for (i = 1; i <= numlabels; i++) tagged[i] = 0;

  /* each arc should be walked once positively and once negatively
   * tag 0
   * tag 1 +
   * tag 2 -
   * tag 3 +-
   */

  while (1)
  {
    for (i = 1; i <= numlabels; i++)
    {
      if (tagged[i] < 3) break;
    }
    if (i > numlabels) break; /* complete tour! everything fine */

    countregions++;
    /* i is an arc that has not been walked twice */
    if (tagged[i] == 0 || tagged[i] == 2)
    {
      if (!tour_of_region (i, 1)) return (0);  /* something went wrong during tour */
    } else {
      if (!tour_of_region (i, -1)) return (0);  /* something went wrong during tour */
    }
  }
  if (countregions != numregions) return (0);
  //printf ("correct number of regions!\n");
  /* we NEED to find a test that works for sure! This is not enough */
  return (1);
}

/*
 * return 0 if something goes wrong
 */

int
tour_of_region (int startlabel, int startvelocity)
{
  int label = startlabel;
  int velocity = startvelocity;

  while (1)
  {
    if (velocity > 0)
    {
      if ((tagged[label] & 1)) return (0); /* should not be already visited */
      tagged[label] += 1;
    } else {
      if ((tagged[label] & 2)) return (0);
      tagged[label] += 2;
    }
    walk_left (&label, &velocity);
    if (label == startlabel && velocity == startvelocity) return (1);
  }
  return (0);
}

/*
 * walk around a region conterclockwise
 */

void
walk_left (int *labelpt, int *velocitypt)
{
  if (*velocitypt > 0)
  {
    if (regionsign[*labelpt] > 0)
    {
      /* velocity remains positive */
      *labelpt = nextlabel(abscode[*labelpt]);
    } else {
      *velocitypt = - *velocitypt;
      *labelpt = abscode[*labelpt];
    }
  } else {
    if (regionsign[prevlabel(*labelpt)] > 0)
    {
      /* velocity remains negative */
      *labelpt = abscode[prevlabel(*labelpt)];
    } else {
      *velocitypt = - *velocitypt;
      *labelpt = nextlabel(abscode[prevlabel(*labelpt)]);
    }
  }
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
