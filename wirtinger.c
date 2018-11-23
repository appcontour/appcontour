#include <assert.h>
#include "contour.h"
#include "readdtcode.h"
#include "fundamental.h"
#include "wirtinger.h"

struct presentation *
wirtingerfromloiv (struct vecofintlist *loiv)
{
  int i, selflink, numnodes, numlabels;
  int over, label, startlabel, signature;
  int node;
  extern int debug;
  int *dtcode, *dt_involution, *dt_realization, *dt_over, *gregionsign;
  int *overgen, *ingen;
  int underpasslabel, overpasslabel, nextnode, nextlabel, nextunder;
  struct presentation *p;
  struct presentationrule *rule;
  int longlength;

  if (debug) printloiv (loiv);
  assert (loiv->next == 0);
  assert (loiv->handedness);
  dtcode = loiv->vec;
  numnodes = loiv->len;
  numlabels = 2*numnodes;
  gregionsign = loiv->handedness;

  dt_involution = (int *) malloc (numlabels * sizeof (int));
  dt_realization = (int *) malloc (numlabels * sizeof (int));
  dt_over = (int *) malloc (numlabels * sizeof (int));
  overgen = (int *) malloc (numnodes * sizeof (int));
  ingen = (int *) malloc (numnodes * sizeof (int));

  p = (struct presentation *) malloc (sizeof (struct presentation));
  p->gennum = numnodes;
  p->elements = 0;
  p->rules = 0;
  p->characteristic = 0;
  p->espected_deficiency = 1;

  selflink = 0;
  for (i = 0; i < numnodes; i++)
  {
    label = abs(dtcode[i]) - 1;
    dt_involution[2*i] = label;
    dt_involution[label] = 2*i;
    dt_over[2*i] = (dtcode[i]>0)?1:0;
    dt_over[label] = 1 - dt_over[2*i];
    dt_realization[2*i] = gregionsign[i];
    dt_realization[label] = -gregionsign[i];
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
    if (debug) printf ("Node: %d involving labels %d %d\n", i + 1, 2*i + 1, dt_involution[2*i] + 1);
    over = 1;
    if (dtcode[i] < 0) over = -1;
    if (debug) {
      printf (" overpass label: %d\n", overpasslabel + 1);
      printf (" overpass generator: %d\n", overgen[i] + 1);
      printf (" incoming generator: %d\n", ingen[i] + 1);
      printf (" outgoing generator: %d\n", i+1);
    }

    rule = (struct presentationrule *) malloc (4*sizeof (int) + sizeof (struct presentationrule));
    rule->length = 4;
    rule->next = p->rules;
    p->rules = rule;
    rule->var[0] = -(i+1);
    rule->var[2] = ingen[i] + 1;
    if (dt_realization[underpasslabel] > 0)
    {
      if (debug) {
        printf (" relator: %c = ", i + 'a');
        printf ("%c", overgen[i] + 'A');
        printf ("%c", ingen[i] + 'a');
        printf ("%c", overgen[i] + 'a');
        printf ("\n");
      }
      rule->var[1] = -(overgen[i] + 1);
      rule->var[3] = overgen[i] + 1;
    } else {
      if (debug) {
        printf (" relator: %c = ", i + 'a');
        printf ("%c", overgen[i] + 'a');
        printf ("%c", ingen[i] + 'a');
        printf ("%c", overgen[i] + 'A');
        printf ("\n");
      }
      rule->var[1] = overgen[i] + 1;
      rule->var[3] = -(overgen[i] + 1);
    }
    signature = over*gregionsign[i];
    assert (signature);
    selflink += signature;
    if (debug) printf (" selflinking: %c\n", (signature > 0)?'+':'-');
  }
  if (debug) printf ("Selflinking number: %d\n", selflink);
  /* the meridian is simply the first generator */
  rule = (struct presentationrule *) malloc (sizeof(int) + sizeof (struct presentationrule));
  rule->length = 1;
  rule->var[0] = 1;
  rule->next = 0;
  p->elements = rule;
  /* computing the longitude starting from the underpass of node 1 */
  longlength = numnodes + abs(selflink);
  rule = (struct presentationrule *) malloc (longlength*sizeof(int) + sizeof(struct presentationrule));
  rule->length = longlength;
  startlabel = 0;
  if (dtcode[0] > 0) startlabel = dt_involution[startlabel];
  if (debug) printf ("longitude starting at label %d\n", startlabel + 1);
  label = startlabel;
  i = 0;
  do {
    label++;
    if (label >= numlabels) label -= numlabels;
    if ( (label % 2) == 0)
    {
      node = label/2;
    } else {
      node = dt_involution[label]/2;
    }
    if (dt_over[label])
    {
      if (debug) printf ("Passing over in node %d\n", node+1);
      continue;
    }
    rule->var[i] = overgen[node] + 1;
    if (dt_realization[label] < 0) rule->var[i] = - overgen[node] - 1;
    i++;
    if (debug) printf ("Passing under generator %d in node %d\n", overgen[node]+1,node+1);
  } while (label != startlabel);

  while (i < longlength)
  {
    rule->var[i] = 1;
    if (selflink < 0) rule->var[i] = -1;
    i++;
  }

  p->elements->next = rule;

  free (dt_involution);
  free (dt_realization);
  free (dt_over);
  free (overgen);
  free (ingen);
  return (p);
}
