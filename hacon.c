/*
 */

#include <assert.h>
#include "contour.h"
#include "hacon.h"

struct hacongraph *
compute_hacon (struct sketch *s)
{
  struct hacongraph *hg;
  extern int verbose;
  int numhaconarcs, numhaconnodes;
  int **data;

  data = init_hacon_strata (s);
  numhaconnodes = tag_hacon_strata (data, s);
  printf ("Number of hacon graph nodes: %d\n", numhaconnodes);
  if (verbose) describe_hacon_nodes (numhaconnodes, data, s);
  numhaconarcs = count_link_components (s);
  printf ("Number of hacon graph arcs: %d\n", numhaconarcs);

  hg = (struct hacongraph *) malloc (sizeof (struct hacongraph));
  return (hg);
}

/*
 * we want to tag each stratum according to
 * the hacon node to which they belong
 * we do this by first tagging all with a negative
 * dummy integer, then 0 is the tag of the first
 * hacon node.
 * in a cicle we look for an untagged stratum and give it
 * a *new* hacon tag, then we loop on all strata and try to
 * locally extend the hacon node by local adjacency
 */

int
tag_hacon_strata (int **data, struct sketch *s)
{
  int k, hacon_tag;
  int *rstrata;
  struct region *r;

  /* reset all tags to -1 */
  for (r = s->regions; r; r = r->next)
  {
    for (k = 0; k < r->f; k++)
    {
      rstrata = data[r->tag];
      rstrata[k] = -1;
    }
  }

  hacon_tag = 0;
  while (single_tag_hacon_strata(hacon_tag, data, s)) hacon_tag++;

  return (hacon_tag);
}

int
single_tag_hacon_strata (int tag, int **data, struct sketch *s)
{
  struct region *r;
  int i, k, c, count, found;
  int *rstrata;

  found = 0;
  for (r = s->regions; r; r = r->next)
  {
    i = r->tag;   /* tags start from zero */
    rstrata = data[i];
    for (k = 0; k < r->f; k++)
    {
      if (rstrata[k] < 0)
      {
        found = 1;
        break;
      }
    }
    if (found) break;
  }

  if (found == 0) return (0);  /* non ho trovato strati non etichettati */
  /* trovata una regione (r) con una strato (k) non etichettata */
  rstrata = data[r->tag];
  rstrata[k] = tag;

  count = 1;
  while ((c = hacon_try_expand_node (tag, data, s)))
  {
    count += c;
  }
  return (count);
}

int
hacon_try_expand_node (int tag, int **data, 
  struct sketch *s)
{
  int count = 0;
  int k;
  struct region *r;
  struct border *bp, *bpstart;
  struct borderlist *bl;
  int *rstrata;

  for (r = s->regions; r; r = r->next)
  {
    rstrata = data[r->tag];
    for (k = 0; k < r->f; k++)
    {
      if (rstrata[k] != tag) continue;
      for (bl = r->border; bl; bl = bl->next)
      {
        bpstart = bl->sponda;
        if (bpstart == 0) continue;
        bp = bpstart;
        do {
          count += local_hacon_try_expand_node (data, bp, k);
        } while (bp = bp->next, bp != bpstart);
      }
    }
  }
  return (count);
}

int
local_hacon_try_expand_node (int **data, 
  struct border *bp, int k)
{
  int count, d, i, dmin, dmax, ori;
  int htag, rtag, stag;
  struct arc *a;
  struct border *btrans;
  int *rdata, *sdata, *r1data, *s1data;

  a = bp->info;
  btrans = gettransborder (bp);
  stag = btrans->border->region->tag;
  rtag = bp->border->region->tag;
  assert (stag != rtag);
  rdata = data[rtag];
  sdata = data[stag];
  htag = rdata[k];
  assert (htag >= 0);
  ori = bp->orientation;
  /* find dmin and dmax */
  dmin = dmax =  a->depths[0];
  for (i = 1; i < a->dvalues; i++)
  {
    if (a->depths[i] < dmin) dmin = a->depths[i];
    if (a->depths[i] > dmax) dmax = a->depths[i];
  }

  if (k < dmin)
  {
    if (sdata[k] < 0)
    {
      sdata[k] = htag;
      return (1);
    }
    assert (sdata[k] == htag);
    return (0);
  }

  if (ori > 0 && k >= dmax + 2)
  {
    if (sdata[k-2] < 0)
    {
      sdata[k-2] = htag;
      return (1);
    }
    assert (sdata[k-2] == htag);
    return (0);
  }

  if (ori < 0 && k >= dmax)
  {
    if (sdata[k+2] < 0)
    {
      sdata[k+2] = htag;
      return (1);
    }
    assert (sdata[k+2] == htag);
    return (0);
  }

  if ((dmin % 2) != (k % 2)) dmin++;
  r1data = rdata;
  s1data = sdata;
  if (ori < 0) {r1data = sdata; s1data = rdata;}
  count = 0;
  for (d = dmin; d < dmax; d += 2)
  {
    if (s1data[d] < 0)
    {
      s1data[d] = htag;
      count++;
    } else assert (s1data[d] == htag);
  }
  for (d = dmin; d < dmax + 2; d += 2)
  {
    if (r1data[d] < 0)
    {
      r1data[d] = htag;
      count++;
    } else assert (r1data[d] == htag);
  }

  return (count);
}

int **
init_hacon_strata (struct sketch *s)
{
  int rnum, i, rtag;
  struct region *r;
  int **data;

  rnum = s->regioncount;
  data = (int **) 
         malloc (rnum*sizeof(int *));

  for (i = 0; i < rnum; i++) data[i] = 0;
  for (r = s->regions; r; r = r->next)
  {
    rtag = r->tag;
    if (r->f > 0) 
      data[rtag] = (int *)
            malloc (r->f*sizeof(int));
  }
  return (data);
}

void
describe_hacon_nodes (int numhaconnodes, int **data, 
  struct sketch *s)
{
  int htag, k;
  struct region *r;
  int *rstrata;

  for (htag = 0; htag < numhaconnodes; htag++)
  {
    printf ("Hacon node %d:\n", htag);
    for (r = s->regions; r; r = r->next)
    {
      rstrata = data[r->tag];
      for (k = 0; k < r->f; k++)
      {
        if (rstrata[k] == htag)
        {
          printf ("  region %d stratum %d\n", r->tag, k);
        }
      }
    }
  }
}

/*
 * display hacon graph
 */

void
print_hacon (struct hacongraph *h)
{
  printf ("not implemented\n");
}
