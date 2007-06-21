/*
 */

#include <assert.h>
#include "contour.h"
#include "hacon.h"

struct hacongraph *
compute_hacon (struct sketch *s)
{
  struct hacongraph *hg;
  struct region *r;
  struct arc *a;
  struct haconnode **hn;
  struct haconarc **ha;
  struct haconnode *hnpt;
  int arcf, st;
  int numhaconarcs;

  numhaconarcs = count_link_components (s);

  hg = (struct hacongraph *) malloc (sizeof (struct hacongraph));
  hg->node = 0;
  hg->arc = 0;
  hg->next = 0;
  hn = (struct haconnode **) malloc (s->regioncount * sizeof (struct haconnode *));
  hg->nodealloc = hn;
  ha = (struct haconarc **) malloc (s->arccount * sizeof (struct haconarc *));
  hg->nodealloc = hn;

  // if (tag_connected_components (s) < 0) return (-1);
  for (r = s->regions; r; r = r->next)
  {
    printf ("region: %d: %d strati\n", r->tag, r->f);
    assert (r->tag < s->regioncount);
    hn[r->tag] = 0;
    if (r->f > 0)
    {
      hnpt = hn[r->tag] = (struct haconnode *) malloc (r->f * sizeof (struct haconnode));
      for (st = 0; st < r->f; st++, hnpt++)
      {
        hnpt->strato = st;
        hnpt->region = r;
      }
    }
  }
  for (a = s->arcs; a; a = a->next)
  {
    arcf = a->regionleft->border->region->f + 
           a->regionright->border->region->f;
    arcf /= 2;
    printf ("arc: %d, %d strati\n", a->tag, arcf);
    assert (a->tag <= s->arccount);
    ha[a->tag] = (struct haconarc *) malloc (arcf * sizeof (struct haconarc));
  }
  return (hg);
}

void
print_hacon (struct hacongraph *h)
{
  printf ("not implemented\n");
}
