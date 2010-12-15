#include "morse2morse.h"

void print_patchendlist (struct patchend *inf);

void
collect_simple_arcs (struct mdesc *contour)
{
  struct line *l;
  struct event *ev;
  struct patchend *infinity, *patchendpt, *endleft, *endright, *endbefore, *endafter;
  struct patch *patch, *patchleft, *patchright, *patchbefore, *patchafter;
  int i, *newd, cusps, cuspsbefore, cuspsafter;

  infinity = (struct patchend *) malloc (sizeof (struct patchend));
  infinity->right = infinity->left = infinity;
  infinity->patch = 0;

  for (l = contour->lastline; l; l = l->prev)
  {
    //printf ("LINE: %d\n", l->tag);

    patchendpt = infinity->right;

    for (ev = l->events; ev; ev = ev->next)
    {
      //printf ("EV: %d: %d\n", ev->tag, ev->type);
      /*
       * operazione preliminare: assegnare il valore d
       * nel caso non fosse stato assegnato dopo un evento
       * CROSS (dvalues *deve* essere 1)
       */
      switch (ev->type)
      {
        case EVENT_VERT:
        case EVENT_CUSP:
          assert (patchendpt != infinity);
          patch = patchendpt->patch;
          if (patch->dvalues == 1 && patch->d[0] == INVALIDINT)
          {
            patch->d[0] = ev->huffman;
            if (ev->type == EVENT_CUSP && ev->orientation == ORIENT_DOWN)
              patch->d[0] -= ev->cuspsign;
          }
          break;

        case EVENT_MAX:
          assert (patchendpt != infinity);
          endleft = patchendpt;
          endright = endleft->right;
          assert (endright != infinity);
          if (endleft->patch->dvalues == 1 && endleft->patch->d[0] == INVALIDINT)
            endleft->patch->d[0] = ev->huffman;
          if (endright->patch->dvalues == 1 && endright->patch->d[0] == INVALIDINT)
            endright->patch->d[0] = ev->huffman;
          break;

        case EVENT_CROSS:
          assert (patchendpt != infinity);
          endleft = patchendpt;
          endright = endleft->right;
          assert (endright != infinity);
          if (endleft->patch->dvalues == 1 && endleft->patch->d[0] == INVALIDINT)
            endleft->patch->d[0] = ev->huffman;
          if (endright->patch->dvalues == 1 && endright->patch->d[0] == INVALIDINT)
            endright->patch->d[0] = ev->huffman2;
          break;
      }
      switch (ev->type)
      {
        case EVENT_VERT:
          patchendpt = patchendpt->right;
          break;

        case EVENT_CUSP:
          assert (patchendpt != infinity);
          assert (ev->orientation);
          patch = patchendpt->patch;
          /* allocating a new vector of dvalues */
          newd = (int *) malloc ((patch->dvalues + 1)*sizeof(int));
          if (ev->orientation == ORIENT_UP)
          {
            assert (patchendpt->isstart == 0);
            /* must append information at the end */
            for (i = 0; i < patch->dvalues; i++) newd[i] = patch->d[i];
            assert (patch->d[patch->dvalues-1] == ev->huffman);
            newd[patch->dvalues] = newd[patch->dvalues-1] + ev->cuspsign;
          } else {
            assert (patchendpt->isstart);
            /* must add information at the beginning */
            for (i = 0; i < patch->dvalues; i++) newd[i+1] = patch->d[i];
            newd[0] = newd[1] - ev->cuspsign;
            assert (newd[0] == ev->huffman);
          }
          free (patch->d);
          patch->d = newd;
          patch->dvalues++;
          patchendpt = patchendpt->right;
          break;

        case EVENT_MIN:
          /*
           * a new patch must be created *before* patchendpt.
           */
          patch = (struct patch *) malloc (sizeof(struct patch));
          patch->dvalues = 1;
          patch->d = (int *) malloc (sizeof (int));
          patch->d[0] = ev->huffman;
          /*
           * two new patchends...
           */
          endleft = (struct patchend *) malloc (sizeof(struct patchend));
          endright = (struct patchend *) malloc (sizeof(struct patchend));

          endleft->patch = endright->patch = patch;
          endleft->right = endright;
          endright->left = endleft;
          endright->right = patchendpt;
          endleft->left = patchendpt->left;
          patchendpt->left->right = endleft;
          patchendpt->left = endright;
          if (ev->orientation == ORIENT_RIGHT)
          {
            endleft->isstart = 1;
            endright->isstart = 0;
            patch->start = endleft;
            patch->end = endright;
          } else {
            endleft->isstart = 0;
            endright->isstart = 1;
            patch->start = endright;
            patch->end = endleft;
          }
//print_patchendlist (infinity);
          break;

        case EVENT_MAX:
          /*
           * merging of two patches, or closing of a single patch
           */
//printf("MAX(before)\n");print_patchendlist (infinity);
          endleft = patchendpt;
          endright = endleft->right;
          assert (endright != infinity);
          patchendpt = endright->right;
          endleft->left->right = endright->right;
          patchendpt->left = endleft->left;
          patchleft = endleft->patch;
          patchright = endright->patch;
          if (ev->orientation == ORIENT_LEFT)
          {
            patchbefore = patchright;
            patchafter = patchleft;
            endbefore = endright;
            endafter = endleft;
          } else {
            patchbefore = patchleft;
            patchafter = patchright;
            endbefore = endleft;
            endafter = endright;
          }
          assert (endafter->isstart);
          assert (endbefore->isstart == 0);
          assert (patchafter->start == endafter);
          assert (patchbefore->end == endbefore);
          if (patchleft == patchright)
          {
            /*
             * closing a patch
             */
            cusps = patchleft->dvalues - 1;
            assert (patchleft->d[0] == patchleft->d[cusps]);
            /*
             * attach whole information at this event
             */
            ev->cusps = patchleft;
            patchleft->start = patchleft->end = 0;
          } else {
            /*
             * two different patches terminate/start here
             */
            cuspsbefore = patchbefore->dvalues - 1;
            cuspsafter = patchafter->dvalues - 1;
            assert (patchbefore->d[cuspsbefore] == patchafter->d[0]);
            newd = (int *) malloc ((cuspsbefore + cuspsafter + 1)*sizeof(int));
            for (i = 0; i < cuspsbefore; i++) newd[i] = patchbefore->d[i];
            for (i = 0; i <= cuspsafter; i++)
              newd[i + cuspsbefore] = patchafter->d[i];
            free (patchbefore->d);
            free (patchafter->d);
            patchbefore->d = newd;
            patchbefore->dvalues = cuspsbefore + cuspsafter + 1;
            patchbefore->end = patchafter->end;
            free (patchafter);
            if (patchbefore->end) patchbefore->end->patch = patchbefore;
            if (patchbefore->start == 0 && patchbefore->end == 0)
            {
              ev->cusps = patchbefore;
            }
          }
          free (endleft);
          free (endright);
//printf("MAX(after)\n");print_patchendlist (infinity);
          break;

        case EVENT_CROSS:
          endleft = patchendpt;
          endright = endleft->right;
          patchendpt = endright->right;
          patchleft = endleft->patch;
          patchright = endright->patch;
          if (ev->orientation == ORIENT_DOWN)
          {
            assert (patchleft->start == endleft);
            patchleft->start = 0;
            if (patchleft->end == 0) ev->cusps = patchleft;
          } else {
            assert (patchleft->end == endleft);
            patchleft->end = 0;
            if (patchleft->start == 0) ev->cusps = patchleft;
          }
          if (ev->orientation2 == ORIENT_DOWN)
          {
            assert (patchright->start == endright);
            patchright->start = 0;
            if (patchright->end == 0) ev->cusps2 = patchright;
          } else {
            assert (patchright->end == endright);
            patchright->end = 0;
            if (patchright->start == 0) ev->cusps2 = patchright;
          }
          endleft->patch = (struct patch *) malloc (sizeof (struct patch));
          endleft->patch->d = (int *) malloc (sizeof(int));
          endleft->patch->dvalues = 1;
          endleft->patch->d[0] = INVALIDINT;
          endright->patch = (struct patch *) malloc (sizeof (struct patch));
          endright->patch->d = (int *) malloc (sizeof(int));
          endright->patch->dvalues = 1;
          endright->patch->d[0] = INVALIDINT;
          if (ev->orientation2 == ORIENT_DOWN)
          {
            endleft->patch->end = 0;
            endleft->patch->start = endleft;
            endleft->isstart = 1;
          } else {
            endleft->patch->start = 0;
            endleft->patch->end = endleft;
            endleft->isstart = 0;
          }
          if (ev->orientation == ORIENT_DOWN)
          {
            endright->patch->end = 0;
            endright->patch->start = endright;
            endright->isstart = 1;
          } else {
            endright->patch->start = 0;
            endright->patch->end = endright;
            endright->isstart = 0;
          }
          break;

        default:
          printf ("caso non gestito: %d\n", ev->type);
          exit (1);
      }
    }
    assert (patchendpt == infinity);
//print_patchendlist(infinity);
  }
}

void
print_patchendlist (struct patchend *inf)
{
  struct patchend *pept;
  struct patch *patch;

  pept = inf->right;
  assert (pept->left == inf);
  printf ("start of list\n");
  while (pept != inf)
  {
    patch = pept->patch;
    printf ("patch-end %p: isstart=%d patch: %p%c%d%c\n", 
      pept, pept->isstart, patch, 
      patch->start?'(':'[',
      patch->dvalues,
      patch->end?')':']');
    pept = pept->right;
    assert (pept == pept->right->left);
    assert (pept == pept->left->right);
  }
  printf ("end of list\n");
}
