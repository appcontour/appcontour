#include <assert.h>
#include <strings.h>
#include "contour.h"

extern int debug;

/*
 * definizione regole di trasformazione per superfici isotope
 */

int
apply_rule (char *rule, struct sketch *sketch)
{
  int res, rcount = 1;
  char *chpt;

  canonify (sketch);
  if (appcontourcheck (sketch, 0) == 0)
  {
    fprintf (stderr, "This sketch is NOT an apparent contour, ");
    fprintf (stderr, "cannot apply rules\n");
    exit (13);
  }

  for (chpt = rule; *chpt; chpt++)
  {
    if (*chpt == ':')
    {
      *chpt++ = 0;
      rcount = atoi (chpt);
    }
  }
  free_connected_components (sketch);
  res = 0;
  if (strcasecmp (rule, "n1") == 0) res = rule_n14 (sketch, 1, rcount);
  else if (strcasecmp (rule, "n2") == 0) res = rule_n14 (sketch, 2, rcount);
  else if (strcasecmp (rule, "n3") == 0) res = rule_n14 (sketch, 3, rcount);
  else if (strcasecmp (rule, "n4") == 0) res = rule_n14 (sketch, 4, rcount);
  else if (strcasecmp (rule, "n5") == 0) res = rule_n5 (sketch, rcount);
  else if (strcasecmp (rule, "n6") == 0) res = rule_n6 (sketch, rcount);
  else if (strcasecmp (rule, "c1") == 0) res = rule_c1 (sketch, rcount);
  else if (strcasecmp (rule, "c2") == 0) res = rule_c2 (sketch, rcount);
  else if (strcasecmp (rule, "a1") == 0) res = rule_a1 (sketch, rcount);
  else if (strcasecmp (rule, "a2") == 0) res = rule_a2 (sketch, rcount);
  else if (strcasecmp (rule, "cn1") == 0) res = rule_cn1 (sketch, rcount);
  else if (strcasecmp (rule, "cn2l") == 0) res = rule_cn2l (sketch, rcount);
  else if (strcasecmp (rule, "cn2r") == 0) res = rule_cn2r (sketch, rcount);
  else if (strcasecmp (rule, "cn2lb") == 0) res = rule_cn2lb (sketch, rcount);
  else if (strcasecmp (rule, "cn2rb") == 0) res = rule_cn2rb (sketch, rcount);
  else if (strcasecmp (rule, "cn3") == 0) res = rule_cn3 (sketch, rcount);
  else printf ("Invalid rule %s\n", rule);

  if (debug) printf ("res = %d\n", res);
  if (res)
  {
    checkconsistency (sketch);
    postprocesssketch (sketch);
    canonify (sketch);
  }
  return (res);
}

int
testallrules (struct sketch *sketch)
{
  int ru14, rcount, exitcode = 0;
  int res;

  printf ("Rules that apply:\n");
  for (ru14 = 1; ru14 <=4; ru14++)
  {
    rcount = 1;
    while (1)
    {
      if (debug) printf ("trying rule N%d:%d\n", ru14, rcount);
      res = rule_n14 (sketch, ru14, -rcount);
      if (res)
      {
        if (exitcode == 1) printf (" ");
        printf ("N%d", ru14);
        if (rcount > 1) printf (":%d", rcount);
        rcount++;
        exitcode = 1;
      } else break;
    }
  }
  exitcode = testsinglerule ("N5", rule_n5, exitcode, sketch);
  exitcode = testsinglerule ("N6", rule_n6, exitcode, sketch);
  exitcode = testsinglerule ("C1", rule_c1, exitcode, sketch);
  exitcode = testsinglerule ("C2", rule_c2, exitcode, sketch);
  exitcode = testsinglerule ("A1", rule_a1, exitcode, sketch);
  exitcode = testsinglerule ("A2", rule_a2, exitcode, sketch);
  exitcode = testsinglerule ("CN1", rule_cn1, exitcode, sketch);
  exitcode = testsinglerule ("CN2L", rule_cn2l, exitcode, sketch);
  exitcode = testsinglerule ("CN2R", rule_cn2r, exitcode, sketch);
  exitcode = testsinglerule ("CN2LB", rule_cn2lb, exitcode, sketch);
  exitcode = testsinglerule ("CN2RB", rule_cn2rb, exitcode, sketch);
  exitcode = testsinglerule ("CN3", rule_cn3, exitcode, sketch);
  printf ("\n");
  return (exitcode);
}

int
testsinglerule (char *rname, int (*rulefunc)(struct sketch *, int), 
                int exitcode, struct sketch *sketch)
{
  int rcount = 1, res;

  while (1)
  {
    if (debug) printf ("trying rule %s:%d\n", rname, rcount);
    res = (*rulefunc) (sketch, -rcount);
    if (res)
    {
      if (exitcode == 1) printf (" ");
      printf ("%s", rname);
      if (rcount > 1) printf (":%d", rcount);
      rcount++;
      exitcode = 1;
    } else break;
  }
  return (exitcode);
}

/*
 * rcount indica in quale situazione applicare
 * rcount = 1 -> applica appena trova una regione ok
 * rcount < 0 significa solo testare se la regola
 * e' applicabile, ma non applicarla
 */

int
rule_n14 (struct sketch *sketch, int rule, int rcount)
{
  struct region *r;
  struct border *bp, *bstart;
  struct border *outnw;
  int ok, count, gdia;
  int isrule23 = 0, req_d_increase = 0, ori = 0;
  int onlytest = 0;

  switch (rule)
  {
    case 1:
    ori = 1;
    break;

    case 2:
    isrule23 = 1;
    req_d_increase = 2;
    break;

    case 3:
    isrule23 = 1;
    break;

    case 4:
    ori = -1;
    break;
  }

  if (rcount < 0) {onlytest = 1; rcount *= -1;}
  if (debug && onlytest) printf ("only testing rule, rc %d\n", rcount);
  if (debug) printf ("rule23 %d, req_d %d, ori %d\n",
             isrule23, req_d_increase, ori);

  /* cerca una regione candidata, che deve avere
   * solo il bordo esterno, composto da due tratti
   * con orientazione positiva
   */

  for (r = sketch->regions; r; r = r->next)
  {
    if (r->border->next) continue;
    bstart = r->border->sponda;
    ok = 1;
    count = 0;
    if (isrule23 && bstart->orientation < 0) bstart = bstart->next;
    bp = bstart;
    do {
      count++;
      if (isrule23 == 0 && ori*bp->orientation < 0) {ok = 0; break;}
      if (count > 2) {ok = 0; break;}
      if (bp->info->cusps > 0) {ok = 0; break;}
      if (bp->info->endpoints != 2) {ok = 0; break;}
      bp = bp->next;
    } while (bp != bstart);
    if ((gdia = get_d_increase_across_node (bstart->info, 1)) !=
         get_d_increase_across_node (bstart->info, -1)) ok = 0;
    if (ok && isrule23)
    {
      if (bstart->orientation < 0 || bstart->next->orientation > 0) ok = 0;
      if (ok && 
        gdia != req_d_increase) ok = 0;
    }
    if (ok && count == 2 && rcount-- <= 1)
    {
      if (onlytest) return (1);
      fprintf (stderr, "Trovata regione candidata: %d\n", r->tag);
      outnw = rimuovi_losanga (r, sketch);
      if (debug) printsketch (sketch);
      taglia_nodo (outnw, sketch);
      postprocesssketch (sketch);
      //canonify (sketch);
      return (1);
    }
  }
  return (0);
}

int
rule_n5 (struct sketch *sketch, int rcount)
{
  struct region *r;
  struct border *bp, *bstart;
  int ok, count, totalin;
  int onlytest = 0;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata, che deve avere
   * esattamente tre lati (solo bordo esterno)
   * di cui uno passa sotto gli altri due
   */

  for (r = sketch->regions; r; r = r->next)
  {
    if (r->border->next) continue;
    bstart = r->border->sponda;
    ok = 1;
    count = 0;
    totalin = 0;
    bp = bstart;
    do {
      count++;
      totalin += abs (get_d_increase_across_node (bp->info, bp->orientation));
      if (count > 3) {ok = 0; break;}
      if (bp->info->cusps > 0) {ok = 0; break;}
      if (bp->info->endpoints != 2) {ok = 0; break;}
      bp = bp->next;
    } while (bp != bstart);
    if (totalin == 0 || totalin == 6) ok = 0;
    if (ok && count == 3 && rcount-- <= 1)
    {
      if (onlytest) return (1);
      fprintf (stderr, "Trovata regione candidata: %d\n", r->tag);
      triple_switch (bstart);
      //canonify (sketch);
      return (1);
    }
  }
  return (0);
}

int
rule_n6 (struct sketch *sketch, int rcount)
{
  return (rule_cn1_n6 (sketch, rcount, 1));
}

int
rule_c1 (struct sketch *sketch, int rcount)
{
  struct region *r;
  struct border *bstart;
  int onlytest = 0;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata, che deve avere
   * per bordo un s1 (senza buchi) con esattamenre
   * due cuspidi
   */

  for (r = sketch->regions; r; r = r->next)
  {
    if (r->border->next) continue;
    bstart = r->border->sponda;
    if (bstart->next != bstart) continue;
    if (bstart->orientation < 0) continue;
    if (bstart->info->endpoints != 0) continue;
    if (bstart->info->cusps != 2) continue;

    if (rcount-- <= 1)
    {
      if (onlytest) return (1);
      fprintf (stderr, "Trovata regione candidata: %d\n", r->tag);
      remove_s1 (bstart, sketch);
      return (1);
    }
  }
  return (0);
}

int
rule_c2 (struct sketch *sketch, int rcount)
{
  int onlytest = 0;
  struct region *r;
  struct border *cusp1, *cusp2;
  int i, j, cusp1pos, cusp2pos, da, db;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata, che deve avere
   * due cuspidi rientranti (orientazione negativa)
   * con gli stessi due valori di d in ordine opposto
   */

  for (r = sketch->regions; r; r = r->next)
  {
    for (i = 0;; i++)             /* ciclo sulle cuspidi */
    {
      if ((cusp1 = get_ith_cusp (r, i, &cusp1pos)) == 0) break;
      if (cusp1->orientation > 0) continue;
      da = cusp1->info->depths[cusp1pos];
      db = cusp1->info->depths[cusp1pos+1];
      for (j = i+1;; j++)
      {
        if ((cusp2 = get_ith_cusp (r, j, &cusp2pos)) == 0) break;
        if (cusp2->orientation > 0) continue;
        if (cusp2->info->depths[cusp2pos] != db) continue;
        if (cusp2->info->depths[cusp2pos+1] != da) continue;
        if (debug) printf ("coppia di cuspidi: %d, %d\n", i, j);
        /* ho trovato una coppia di cuspidi abbinabile */
        if (rcount-- <= 1)
        {
          if (onlytest) return (1);
          fprintf (stderr, "Trovata regione: %d c1 = %d c2 = %d\n", 
                           r->tag, i, j);
          if (debug) printf ("prima di chiamare join_cusps\n");
          if (debug) printsketch (sketch);
          join_cusps (cusp1, cusp1pos, cusp2, cusp2pos, sketch);
          if (debug) checkconsistency (sketch);
          taglia_nodo (cusp1->next, sketch);
          return (1);
        }
      }
    }
  }
  return (0);
}

int
rule_a1 (struct sketch *sketch, int rcount)
{
  return (rule_a12 (sketch, rcount, 1));
}

int rule_a2 (struct sketch *sketch, int rcount)
{
  return (rule_a12 (sketch, rcount, -1));
}

int
rule_a12 (struct sketch *sketch, int rcount, int ddiff)
{
  struct region *r;
  struct borderlist *extbl, *inbl;
  struct border *extb, *inb;
  int onlytest = 0;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata, che deve essere
   * un anello con orientazione positiva, senza
   * cuspidi, d_in = d_out + 1
   */

  for (r = sketch->regions; r; r = r->next)
  {
    extbl = r->border;
    if (extbl->sponda == 0) continue;                 /* regione esterna */
    if ((inbl = extbl->next) == 0) continue;  /* non ci sono buchi */
    if (inbl->next != 0) continue;            /* ci sono piu' buchi */
    extb = extbl->sponda;
    inb = inbl->sponda;
    if (extb->info->endpoints != 0) continue;
    if (inb->info->endpoints != 0) continue;  /* non sono degli S1 */
    if (extb->orientation <= 0) continue;
    if (inb->orientation <= 0) continue;
    if (extb->info->cusps > 0) continue;
    if (inb->info->cusps > 0) continue;
    if (inb->info->depths[0] != extb->info->depths[0] + ddiff) continue;
    
    if (rcount-- <= 1)
    {
      if (onlytest) return (1);
      fprintf (stderr, "Trovata regione candidata: %d\n", r->tag);
      remove_annulus (r, sketch);
      return (1);
    }
  }
  return (0);
}

int
rule_cn1 (struct sketch *sketch, int rcount)
{
  return (rule_cn1_n6 (sketch, rcount, 0));
}

int
rule_cn1_n6 (struct sketch *sketch, int rcount, int isn6)
{
  struct region *r;
  struct borderlist *extbl;
  struct border *extb;
  struct arc *arc;
  int ori, diff, onlytest = 0;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata, che deve essere senza buchi
   * circondata da una curva unica con orientazione
   * positiva, un endpoint e due cuspidi
   * con valori opportuni di "d"
   * per la regola n6 non ci devono essere cuspidi!
   */

  ori = 1;
  if (isn6) ori = -1;
  for (r = sketch->regions; r; r = r->next)
  {
    extbl = r->border;
    if (extbl->sponda == 0) continue;         /* regione esterna */
    if (extbl->next != 0) continue;           /* non deve avere buchi */
    extb = extbl->sponda;
    if (ori*extb->orientation <= 0) continue;
    if (extb->next != extb) continue;
    arc = extb->info;
    if (arc->endpoints != 1) continue;
    if (isn6)
    {
      if (arc->cusps != 0) continue;
      diff = 0;
    } else {
      if (arc->cusps != 2) continue;        /* devono esserci due cuspidi */
      diff = arc->depths[2] - arc->depths[0];
      if (diff != 2 && diff != -2) continue;  /* d crescente o decrescente */
      if (get_d_increase_across_node (arc, -1) != 
          diff + get_d_increase_across_node (arc, 1)) continue;
                                            /* incosistent d values at node */
    }

    if (rcount-- <= 1)
    {
      if (onlytest) return (1);
      fprintf (stderr, "Trovata regione candidata: %d\n", r->tag);
      remove_ear (r, sketch);
      return (1);
    }
  }
  return (0);
}

int
rule_cn2l (struct sketch *sketch, int rcount)
{
  return (rule_cn2lr (sketch, rcount, 1, 0));
}

int
rule_cn2r (struct sketch *sketch, int rcount)
{
  return (rule_cn2lr (sketch, rcount, 0, 0));
}

int
rule_cn2lb (struct sketch *sketch, int rcount)
{
  return (rule_cn2lr (sketch, rcount, 1, 1));
}

int
rule_cn2rb (struct sketch *sketch, int rcount)
{
  return (rule_cn2lr (sketch, rcount, 0, 1));
}

int
rule_cn2lr (struct sketch *sketch, int rcount, int isleft, int isback)
{
  struct region *r;
  struct borderlist *bl;
  struct border *bpstart, *bp, *bpn;
  struct arc *arcin, *arcout;
  int ori, orib, dincr, onlytest = 0;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata: due archi consecutivi
   * (eventualmente coincidenti) orientati positivamente.
   * "d" diminuisce (aumenta se "isback") di due unita', 
   * il secondo arco deve avere una prima cuspide con d 
   * che aumenta (diminuisce) di uno
   */

  ori = orib = 1;
  if (! isleft) ori = -1;
  if (isback) orib = -1;
  for (r = sketch->regions; r; r = r->next)
  {
    for (bl = r->border; bl; bl = bl->next)
    {
      bpstart = bl->sponda;
      if (bpstart == 0) continue;
      if (bpstart->info->endpoints == 0) continue;
      bp = bpstart;
      do {
        if (bp->orientation < 0) continue;
        bpn = bp->next;
        if (bpn->orientation < 0) continue;
        arcin = bp->info;
        arcout = bpn->info;
        dincr = arcout->depths[0] - arcin->depths[arcin->dvalues - 1];
        if (!isleft)     /* scambio i ruoli */
        {
          arcin = bpn->info;
          arcout = bp->info;
        }
        if (dincr != - orib*ori*2) continue;
        if (arcout->cusps < 1) continue;
        if (isleft)
        {
          if (arcout->depths[1] - arcout->depths[0] != orib) continue;
        } else {
          if (arcout->depths[arcout->dvalues - 2] -
              arcout->depths[arcout->dvalues - 1] != orib) continue;
        }
        if (get_d_increase_across_node (arcin, ori) != -1-orib) continue;
        if (get_d_increase_across_node (arcout, -ori) != orib-1) continue;
        if (rcount-- <= 1)
        {
          if (onlytest) return (1);
          fprintf (stderr, "Trovata regione: %d, arcin %d, arcout %d\n", 
                           r->tag, arcin->tag, arcout->tag);
          // chiama la funzione che applica la regola
          applyrulecn2 (bpn, arcout, ori, orib, sketch);
          return (1);
        }
      } while (bp = bp->next, bp != bpstart);
    }
  }
  return (0);
}

int
rule_cn3 (struct sketch *sketch, int rcount)
{
  struct region *r;
  struct borderlist *extbl;
  struct border *b1, *b2, *btemp;
  struct arc *arc1, *arc2;
  int onlytest = 0;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata: solo bordo esterno, con due archi,
   * uno senza cuspidi e uno con una cuspide e orientato positivamente
   * valori uscenti di "d" dai due nodi devono essere compatibili
   */

  for (r = sketch->regions; r; r = r->next)
  {
    extbl = r->border;
    if (extbl->sponda == 0) continue;         /* regione esterna */
    if (extbl->next != 0) continue;           /* non deve avere buchi */
    b1 = extbl->sponda;
    b2 = b1->next;
    if (b2->next != b1) continue;
    if (b1->info->cusps != 0)              /* scambio le sponde */
    {
      btemp = b1;
      b1 = b2;          /* b1 senza cuspide, b2 una cuspide */
      b2 = btemp;
    }
    arc1 = b1->info;
    if (arc1->cusps != 0) continue;
    if (get_d_increase_across_node (arc1, 1) !=
        get_d_increase_across_node (arc1, -1) ) continue;
    if (b2->orientation < 0) continue;
    arc2 = b2->info;
    if (arc2->cusps != 1) continue;
    if (get_d_increase_across_node (arc2, -1) !=
        get_d_increase_across_node (arc2, 1) ) continue;

    if (rcount-- <= 1)
    {
      if (onlytest) return (1);
      fprintf (stderr, "Trovata regione candidata: %d\n", r->tag);
      // chiama la funzione che esegue la "rule"
      remove_cusp (r, sketch);
      return (1);
    }
  }
  return (0);
}

/*
 * auxiliary functions used by rules
 */

struct border *
get_ith_cusp (struct region *r, int i, int *cusppos)
{
  struct borderlist *bl;
  struct border *bp;
  int ii = 0, j;

  for (bl = r->border; bl; bl = bl->next)
  {
    if (bl->sponda == 0) continue;
    bp = bl->sponda;
    do {
      for (j = 0; j < bp->info->cusps; j++)
      {
        if (i == ii++)          /* trovata */
        {
          *cusppos = j;
          return (bp);
        }
      }
      bp = bp->next;
    } while (bp != bl->sponda);
  }
  return (0);
}

/*
 * trasformazione topologica del tipo
 * 
 *   \ /      \   /
 *    X        \ /
 *   / \  ==>   X
 *   \ /       / \
 *    X       /   \
 *   / \     /     \
 *
 */

struct border *
rimuovi_losanga (struct region *r, struct sketch *sketch)
{
  struct border *b1, *b2, *b1trans, *b2trans, *b1transn, *b2transn;
  struct arc *arc1, *arc2;

  b1 = r->border->sponda;
  b2 = b1->next;
  arc1 = b1->info;
  arc2 = b2->info;
  //assert (arc1->regionleft == b1);
  //assert (arc2->regionleft == b2);
  b1trans = gettransborder (b1);
  b2trans = gettransborder (b2);
  arc1->regionright = arc1->regionleft = 0;
  arc2->regionright = arc2->regionleft = 0;
  removearc (arc1, sketch);
  removearc (arc2, sketch);
  b1trans->info = b2trans->info = 0;
  freeborder (b1);
  r->border->sponda = 0;
  assert (r->border->next == 0);
  free (r->border);
  r->border = 0;
  removeregion (r, sketch);
  b1transn = b1trans->next;
  b2transn = b2trans->next;
  assert (b1transn != b1trans);
  assert (b2transn != b2trans);
  if (b1transn->next == b1trans) b1transn->info->endpoints = 1;
  if (b2transn->next == b2trans) b2transn->info->endpoints = 1;
  removeborder (b1trans);
  removeborder (b2trans);
  return (b1transn);
}

/*
 * trasformazione topologica del tipo
 * 
 *      \   /    \     /
 *    b1n\ /      \   /
 *        X  ===>  \ /
 *       / \       / \
 *      /   \     /   \
 *     /     \   /     \
 *
 * assumiamo una orientazione consistente
 * degli archi coinvolti
 */

void
taglia_nodo (struct border *b1n, struct sketch *sketch)
{
  struct border *b1p, *b2p, *b3p, *b4p, *b2n;
  struct region *r2, *r4;
  struct arc *arc12, *arc23, *arc34, *arc41, *arcleft, *arcright;
  int changes, r2ext, r4ext;

  if (debug) printf ("entering taglia_nodo\n");
  if (debug) printsketch (sketch);
  if (debug) printf ("b1n e' nella regione %d\n", b1n->border->region->tag);
  assert (b1n->info);
  
  b2p = gettransborder (b1n);
  b2n = b2p->next;
  b3p = gettransborder (b2p->next);
  b4p = gettransborder (b3p->next);
  b1p = gettransborder (b4p->next);

  r2 = b2p->border->region;
  r4 = b4p->border->region;

  if (r2 == r4)
  {
    if (debug) printf ("r2 e r4 sono la stessa regione\n");
    /* le regioni 2 e 4 sono la stessa regione.
     * ci sono tre casi, ma li tratto tutti allo stesso modo
     * aggiustando alla fine con una chiamata a adjust_isexternalinfo
     * ad esempio se (r2->border != b2p->border)
     * significa che sono su un'isola, che si spezza in due
     */
    assert (b2p->border == b4p->border);
    assert (b2p != b4p);
    if (topo_change (b2p, b4p) == 0) assert (0);
    changes = adjust_isexternalinfo (sketch);
    assert (changes <= 1);
    if (debug) printf ("changes = %d\n", changes);
    /* purtroppo non e' facile sapere quale delle due componenti
     * connesse sara' un buco, tiro a indovinare e poi utilizzo
     * la funzione adjust_isexternalinfo per sistemare
     */
  } else {
    /* la 2 e 4 sono regioni diverse, che ora
     * diventano la stessa!
     */
    r2ext = r4ext = 0;
    if (r2->border == b2p->border) r2ext = 1;
    if (r4->border == b4p->border) r4ext = 1;
    assert (r2ext || r4ext);     /* almeno per una sono sul bordo esterno */
    /* devo rimappare la regione r4 nella r2 */
    if (regionunion (r2, r4, sketch) != r2)
    {
      printf ("errore fatale in regionunion\n");
      assert (0);
    }
    assert (b2p->border != b4p->border);
    if (debug) printf ("redefining border\n");
    if (topo_change (b2p, b4p) == 0) assert (0);
    if (debug) printf ("adjusting...\n");

    changes = adjust_isexternalinfo (sketch);
    assert (changes <= 1);
  }
  if (debug) checkconsistency (sketch);

  /*
   * ora devo fare il merge degli archi
   * ed eliminare alcune sponde
   */

  assert (b2p->orientation == b2p->next->orientation);
  assert (b4p->orientation == b4p->next->orientation);
  arc12 = b2p->info;
  arc41 = b1p->info;
  if (debug)
  {
    printf ("before merging: arc12 = %d\n", arc12->tag);
    printf ("before merging: arc41 = %d\n", arc41->tag);
  }
  if (b2p->orientation == 1)
    arcleft = mergearcs (arc12, arc41, sketch);
  else
    arcleft = mergearcs (arc41, arc12, sketch);

  if (debug) printf ("after one merge, %d, %d\n", arcleft->tag,
             arcleft->endpoints);
  b1p->info = b2p->info = arcleft;
  // b4p->next->info = b1p->next->info = arcleft;   // serve?

  if (b1p->next != b1p)
  {
    if (b1p->next == b3p) b3p = b1p;
    removeborder (b1p->next);
  }
  if (b2p->next != b2p)
  {
    if (b2p->next == b4p) b4p = b2p;
    removeborder (b2p->next);
  }
  if (b1p->orientation < 0)      /* anticipo causa assert in mergearcs */
      arcleft->regionright = b1p;
    else
      arcleft->regionright = b2p;

  if (debug) checkconsistency (sketch);

  arc23 = b3p->info;
  arc34 = b4p->info;
  if (debug)
  {
    printf ("before merging: arc23 = %d\n", arc23->tag);
    printf ("before merging: arc34 = %d\n", arc34->tag);
  }
  assert (arc23->regionright->info == arc23);
  assert (arc34->regionright->info == arc34);
  if (b4p->orientation == 1)
    arcright = mergearcs (arc34, arc23, sketch);
  else
    arcright = mergearcs (arc23, arc34, sketch);
  if (debug)
  {
    printf ("after merging: arcright = %d\n", arcright->tag);
  }

  b3p->info = b4p->info = arcright;

  if (b3p->next != b3p) removeborder (b3p->next);
  if (b4p->next != b4p) removeborder (b4p->next);
  if (b3p->orientation > 0)
  {
    arcright->regionleft = b3p;
    arcright->regionright = b4p;
  } else {
    arcright->regionright = b3p;
    arcright->regionleft = b4p;
  }
  assert (arcright->regionleft == b3p || arcright->regionleft == b4p);
  assert (arcright->regionright == b3p || arcright->regionright == b4p);

  if (debug) checkconsistency (sketch);
  if (debug) printsketch (sketch);
  if (debug) printf ("esco da taglia_nodo\n");
}

/*
 * trasformazione topologica corrispondente ad una retta
 * che attraversa un incrocio
 */

void 
triple_switch (struct border *b1)
{
  struct border *b[4], *bt[4];
  struct border *btn[4];
  struct border *btnt[4];
  struct arc *a[4];
  int i;

  b[1] = b1;
  b[2] = b[1]->next;
  b[3] = b[2]->next;

  assert (b[1] == b[3]->next);

  for (i = 1; i <= 3; i++)
  {
    a[i] = b[i]->info;
    assert (a[i]->depthsdim == 1);
    a[i]->depths[0] += get_d_increase_across_node (a[i], 1) +
                       get_d_increase_across_node (a[i], -1);
    bt[i] = gettransborder (b[i]);
    btn[i] = bt[i]->next;
    btnt[i] = gettransborder (btn[i]);
  }

  btnt[0] = btnt[3];

  for (i = 1; i <= 3; i++)
  {
    ensurecanremoveborder (bt[i]);
    topo_change_l (prevborder(bt[i]), bt[i]);
    b[i]->orientation *= -1;
    bt[i]->orientation *= -1;
    bt[i]->border = btnt[i-1]->border;
    topo_change_l (btnt[i-1], bt[i]);
  }
}

/*
 * rimozione di una regione circondata da un arco senza nodi
 * (un s1)
 */

void 
remove_s1 (struct border *b, struct sketch *sketch)
{
  struct region *rin, *rout;
  struct borderlist *bl, *btl, *bltemp;
  struct border *bt;
  struct arc *arc;

  assert (b->next == b);
  bl = b->border;
  assert (bl->next == 0);
  rin = bl->region;
  assert (rin->border == bl);

  bt = gettransborder (b);
  assert (bt->next == bt);

  btl = bt->border;
  rout = btl->region;
  assert (rout->border != btl);    /* non puo' essere il bordo esterno */
  assert (rout->border->next);     /* deve esserci almeno un buco */

  arc = b->info;
  assert (arc->endpoints == 0);

  rin->border = 0;
  removeregion (rin, sketch);
  freeborderlist (bl);

  for (bltemp = rout->border; bltemp->next; bltemp = bltemp->next)
  {
    if (bltemp->next == btl)
    {
      bltemp->next = btl->next;    /* rimuovo btl dalla lista */
      btl->next = 0;
      freeborderlist (btl);
      break;
    }
  }

  removearc (arc, sketch);
}

/*
 * unisci una coppia di cuspidi in un nuovo nodo
 */

void
join_cusps (struct border *cusp1, int cusp1pos, 
            struct border *cusp2, int cusp2pos, 
            struct sketch *sketch)
{
  struct borderlist *bl1, *bl2, *newbl;
  struct region *newr;
  struct border *btemp;

  spezza_bordo (cusp1, cusp1pos, sketch);
  if (debug) printf ("dopo il primo spezza_bordo\n");
  if (debug) printsketch (sketch);
  if (cusp1 == cusp2)
  {
    /* devo ricalcolare cusp2 e cusp2pos tenendo conto
     * di quanto fatto da 'spezza_bordo'
     */
    if (cusp1->info->endpoints == 1)
    {
      if (debug) printf ("era un S1\n");
      if (cusp2pos > cusp1pos) cusp2pos -= cusp1pos + 1;
        else cusp2pos += cusp1->info->cusps - cusp1pos;
    } else {
      if (cusp2pos > cusp1pos)
      {
        cusp2pos -= cusp1pos + 1;
      } else {
        cusp2 = cusp2->next;
      }
    }
  }
  if (debug)
  {
    printf ("cusp2pos = %d\n", cusp2pos);
  }
  spezza_bordo (cusp2, cusp2pos, sketch);
  if (cusp1 == cusp2) cusp1 = cusp2->next;
  if (debug) printf ("dopo il secondo spezza_bordo\n");
  if (debug) printsketch (sketch);

  bl1 = cusp1->border;
  bl2 = cusp2->border;
  if (bl1 == bl2)         /* la regione si spezza in due */
  {
    if (debug) printf ("la regione si spezza in due\n");
    newr = newregion (sketch);
    newbl = newborderlist (newr);
    topo_change_l (cusp1, cusp2);
    if (cusp1->next == cusp1) cusp1->info->endpoints = 1;
    if (cusp2->next == cusp2) cusp2->info->endpoints = 1;
    redefineborder (cusp2, newbl);
    bl1->sponda = cusp1;
    newbl->sponda = cusp2;
    if (debug)
    {
      printf ("arco di cusp2: %d\n", cusp2->info->tag);
      btemp = gettransborder (cusp2);
      printf ("arco di futuro b2p: %d\n", btemp->info->tag);
      btemp = gettransborder (cusp2->next);
      printf ("arco di futuro b3p: %d\n", btemp->info->tag);
    }
  } else {                /* rimane una sola regione */
    if (bl2 != bl2->region->border)   /* e' un buco */
    {
      redefineborder (cusp2, bl1);
      bl1 = bl2;
    } else {
      redefineborder (cusp1, bl2);
    }
    topo_change_l (cusp1, cusp2);
    if (cusp1->next == cusp1) cusp1->info->endpoints = 1;
    if (cusp2->next == cusp2) cusp2->info->endpoints = 1;
    extractborderlist (bl1);
    bl1->sponda = 0;
    freeborderlist (bl1);
  }
}

/*
 * suddivido un arco in corrispondenza di una cuspide,
 * eventualmente creo i border necessari
 */

void
spezza_bordo (struct border *cusp, int cusppos, struct sketch *sketch)
{
  int *newdepths, i, j;
  struct arc *arc, *arcnew;
  struct border *bt, *btnew;
  struct border *cuspnew;

  arc = cusp->info;
  assert (cusp->orientation < 0);
  if (arc->endpoints == 0)  /* questo e' un s1 */
  {
    /* in questo caso non si creano altri archi e bordi
     * anche il numero di valori di d rimane uguale, ma non
     * serve ripetere il primo valore
     */
    assert (arc->depthsdim == arc->cusps + 1);
    arc->cusps--;
    arc->depthsdim--;
    newdepths = (int *) malloc (arc->depthsdim * sizeof (int));
    for (i = 0; i < arc->depthsdim; i++)
    {
      j = i + cusppos + 1;
      if (j >= arc->depthsdim) j -= arc->depthsdim;
      newdepths[i] = arc->depths[j];
    }
    free (arc->depths);
    arc->depths = newdepths;
    arc->endpoints = 1;
  } else {
    /* in questo caso si crea un nuovo arco e due nuove sponde */
    arcnew = newarc (sketch);
    bt = gettransborder (cusp);
    btnew = newborder (bt->border);
    btnew->next = bt->next;
    bt->next = btnew;
    btnew->info = arcnew;
    arcnew->regionleft = btnew;
    btnew->orientation = bt->orientation;

    cuspnew = newborder (cusp->border);
    cuspnew->next = cusp->next;
    cusp->next = cuspnew;
    cusp->info = arcnew;
    arcnew->regionright = cusp;
    cuspnew->info = arc;
    arc->regionright = cuspnew;
    cuspnew->orientation = cusp->orientation;
    arcnew->depthsdim = arc->depthsdim - cusppos - 1;
    arcnew->cusps = arc->cusps - cusppos - 1;
    arcnew->dvalues = arcnew->depthsdim;
    newdepths = (int *) malloc (arcnew->depthsdim * sizeof (int));
    arcnew->depths = newdepths;
    /* devo definire endpoints, che ora e' 2 */
    arc->endpoints = arcnew->endpoints = 2;
    for (i = 0; i < arcnew->depthsdim; i++)
    {
      newdepths[i] = arc->depths[i + cusppos + 1];
    }
    arc->depthsdim = cusppos + 1;
    arc->cusps = cusppos;
    arc->dvalues = arc->depthsdim;
    newdepths = (int *) malloc (arc->depthsdim * sizeof (int));
    for (i = 0; i < arc->depthsdim; i++)
    {
      newdepths[i] = arc->depths[i];
    }
    free (arc->depths);
    arc->depths = newdepths;
  }
  if (debug) checkconsistency (sketch);
}

/*
 * rimuovo una regione anulare
 * eventuali buchi della regione interna
 * vengono aggiunti ai buchi della regione esterna
 */

void
remove_annulus (struct region *r, struct sketch *sketch)
{
  struct region *extr, *inr;
  struct arc *extarc, *inarc;
  struct borderlist *extbl, *inbl, *extblt, *inblt;
  struct borderlist *inbltholes, *bl;
  struct border *extb, *inb, *extbt, *inbt;

  if (debug) printf ("entering remove_annulus, region %d\n", r->tag);
  if (debug) printsketch (sketch);
  extbl = r->border;
  extb = extbl->sponda;
  inbl = extbl->next;
  inb = inbl->sponda;
  assert (extbl->next);
  extarc = extb->info;
  inarc = inb->info;
  extbt = gettransborder (extb);
  inbt = gettransborder (inb);
  freeborderlist (extbl);        /* cancella anche il buco */

  extblt = extbt->border;
  extr = extblt->region;
  extblt = extractborderlist (extblt);  /* rimuovo la sponda esterna */
  freeborderlist (extblt);

  inblt = inbt->border;
  inr = inblt->region;
  assert (inr->border == inblt);  /* questo deve essere il bordo esterno */
  inbltholes = inblt->next;
  redefineregion (inr, extr);
  inblt->next = 0;
  freeborderlist (inblt);
  inr->border = r->border = 0;
  removeregion (r, sketch);
  removeregion (inr, sketch);

  if (debug) printsketch (sketch);
  if (debug) checkconsistency (sketch);
  removearc (inarc, sketch);
  removearc (extarc, sketch);

  assert (extr->border);
  for (bl = extr->border; bl; bl = bl->next)
  {
    if (bl->next == 0)
    {
      bl->next = inbltholes;
      break;
    }
  }
}

/*
 * rimuovo una regione a "coda di rondine"
 * oppure senza cuspidi (regola n6)
 */

void
remove_ear (struct region *r, struct sketch *sketch)
{
  struct borderlist *bl, *blext;
  struct border *b, *bext, *bextp, *bextt;
  struct arc *arcear, *arc1, *arc2;
  int cnum, ori, dincr;

  bl = r->border;
  assert (bl->next == 0);
  b = bl->sponda;
  assert (b->next == b);
  ori = b->orientation;
  arcear = b->info;
  cnum = arcear->cusps;
  assert ((ori > 0?(cnum == 2):(cnum == 0)));
  bext = gettransborder (b);
  bextp = prevborder (bext);
  assert (bextp != bext);
  blext = bext->border;

  freeborderlist (bl);
  r->border = 0;
  removeregion (r, sketch );

  bext->info = 0;
  if (ori > 0) arcear->regionright = 0;
    else arcear->regionleft = 0;
  bextp = removeborder (bext);
  removearc (arcear, sketch);

  arc1 = bextp->info;
  bext = bextp->next;
  arc2 = bext->info;
  bextt = gettransborder (bext);

  assert (ori*bext->orientation > 0);
  assert (ori*bextp->orientation > 0);
  if (ori > 0)
  {
    arc1 = mergearcsc (arc1, arc2, 0, sketch);
  } else {
    dincr = arc1->depths[0] - arc2->depths[arc2->dvalues - 1];
    arc1 = mergearcsc (arc2, arc1, dincr, sketch);
  }
  bextp->info = bext->info = 0;
  bextt->info = bextt->next->info = 0; 
  if (bextt->next != bextt) bextt = removeborder (bextt->next);
  if (bextp != bext) bextp = removeborder (bext); 
  bextp->info = bextt->info = arc1;
  if (ori > 0)
  {
    arc1->regionleft = bextp;
    arc1->regionright = bextt;
  } else {
    arc1->regionright = bextp;
    arc1->regionleft = bextt;
  }
  if (debug) checkconsistency (sketch);
}

/*
 * esegui le modifiche necessarie per le regole cn2l/cn2r
 * bisogna cambiare una profondita' e richiamare taglia_nodo
 * 
 * orib = 1  -> cuspide sul davanti della falda
 * orib = -1 -> cuspide dietro alla falda
 */

void
applyrulecn2 (struct border *b1n, struct arc *arc, 
              int ori, int orib, struct sketch *sketch)
{
  int index;

  if (debug) checkconsistency (sketch);
  assert (arc->cusps > 0);
  index = 0;
  if (ori < 0) index = arc->dvalues - 1;
  arc->depths[index] += 2*orib;
  taglia_nodo (b1n, sketch);
}

/*
 * esecuzione mossa cn3: rimozione di una cuspide
 */

void
remove_cusp (struct region *r, struct sketch *sketch)
{
  struct arc *arc, *arccusp;
  struct borderlist *bl;
  struct border *bcusp, *bseg, *bsegt, *bsegtn, *outnw;
  int i, *newdepths;

  bl = r->border;
  assert (bl->next == 0);
  assert (bl->sponda);

  bcusp = bl->sponda;
  if (bcusp->info->cusps < 1) bcusp = bcusp->next;
  arccusp = bcusp->info;
  assert (arccusp->cusps == 1);
  bseg = bcusp->next;
  assert (bseg->next == bcusp);

  bsegt = gettransborder (bseg);
  bsegtn = bsegt->next;
  assert (bsegtn->orientation == bcusp->orientation);
  arc = bsegtn->info;
  /* ora devo spostare le informazioni sulla cuspide sull'arco arc */
  newdepths = (int *) malloc ((arc->depthsdim + 1)*sizeof (int));
  if (bsegtn->orientation > 0)
  {
    newdepths[0] = arc->depths[0] - arccusp->depths[1]
                                  + arccusp->depths[0];
    for (i = 0; i < arc->depthsdim; i++)
      newdepths[i+1] = arc->depths[i];
  } else {
    for (i = 0; i < arc->depthsdim; i++)
      newdepths[i] = arc->depths[i];

    newdepths[arc->depthsdim] = newdepths[arc->depthsdim - 1]
       - arccusp->depths[0] + arccusp->depths[1];
  }
  arc->depthsdim++;
  arc->cusps++;
  arc->dvalues++;
  free (arc->depths);
  arc->depths = newdepths;
  /* devo togliere la cuspide da arccusp, non sto a liberare/allocare */
  arccusp->cusps = 0;
  arccusp->depthsdim = 1;
  arccusp->dvalues = 1;
  if (bsegtn->orientation < 0)
    arccusp->depths[0] = arccusp->depths[1];

  if (debug) checkconsistency (sketch);

  outnw = rimuovi_losanga (r, sketch);
  taglia_nodo (outnw, sketch);
}
