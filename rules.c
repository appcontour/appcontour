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
  char *chpt, *endrule;
  extern int heisemberg;
  extern int quiet;

  canonify (sketch);
  if (appcontourcheck (sketch, 1, 0) == 0)
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
      rcount = strtol (chpt, &endrule, 10);
      if (rcount < 1) rcount = 1;
      if (*endrule == ':')
      {
        *endrule++ = 0;
        heisemberg = atoi (endrule);
        if (!quiet) fprintf (stderr, "transfer of islands requested: %d\n", heisemberg);
      }
    }
  }
  free_connected_components (sketch);
  res = 0;
  if (strcasecmp (rule, "n1") == 0) res = rule_n14 (sketch, 1, rcount);
  else if (strcasecmp (rule, "n2") == 0) res = rule_n14 (sketch, 2, rcount);
  else if (strcasecmp (rule, "n3") == 0) res = rule_n14 (sketch, 3, rcount);
  else if (strcasecmp (rule, "n4") == 0) res = rule_n14 (sketch, 4, rcount);
  else if (strcasecmp (rule, "n5") == 0) res = rule_n5 (sketch, rcount);
  else if (strcasecmp (rule, "cr2") == 0) res = rule_cr2 (sketch, rcount);
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
  else if (strcasecmp (rule, "cr3l") == 0) res = rule_cr3l (sketch, rcount);
  else if (strcasecmp (rule, "cr3r") == 0) res = rule_cr3r (sketch, rcount);
  else if (strcasecmp (rule, "cr1") == 0) res = rule_cr1 (sketch, rcount);
  else if (strcasecmp (rule, "cr1b") == 0) res = rule_cr1b (sketch, rcount);
  else if (strcasecmp (rule, "cr4l") == 0) res = rule_cr4l (sketch, rcount);
  else if (strcasecmp (rule, "cr4r") == 0) res = rule_cr4r (sketch, rcount);
  else if (strcasecmp (rule, "cr4lb") == 0) res = rule_cr4lb (sketch, rcount);
  else if (strcasecmp (rule, "cr4rb") == 0) res = rule_cr4rb (sketch, rcount);
  else if (strcasecmp (rule, "t1") == 0) res = rule_t1 (sketch, rcount);
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
  exitcode = testsinglerule ("CR2", rule_cr2, exitcode, sketch);
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
  exitcode = testsinglerule ("CR3L", rule_cr3l, exitcode, sketch);
  exitcode = testsinglerule ("CR3R", rule_cr3r, exitcode, sketch);
  exitcode = testsinglerule ("CR1", rule_cr1, exitcode, sketch);
  exitcode = testsinglerule ("CR1B", rule_cr1b, exitcode, sketch);
  exitcode = testsinglerule ("CR4L", rule_cr4l, exitcode, sketch);
  exitcode = testsinglerule ("CR4R", rule_cr4r, exitcode, sketch);
  exitcode = testsinglerule ("CR4LB", rule_cr4lb, exitcode, sketch);
  exitcode = testsinglerule ("CR4RB", rule_cr4rb, exitcode, sketch);
  exitcode = testsinglerule ("T1", rule_t1, exitcode, sketch);
  /* commented out because there is an infinite loop in some cases */
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
      taglia_nodo (outnw, sketch, 0, 0);
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
rule_cr2 (struct sketch *sketch, int rcount)
{
  return (rule_cn1_cr2 (sketch, rcount, 1));
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
          /* OLD join_cusps (cusp1, cusp1pos, cusp2, cusp2pos, sketch); */
          pinch_arcs (&cusp1, cusp1pos, &cusp2, cusp2pos, sketch, 1);
          if (debug) checkconsistency (sketch);
          taglia_nodo (cusp1->next, sketch, 0, 0);
          return (1);
        }
      }
    }
  }
  return (0);
}

int
rule_t1 (struct sketch *sketch, int rcount)
{
  int onlytest = 0;
  struct region *r, *r1, *r2;
  struct border *sp1, *sp2;
  struct borderlist *bl1, *bl2;
  struct arc *a1, *a2;
  int i, ok, d1min, d1max, changes;
  extern int heisemberg;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* nel caso piu' frequente di una isola, la sposta nell'altra
   * regione
   */
  if (heisemberg < 0) heisemberg = 1;

  /* cerca una regione candidata, che deve avere
   * due archi orientati positivamente con regioni affacciate
   * diverse e |d1 - d2| = 1 per almeno un paio di archi tra due
   * cuspidi sui rispettivi archi.  Inoltre ci deve essere almeno una
   * isola da trasferire nelle regioni affacciate.
   */

  for (r = sketch->regions; r; r = r->next)
  {
    //fprintf (stderr, "testing region %d\n", r->tag);
    if (r->border->sponda == 0) continue;  /* questa e' la regione esterna */
    for (bl1 = r->border; bl1; bl1 = bl1->next)
    {
      for (bl2 = bl1; bl2; bl2 = bl2->next)
      {
        sp1 = bl1->sponda;
        do {
          if (sp1->orientation <= 0) continue;
          a1 = sp1->info;
          d1min = d1max = a1->depths[0];
          for (i = 1; i < a1->dvalues; i++)
          {
            if (a1->depths[i] < d1min) d1min = a1->depths[i];
            if (a1->depths[i] > d1max) d1max = a1->depths[i];
          }
          r1 = a1->regionright->border->region;
          if (bl2 == bl1) sp2 = sp1; else sp2 = bl2->sponda;
          do {
            //fprintf (stderr, "region %d, d1min=%d d1max=%d\n", r->tag, d1min, d1max);
            //fprintf (stderr, "d = %d\n", sp2->info->depths[0]);
            if (sp2->orientation <= 0) continue;
            if (sp2 == sp1) continue;
            a2 = sp2->info;
            //fprintf (stderr, "a1: %d, a2: %d\n", a1->tag, a2->tag);
            ok = 0;
            for (i = 0; i < a2->dvalues; i++)
            {
              if (a2->depths[i] <= d1max + 1 && a2->depths[i] >= d1min - 1)
              {
                if (d1min < d1max || a2->depths[i] != d1max)
                {
                  ok = 1;
                  break;
                }
              }
            }
            if (ok == 0) continue;
            r2 = a2->regionright->border->region;
            if (r2 == r1) continue;
            if (r1->border->next == 0 && r2->border->next == 0) continue;
            if (debug) printf ("spostamento di isole, regione %d\n", r->tag);
            /* ho trovato un possibile trasferimento di isole */
            if (rcount-- <= 1)
            {
              if (onlytest) return (1);
              /* sposta isole... */
              if (heisemberg < 0)
              {
                fprintf (stderr, "You need to indicate which isles to transfer, use\n");
                fprintf (stderr, "switch \"--ti <int>\", or the special rule names ");
                fprintf (stderr, "T1::<int>, T1:2:<int>, ...\n");
                fprintf (stderr, "<int> is interpreted as a bit mask\n");
                return (1);
              }
              transfer_isles (a1->regionright->border, a2->regionright->border,
                              heisemberg);
              changes = adjust_isexternalinfo (sketch);
              if (changes) fprintf (stderr, "changes: %d\n", changes);
              return (1);
            }
          } while (sp2 = sp2->next, sp2 != bl2->sponda);
        } while (sp1 = sp1->next, sp1 != bl1->sponda);
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
  return (rule_cn1_cr2 (sketch, rcount, 0));
}

int
rule_cn1_cr2 (struct sketch *sketch, int rcount, int iscr2)
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
   * per la regola cr2 non ci devono essere cuspidi!
   */

  ori = 1;
  if (iscr2) ori = -1;
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
    if (iscr2)
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

int
rule_cr1 (struct sketch *sketch, int rcount)
{
  return (rule_cr11b (sketch, rcount, 0));
}

int
rule_cr1b (struct sketch *sketch, int rcount)
{
  return (rule_cr11b (sketch, rcount, 1));
}

int
rule_cr11b (struct sketch *sketch, int rcount, int isback)
{
  struct region *r;
  struct borderlist *bl, *blp;
  struct border *bpstart, *bppstart, *bp, *bpp;
  int i, j, catafound, onlytest = 0, orib, db;
  struct arc *arc, *arct = 0;
  int ip2, numchecks;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* per applicare la regola cr1 ci vuole una sequenza di tre valori
   * di "d" (due cuspidi) del tipo d-+  (es 1 0 1) oppure d+- nel caso
   * di "back", con orientazione negativa; in entrambi i casi serve un 
   * arco positivo con valore "d-1" (distinto dall'eventuale valore uguale tra
   * le due cuspidi)
   */

  orib = 1;
  if (isback) orib = -1;

  for (r = sketch->regions; r; r = r->next)
  {
    for (bl = r->border; bl; bl = bl->next)
    {
      bpstart = bl->sponda;
      if (bpstart == 0) continue;
      bp = bpstart;
      do {
        if (bp->orientation > 0) continue;
        arc = bp->info;
        if (arc->cusps < 2) continue;
        numchecks = arc->cusps - 1;
        if (arc->endpoints == 0) numchecks = arc->cusps;
        for (i = 0; i < numchecks; i++)
        {
          db = arc->depths[i];
          if (arc->depths[i+1] != db - orib) continue;
          ip2 = i+2;
          if (ip2 > arc->cusps) ip2 = 1;
          if (arc->depths[ip2] != db) continue;
          /* trovata sequenza di cuspidi, devo vedere se c'e' un arco
           * catalizzante (con d = db - 1)
           */
          catafound = 0;
          for (blp = r->border; blp; blp = blp->next)
          {
            bppstart = blp->sponda;
            if (bppstart == 0) continue;
            bpp = bppstart;
            do {
              if (bpp->orientation < 0) continue;
              arct = bpp->info;
              for (j = 0; j < arct->dvalues; j++)
              {
                if (arct->depths[j] != db - 1) continue;
                if (bpp == bp && j == i + 1) continue;
                catafound = 1;
                break;
              }
              if (catafound) break;
            } while (bpp = bpp->next, bpp != bppstart);
            if (catafound) break;
          }
          if (! catafound) continue;
          if (rcount-- <= 1)
          {
            if (onlytest) return (1);
            fprintf (stderr, "Trovata regione: %d, arc %d, cataarc %d\n",
                             r->tag, arc->tag, arct->tag);
            /* applica la regola, che in questo caso e' molto semplice */
            arc->depths[i+1] += 2*orib;
            if (arc->endpoints == 0 && i+1 == arc->cusps) arc->depths[0] += 2*orib;
            return (1);
          }
        }
      } while (bp = bp->next, bp != bpstart);
    }
  }
  return (0);
}

int
rule_cr4l (struct sketch *sketch, int rcount)
{
  return (rule_cr4lr (sketch, rcount, 1, 0));
}

int
rule_cr4r (struct sketch *sketch, int rcount)
{
  return (rule_cr4lr (sketch, rcount, 0, 0));
}

int
rule_cr4lb (struct sketch *sketch, int rcount)
{
  return (rule_cr4lr (sketch, rcount, 1, 1));
}

int
rule_cr4rb (struct sketch *sketch, int rcount)
{
  return (rule_cr4lr (sketch, rcount, 0, 1));
}

int
rule_cr4lr (struct sketch *sketch, int rcount, int isleft, int isback)
{
  int ori, orib, onlytest = 0;
  struct region *r, *newr;
  struct borderlist *bl, *newbl;
  struct border *bp, *bpstart, *newb;
  struct arc *arc, *arcs1;
  int i, j, db, testnum, ip2, ip3;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}
  ori = orib = 1;
  if (! isleft) ori = -1;     /* ori > 0 if 'isleft' */
  if (isback) orib = -1;

  for (r = sketch->regions; r; r = r->next)
  {
    for (bl = r->border; bl; bl = bl->next)
    {
      bpstart = bl->sponda;
      if (bpstart == 0) continue;
      bp = bpstart;
      do {
        if (bp->orientation > 0) continue;
        arc = bp->info;
        if (arc->cusps < 3) continue;
        testnum = arc->cusps - 2;
        if (arc->endpoints == 0) testnum = arc->cusps;
        for (i = 0; i < testnum; i++)
        {
          db = arc->depths[i];
          if (arc->depths[i+1] != db + orib) continue;
          ip2 = i + 2;
          if (ip2 > arc->cusps) ip2 -= arc->cusps;
          if (arc->depths[ip2] != db + orib - ori) continue;
          ip3 = i + 3;
          if (ip3 > arc->cusps) ip3 -= arc->cusps;
          if (arc->depths[ip3] != db - ori) continue;
          /* 
	   * trovata sequenza di cuspidi, si tratta in sostanza
	   * di applicare una regola C2 che creera' una isola S^1
	   * cui posso far attraversare l'arco oltre la cuspide
	   * non coinvolta: CR4R = C2 N4^{-1} N3.
           */
          if (rcount-- <= 1)
          {
            if (onlytest) return (1);
            fprintf (stderr, "Trovata regione: %d, arc %d, i %d\n",
                             r->tag, arc->tag, i);
            /* devo creare un buco a forma di s1 nella regione e
             * accorciare il vettore dei valori di d
             */
            arcs1 = newarc (sketch);
            arcs1->depths = (int *) malloc (sizeof (int));
            arcs1->depthsdim = 1;
            arcs1->cusps = 0;
            arcs1->dvalues = 1;
            arcs1->endpoints = 0;
            arcs1->depths[0] = db + (orib - ori)/2 - 1;
            newr = newregion (sketch);
            newbl = newborderlist (newr);
            newb = newborder (newbl);
            newbl->sponda = newb;
            newb->orientation = -1;
            newb->info = arcs1;
            arcs1->regionright = newb;

            newbl = newborderlist (r);
            newb = newborder (newbl);
            newbl->sponda = newb;
            newb->orientation = 1;
            newb->info = arcs1;
            arcs1->regionleft = newb;

            /* ora mi occupo delle cuspidi */
            /* devo eliminare i valori di d in posizione i+1 e i+2 */

            assert (i <= arc->cusps - 1);
            if (i < arc->cusps - 2)
            {
              for (j = i + 1; j <= arc->cusps; j++)
                arc->depths[j] = arc->depths[j+2];
            } else {
              assert (arc->endpoints == 0);
              arc->depths[0] = arc->depths[arc->cusps - 2];
              if (i == arc->cusps - 1) arc->depths[1] = arc->depths[arc->cusps - 1];
            }
            arc->cusps -= 2;
            arc->dvalues -= 2;
            return (1);
          }
        }
      } while (bp = bp->next, bp != bpstart);
    }
  }
  return (0);
}

int
rule_cr3l (struct sketch *sketch, int rcount)
{
  return (rule_cr3lr (sketch, rcount, 1));
}

int
rule_cr3r (struct sketch *sketch, int rcount)
{
  return (rule_cr3lr (sketch, rcount, -1));
}

/* ori = 1 -> cr3 del paper, altrimenti e' quella simmetrica */

int
rule_cr3lr (struct sketch *sketch, int rcount, int ori)
{
  int onlytest = 0, d, j;
  struct region *r;
  struct borderlist *bl1, *bl2;
  struct border *b1, *b2, *b1n;
  struct arc *arc1, *arc1n, *arc2;

  if (rcount < 0) {onlytest = 1; rcount *= -1;}

  /* cerca una regione candidata, che deve avere
   * due archi consecutivi orientati positivamente con
   * d, d+2 e un terzo arco (positivo) con d+1
   */

  for (r = sketch->regions; r; r = r->next)
  {
    if (debug) printf ("region: %d\n", r->tag);
    for (bl1 = r->border; bl1; bl1 = bl1->next)
    {
      b1 = bl1->sponda;
      if (b1 == 0) continue;
      do {
        if (b1->orientation < 0) continue;
        if ((b1n = b1->next) == b1) continue;
        if (b1n->orientation < 0) continue;
        arc1 = b1->info;
        arc1n = b1n->info;
        d = arc1->depths[arc1->dvalues - 1];
        if (arc1n->depths[0] != d + ori*2) continue;
        if (get_d_increase_across_node (arc1, 1) != ori - 1) continue;
        if (get_d_increase_across_node (arc1n, -1) != - ori - 1) continue;

        /* we found a good candidate for the consecutive arcs */

        for (bl2 = r->border; bl2; bl2 = bl2->next)
        {
          b2 = bl2->sponda;
          if (b2 == 0) continue;
          do {
	    if (b2->orientation < 0) continue;
            arc2 = b2->info;
            for (j = 0; j < arc2->dvalues; j++)
            {
              if (arc2->depths[j] != d + ori) continue;

              /* ho trovato una posizione di applicabilita' */
              if (rcount-- <= 1)
              {
                if (onlytest) return (1);
                fprintf (stderr, "Trovata regione: %d, arc1/n/2 = ", r->tag);
                fprintf (stderr, "%d %d %d, j = %d\n", 
                  arc1->tag, arc1n->tag, arc2->tag, j);
                applyrulecr3 (b1, b2, j, sketch);
                if (debug) checkconsistency (sketch);
                return (1);
              }
            }
          } while (b2 = b2->next, b2 != bl2->sponda);
        }
      } while (b1 = b1->next, b1 != bl1->sponda);
    }
  }
  return (0);
}

/*
 * Questa funzione serve ad applicare una mossa inversa
 * per il merge di due archi (inversa di N1-N4 oppure C2)
 */

/* local prototype */
int common_work_mergearcs (struct sketch *s, 
			   struct border *bp1, struct border *bp2, 
                           int a1pos, int a2pos, int rule);

#define INV_N1 1
#define INV_N2 2
#define INV_N3 3
#define INV_N3bis 4
#define INV_N4 5
#define INV_C2 6
static char *invmergerules[] = {
  "",
  "INV N1",
  "INV N2",
  "INV N3",
  "INV N3bis",
  "INV N4",
  "INV C2",
  0};

int
list_mergearcs (struct sketch *s, struct region *r,
	struct arc *a1, struct arc *a2, int a1l, int a2l)
{
  struct borderlist *bl;
  struct border *bp;
  int imax, i, res;

  int count = 0;
  if (r == 0) {
    for (r = s->regions; r; r = r->next) {
      printf ("Region %d:\n", r->tag);
      count += list_mergearcs (s, r, 0, 0, -1, -1);
    }
    return (count);
  }

  if (a1 == 0 || (a1l >= 0 && a2 == 0)) {
    for (bl = r->border; bl; bl = bl->next)
    {
      if (bl->sponda == 0) continue;
      bp = bl->sponda;
      do {
        count += list_mergearcs (s, r, 
		(a1)?a1:bp->info, (a1)?bp->info:0, a1l, -1);
	bp = bp->next;
      } while (bp != bl->sponda);
    }
    return (count);
  }

  assert (s && r && a1);
  if (a1l < 0) {
    imax = a1->dvalues;
    for (i = 0; i < imax; i++) {
      count += list_mergearcs (s, r, a1, 0, i, -1);
    }
    return (count);
  }

  assert (a2 && a1l >= 0);
  if (a2l < 0) {
    imax = a2->dvalues;
    for (i = 0; i < imax; i++) {
      count += list_mergearcs (s, r, a1, a2, a1l, i);
    }
    return (count);
  }

  assert (a2l >= 0);
  res = apply_mergearcs (s, r, a1, a2, a1l, a2l, 1);
  if (res) {
    printf ("-r %d -a %d:%d -a %d:%d (%s)\n", 
      r->tag, a1->tag, a1l, a2->tag, a2l, invmergerules[res]);
    return (1);
  }

  return (0);
}

int
apply_mergearcs (struct sketch *s, struct region *r,
	struct arc *a1, struct arc *a2, int a1l, int a2l, int test)
{
  struct borderlist *bl;
  struct border *bp;
  struct border *bp1 = 0;
  struct border *bp2 = 0;
  int d1 = -1000;
  int d2 = -1000;
  int deltad, res;
  extern int verbose;

  assert (s && r && a1 && a2 && a1l >= 0 && a2l >= 0);

  /* consistency check about the number of dvalues */
  if (a1l >= 0 && a1l < a1->dvalues) d1 = a1->depths[a1l];
  if (a2l >= 0 && a2l < a2->dvalues) d2 = a2->depths[a2l];

  if (d1 < 0) fprintf (stderr, "invalid subarc (%d) for first arc\n", a1l);
  if (d2 < 0) fprintf (stderr, "invalid subarc (%d) for second arc\n", a2l);
  if (d1 > d2) { /* ensure d1 <= d2 */
    if (test == 0) fprintf (stderr, 
	"d value of first arc cannot exceed d value of the second\n");
    return (0);
  }

  /* check if given arcs bound the given region */
  for (bl = r->border; bl; bl = bl->next)
  {
    if (bl->sponda == 0) continue;
    bp = bl->sponda;
    do {
      if (bp->info == a1) bp1 = bp;
      if (bp->info == a2) bp2 = bp;
      bp = bp->next;
    } while (bp != bl->sponda);
  }
  if (bp1 == 0) 
    fprintf (stderr, "Arc %d does not bound given region\n", a1->tag);
  if (bp2 == 0) 
    fprintf (stderr, "Arc %d does not bound given region\n", a2->tag);

  if (! (bp1 && bp2 && d1 >= 0 && d2 >= d1)) return (0);

  if (debug) {
    if (bp1->border == bp2->border) {
      fprintf (stderr, 
	"Arcs belong to the same c.c. of the boundary of region.\n");
      fprintf (stderr, "The region will be splitted in two\n");
    } else {
      fprintf (stderr,
	"Arcs belong to different c.c. of the boundary of region.\n");
      fprintf (stderr, "The number of c.c. will decrease by one\n");
    }
  }

  assert (d1 <= d2);
  deltad = d2 - d1;
  if (bp1->orientation < 0 && bp2->orientation < 0)
  {
    if (verbose) {
      fprintf (stderr, "Negative orientations: can apply inv N1\n");
      fprintf (stderr, "with first arc above second arc\n");
    }
    if (test) return (INV_N1);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N1);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation > 0 && deltad >= 2)
  {
    if (verbose) {
      fprintf (stderr, "Positive orientations, delta d >= 2,\n");
      fprintf (stderr, "can apply inv N4\n");
    }
    if (test) return (INV_N4);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N4);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation > 0 && deltad == 1)
  {
    if (verbose) {
      fprintf (stderr, "Positive orientations, delta d = 1,\n");
      fprintf (stderr, "can apply inv C2\n");
    }
    if (test) return (INV_C2);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_C2);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation > 0 && deltad == 0)
  {
    if (test == 0) {
      fprintf (stderr, "Positive orientations with delta d = 0\n");
      fprintf (stderr, "   Cannot merge arcs\n");
    }
    return (0);
  }

  if (bp1->orientation > 0 && bp2->orientation < 0 && deltad == 0)
  {
    if (verbose) {
      fprintf (stderr, "Pos-neg orientations with delta d = 0\n");
      fprintf (stderr, "can apply inv N3\n");
    }
    /* if (test) return (INV_N3bis);
     */
    if (test) return (0); /* do not report a duplicate of INV_N3 */
    res = common_work_mergearcs (s, bp2, bp1, a2l, a1l, INV_N3);
    return (res);
  }

  if (bp1->orientation > 0 && bp2->orientation < 0 && deltad == 1)
  {
    if (test == 0) {
      fprintf (stderr, "Pos-neg orientations with delta d = 1\n");
      fprintf (stderr, "   Cannot merge arcs\n");
    }
    return (0);
  }

  if (bp1->orientation > 0 && bp2->orientation < 0 && deltad >= 2)
  {
    if (verbose) {
      fprintf (stderr, "Pos-neg orientations with delta d >= 2\n");
      fprintf (stderr, "can apply inv N2\n");
    }
    if (test) return (INV_N2);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N2);
    return (res);
  }

  if (bp1->orientation < 0 && bp2->orientation > 0)
  {
    if (verbose) {
      fprintf (stderr, "Neg-pos orientations\n");
      fprintf (stderr, "can apply inv N3\n");
    }
    if (test) return (INV_N3);
    res = common_work_mergearcs (s, bp1, bp2, a1l, a2l, INV_N3);
    return (res);
  }

  return (0);
}

int
common_work_mergearcs (struct sketch *s, 
		struct border *bp1, struct border *bp2, 
		int a1l, int a2l, int rule)
{
  int res, changes;

  pinch_arcs (&bp1, a1l, &bp2, a2l, s, (rule == INV_C2)?(-1):0);

  if (debug) printf ("dopo topo_change\n");
  if (debug) printsketch (s);

  assert (bp1 != bp2);

  if (rule == INV_C2) {
    taglia_nodo (bp1->next, s, 0, 0);
    res = 1;
  } else {
    res = aggiungi_losanga (bp1, bp2, s);
    changes = adjust_isexternalinfo (s);
    if (debug && changes) 
      printf ("%d changes in adjust_isexternalinfo\n", changes);
    if (debug) printf ("dopo aggiungi_losanga\n");
    if (debug) printsketch (s);
  }
  if (res)
  {
    checkconsistency (s);
    postprocesssketch (s);
    canonify (s);
  }
  return (res);
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
 * contrario di rimuovi_losanga:
 *
 * \bp2  /    \   /
 *  \   /      \ /
 *   \ /        X
 *    X   ==>  / \
 *   / \       \ /
 *  /   \       X
 * /  bp1\     / \
 *
 * si conviene che il valore piu' basso di d
 * venga assegnato alla parte indicata da bp1
 */

int
aggiungi_losanga (struct border *bp1, struct border *bp2,
		struct sketch *s)
{
  struct region *newr;
  struct borderlist *newbl;
  struct border *bp1nt, *bp2nt, *newbp1, *newbp2;
  struct border *newbpt1, *newbpt2;
  struct arc *newa1, *newa2;
  int o1, o2, d1, d2;

  o1 = bp1->orientation;
  assert (o1 == 1 || o1 == -1);
  o2 = bp2->orientation;
  assert (bp2->next->orientation == o1);
  assert (bp1->next->orientation == o2);

  /*
   * la creazione della losanga richiede numerose
   * nuove strutture dati: due archi, una regione
   * e 4 nuove sponde.
   */

  newa1 = newarc (s);
  newa2 = newarc (s);
  newa1->depths = (int *) malloc (sizeof (int));
  newa2->depths = (int *) malloc (sizeof (int));
  newa1->depthsdim = newa2->depthsdim = 1;  /* niente cuspidi */
  newa1->cusps = newa2->cusps = 0;
  newa1->dvalues = newa2->dvalues = 1;
  newa1->endpoints = newa2->endpoints = 2;
  if (o1 > 0) d1 = bp1->info->depths[bp1->info->dvalues - 1];
    else d1 = bp1->info->depths[0];
  if (o2 > 0) d2 = bp2->info->depths[bp2->info->dvalues - 1];
    else d2 = bp2->info->depths[0];
  /* calcolo i due valori di d. Arco 1 passa davanti */
  newa1->depths[0] = d1;   /* questo e' giusto perche' passa davanti */
  newa2->depths[0] = d2 - 2*o1;  /* varia di 2 con segno che dipende */

  newr = newregion (s);
  newbl = newborderlist (newr);
  newbp1 = newborder (newbl);
  newbl->sponda = newbp1;
  newbp2 = newborder (newbl);
  newbp2->next = newbp1;  /* o anche = newbp1->next */
  newbp1->next = newbp2;

  bp1nt = gettransborder (bp1->next);  /* sponda opposta al next */
  bp2nt = gettransborder (bp2->next);
  newbpt1 = newborder (bp1nt->border); /* e' la sponda opposta a bpn1 */
  newbpt1->next = bp1nt->next;  /* da inserire dopo bp1nt */
  bp1nt->next = newbpt1;
  newbp1->info = newbpt1->info = newa1;
  newbp1->orientation = -o1;
  newbpt1->orientation = o1;
  newbpt2 = newborder (bp2nt->border); /* e' la sponda opposta a bpn1 */
  newbpt2->next = bp2nt->next;
  bp2nt->next = newbpt2;
  newbp2->info = newbpt2->info = newa2;
  newbp2->orientation = -o2;
  newbpt2->orientation = o2;

  /* infine sistemo i puntatori alle regioni right/left */
  if (o1 > 0) {
    newa1->regionright = newbp1;
    newa1->regionleft = newbpt1;
  } else {
    newa1->regionleft = newbp1;
    newa1->regionright = newbpt1;
  }
  if (o2 > 0) {
    newa2->regionright = newbp2;
    newa2->regionleft = newbpt2;
  } else {
    newa2->regionleft = newbp2;
    newa2->regionright = newbpt2;
  }
  return (1);
}

/*
 * trasformazione topologica del tipo
 * 
 *      \   /        \     /
 *    b1n\ /          \   /
 *        X  ===>      \ /
 *       / \      bleft/ \bright
 *      /   \         /   \
 *     /     \       /     \
 *
 * assumiamo una orientazione consistente
 * degli archi coinvolti
 */

void
taglia_nodo (struct border *b1n, struct sketch *sketch,
             struct border **bleftpt, struct border **brightpt)
{
  struct border *b1p, *b2p, *b3p, *b4p, *b2n;
  struct region *r2, *r4;
  struct arc *arc12, *arc23, *arc34, *arc41, *arcleft, *arcright;
  int changes, allowed_changes;
  extern int heisemberg;

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

  /* Se le regioni 2 e 4 sono la stessa regione
   * ci sono tre casi, ma li tratto tutti allo stesso modo
   * aggiustando alla fine con una chiamata a adjust_isexternalinfo
   * ad esempio se (r2->border != b2p->border)
   * significa che sono su un'isola, che si spezza in due
   */
  topo_change_g (b2p, b4p, TC_DILATION, sketch);
  changes = adjust_isexternalinfo (sketch);
  allowed_changes = 1;
  if (heisemberg >= 0) allowed_changes++;
  /* heisemberg procedure can mix connected components */
  if (changes > allowed_changes) fprintf (stderr, "changes: %d\n", changes);
  assert (changes <= allowed_changes);
  if (debug) printf ("changes = %d\n", changes);
  /* purtroppo non e' facile sapere quale delle due componenti
   * connesse sara' un buco, tiro a indovinare e poi utilizzo
   * la funzione adjust_isexternalinfo per sistemare
   */

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
    b1p = removeborder (b1p->next);
  }
  if (b2p->next != b2p)
  {
    if (b2p->next == b4p) b4p = b2p;
    b2p = removeborder (b2p->next);
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

  if (b3p->next != b3p) b3p = removeborder (b3p->next);
  if (b4p->next != b4p) b4p = removeborder (b4p->next);
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

  if (bleftpt) *bleftpt = b1p;
  if (brightpt) *brightpt = b3p;
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
    assert (a[i]->cusps == 0);
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
 * unisci una coppia di archi in un nuovo nodo.
 * viene chiamata per unire due cuspidi (rule C2)
 * (con "removecusps = 1")
 * e come fase intermedia per l'inverso delle regole
 * Nx (con "removecusps = 0") [era "join_cusps"]
 *
 * NOTA IMPORTANTE: alla fine, i valori cusp?pos non hanno
 * piu' alcun significato (sostanzialmente il nuovo valore e' 0)
 * i puntatori bp1 e bp2 POSSONO richiedere una modifica,
 * in particolare quando in ingresso sono uguali!
 *
 * per ora e' fatta solo per il caso di orientazione negativa
 */

void
pinch_arcs (struct border **bp1pt, int cusp1pos, 
            struct border **bp2pt, int cusp2pos, 
            struct sketch *sketch, int removecusps)
{
  struct borderlist *bl1, *bl2;
  struct border *bp1, *bp2;
  int ori;

  assert (removecusps == 0 || removecusps == 1);
  bp1 = *bp1pt; bp2 = *bp2pt;
  spezza_bordo (bp1, cusp1pos, sketch, removecusps);
  if (debug) printf ("dopo il primo spezza_bordo\n");
  if (debug) printsketch (sketch);
  if (bp1 == bp2)
  {
    /* devo ricalcolare bp2 e cusp2pos tenendo conto
     * di quanto fatto da 'spezza_bordo'
     */
    ori = bp1->orientation;
    if (bp1->info->endpoints == 1)
    {
      if (debug) printf ("era un S1\n");
      if (cusp2pos > cusp1pos) cusp2pos -= cusp1pos + removecusps;
        else cusp2pos += bp1->info->cusps - cusp1pos;
    } else {
      if (cusp2pos > cusp1pos)
      {
        cusp2pos -= cusp1pos + removecusps;
        if (ori > 0) bp2 = bp2->next;   /* cambiera' bp1 se ori < 0 */
      } else {
        if (ori < 0) bp2 = bp2->next;   /* cambiera' bp1 se ori > 0 */
      }
    }
  }
  spezza_bordo (bp2, cusp2pos, sketch, removecusps);
  if (bp1 == bp2) bp1 = bp2->next;
  if (debug) printf ("dopo il secondo spezza_bordo\n");
  if (debug) printsketch (sketch);

  bl1 = bp1->border;
  bl2 = bp2->border;
  if (bl1 == bl2)         /* la regione si spezza in due */
  {
    if (debug) printf ("la regione si spezza in due\n");

    topo_change_g (bp1, bp2, TC_EROSION, sketch);
    /* attenzione, si puo' perdere l'informazione isexternalinfo */

    if (bp1->next == bp1) bp1->info->endpoints = 1;
    if (bp2->next == bp2) bp2->info->endpoints = 1;
  } else {                /* rimane una sola regione */
    if (bl2 != bl2->region->border)   /* e' un buco */
    {
      redefineborder (bp2, bl1);
      bl1 = bl2;
    } else {
      redefineborder (bp1, bl2);
    }
    topo_change_l (bp1, bp2);
    if (bp1->next == bp1) bp1->info->endpoints = 1;
    if (bp2->next == bp2) bp2->info->endpoints = 1;
    extractborderlist (bl1);
    bl1->sponda = 0;
    freeborderlist (bl1);
  }
  *bp1pt = bp1; *bp2pt = bp2;
}

/*
 * suddivido un arco in corrispondenza di una cuspide,
 * (se "removecusp" = 1)
 * eventualmente creo i border necessari
 * se removecusp = 0 l'arco viene tagliato tra due cuspidi
 */

void
spezza_bordo (struct border *bp, int cusppos, struct sketch *sketch,
		int removecusp)
{
  int *newdepths, i, j;
  struct arc *arc, *arcnew;
  struct border *bt, *btnew;
  struct border *bpnew;

  assert (removecusp == 0 || removecusp == 1);
  arc = bp->info;
  if (removecusp == 1) assert (bp->orientation < 0);
  if (arc->endpoints == 0)  /* questo e' un s1 */
  {
    /* in questo caso non si creano altri archi e bordi.
     * l'apertura di un S1 aumenta di uno il numero di
     * valori di d, anche se il dimensionamento rimane
     * uguale perche' per gli S1 i valori estremi si
     * ripetono uguali.
     * se si rimuove una cuspide il numero di valori di 
     * d rimane uguale, ma non serve ripetere il primo valore
     * TODO: Mmm, probabilmente non conviene riallocare,
     * ma si puo' continuare ad usare lo spazio vecchio
     */
    arc->cusps -= removecusp;
    arc->depthsdim -= removecusp;
    arc->dvalues += 1 - removecusp;
    //assert (arc->dvalues <= arc->depthsdim);  // dvalues non affidabile 
    newdepths = (int *) malloc (arc->depthsdim * sizeof (int));
    for (i = 0; i < arc->depthsdim; i++)
    {
      j = i + cusppos + removecusp;
      if (j >= arc->depthsdim) j -= arc->depthsdim - 1 + removecusp;
      newdepths[i] = arc->depths[j];
    }
    free (arc->depths);
    arc->depths = newdepths;
    arc->endpoints = 1;
  } else {
    /* in questo caso si crea un nuovo arco e due nuove sponde */
    arcnew = newarc (sketch);
    bt = gettransborder (bp);
    btnew = newborder (bt->border);
    btnew->next = bt->next;
    bt->next = btnew;
    btnew->orientation = bt->orientation;
    bpnew = newborder (bp->border);
    bpnew->next = bp->next;
    bp->next = bpnew;
    bpnew->orientation = bp->orientation;

    if (bp->orientation < 0) {
      bp->info = arcnew;
      btnew->info = arcnew;
      arcnew->regionleft = btnew;
      arcnew->regionright = bp;
      bpnew->info = arc;
      arc->regionright = bpnew;
    } else {
      bt->info = arcnew;
      btnew->info = arc;
      arcnew->regionleft = bpnew;
      arcnew->regionright = bt;
      bpnew->info = arcnew;
      arc->regionright = btnew;
    }

    arcnew->cusps = arc->cusps - cusppos - removecusp;
    arcnew->depthsdim = arcnew->cusps + 1;
    arcnew->dvalues = arcnew->depthsdim;
    newdepths = (int *) malloc (arcnew->depthsdim * sizeof (int));
    arcnew->depths = newdepths;
    /* devo definire endpoints, che ora e' 2 */
    arc->endpoints = arcnew->endpoints = 2;
    for (i = 0; i <= arcnew->cusps; i++)
    {
      newdepths[i] = arc->depths[i + cusppos + removecusp];
    }
    arc->cusps = cusppos;
    arc->dvalues = cusppos + 1;
    newdepths = (int *) malloc ((cusppos + 1) * sizeof (int));
    for (i = 0; i < arc->cusps + 1; i++)
    {
      newdepths[i] = arc->depths[i];
    }
    free (arc->depths);
    arc->depths = newdepths;
    arc->depthsdim = cusppos + 1;
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
 * oppure senza cuspidi (regola cr2)
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
  taglia_nodo (b1n, sketch, 0, 0);
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
  newdepths = (int *) malloc ((arc->cusps + 2)*sizeof (int));
  if (bsegtn->orientation > 0)
  {
    newdepths[0] = arc->depths[0] - arccusp->depths[1]
                                  + arccusp->depths[0];
    for (i = 0; i <= arc->cusps; i++)
      newdepths[i+1] = arc->depths[i];
  } else {
    for (i = 0; i <= arc->cusps; i++)
      newdepths[i] = arc->depths[i];

    newdepths[arc->cusps + 1] = newdepths[arc->cusps]
       - arccusp->depths[0] + arccusp->depths[1];
  }
  arc->cusps++;
  arc->dvalues++;
  free (arc->depths);
  arc->depths = newdepths;
  arc->depthsdim = arc->cusps + 1;
  /* devo togliere la cuspide da arccusp, non sto a liberare/allocare */
  arccusp->cusps = 0;
  // arccusp->depthsdim = 1;
  arccusp->dvalues = 1;
  if (bsegtn->orientation < 0)
    arccusp->depths[0] = arccusp->depths[1];

  if (debug) checkconsistency (sketch);

  outnw = rimuovi_losanga (r, sketch);
  taglia_nodo (outnw, sketch, 0, 0);
}

/*
 * applicazione mossa CR3 = C2^{-1} + CN2R
 * rimozione di un nodo a favore di due cuspidi;
 * la struttura e' piuttosto complicata...
 */

void
applyrulecr3 (struct border *b1, struct border *b2, int dindex,
              struct sketch *sketch)
{
  struct border *b1nt, *b1ntnt, *b2t, *b1t, *b1tp, *b1n;
  struct arc *arc1, *arc1n, *arc2, *arc3, *arc3n;
  int isclosed1, isclosed2, size1, size2, size12, i, j, count;
  int *newdepths1, *newdepths2, *newdepths12;

  /* devo stare attento a non perdere i dati di b2 e dindex */

  if (debug) printsketch (sketch);
  b1n = b1->next;
  b1nt = gettransborder (b1n);
  b1ntnt = gettransborder (b1nt->next);
  b2t = gettransborder (b2);
  b1t = gettransborder (b1);
  b1tp = prevborder (b1t);

  arc1 = b1->info;
  arc1n = b1->next->info;
  arc2 = b2->info;
  arc3 = b1tp->info;
  arc3n = b1nt->next->info;
  if (debug) printf ("arc1 %d, arc1n %d, arc2 %d, arc3 %d, arc3n %d\n",
              arc1->tag, arc1n->tag, arc2->tag, arc3->tag, arc3n->tag);
  arc3 = mergearcs (arc3, arc3n, sketch);
  arc3->regionleft = b1tp;
  arc3->regionright = b1ntnt;
  b1tp->info = b1ntnt->info = arc3;

  if (arc2->endpoints == 0)
  {
    size12 = arc1->cusps + arc2->cusps + 3;
    if (arc1 != arc1n) size12 += arc1n->cusps;
    newdepths12 = (int *) malloc (size12 * sizeof (int));
    for (i = 0, j = 0; i <= arc1->cusps; i++)
      newdepths12[j++] = arc1->depths[i];
    for (i = dindex; i <= arc2->cusps; i++)
      newdepths12[j++] = arc2->depths[i];
    for (i = 1; i <= dindex; i++)
      newdepths12[j++] = arc2->depths[i];
    newdepths12[j++] = arc1n->depths[0];
    if (arc1 != arc1n)
      for (i = 1; i <= arc1n->cusps; i++)
        newdepths12[j++] = arc1n->depths[i];
    free (arc1->depths);
    arc1->depths = 0;
    if (arc1n->depths) free (arc1n->depths);
    arc1n->depths = 0;
    if (arc2->depths) free (arc2->depths);
    arc2->depths = 0;
    arc1->depths = newdepths12;
    arc1->depthsdim = arc1->dvalues = arc1->cusps = size12;
    arc1->cusps--;
    arc1->endpoints = 1;
    if (b1->next->next != b1) arc1->endpoints = 2;
    if (arc1 == arc1n) {arc1->dvalues--; arc1->endpoints = 0;}
    if (arc1n->depths == 0)
    {
      arc1n->regionleft = arc1n->regionright = 0;
      removearc (arc1n, sketch);
    }
    arc1n = 0;
  } else if (arc1 == arc1n)
  {
    size12 = arc1->cusps + arc2->cusps + 3;
    newdepths12 = (int *) malloc (size12 * sizeof (int));
    for (i = 0, j = 0; i <= dindex; i++)
      newdepths12[j++] = arc2->depths[i];
    for (i = 0; i <= arc1->cusps; i++)
      newdepths12[j++] = arc1->depths[i];
    for (i = dindex; i <= arc2->cusps; i++)
      newdepths12[j++] = arc2->depths[i];
    free (arc1->depths);
    arc1->depths = 0;
    if (arc1n->depths) free (arc1n->depths);
    arc1n->depths = 0;
    if (arc2->depths) free (arc2->depths);
    arc2->depths = 0;
    arc1->depths = newdepths12;
    arc1->depthsdim = arc1->dvalues = arc1->cusps = size12;
    arc1->cusps--;
    arc1->endpoints = 2;
    if (b2t->next == b2t) arc1->endpoints = 1;
    arc1n = 0;
  } else {
    assert (arc1 != arc1n);
    isclosed1 = isclosed2 = 0;
    if (arc1 != arc2)
    {
      size1 = arc2->cusps - dindex + arc1->cusps + 2;
      newdepths1 = (int *) malloc (size1 * sizeof (int));
      for (i = 0, j = 0; i <= arc1->cusps; i++) 
        newdepths1[j++] = arc1->depths[i];
      assert (abs(arc2->depths[dindex] - newdepths1[j-1]) == 1);
      for (i = dindex; i <= arc2->cusps; i++)
        newdepths1[j++] = arc2->depths[i];
    } else {
      isclosed1 = 1;
      size1 = arc1->cusps - dindex + 2;
      newdepths1 = (int *) malloc (size1 * sizeof (int));
      for (i = dindex, j = 0; i <= arc1->cusps; i++)
        newdepths1[j++] = arc1->depths[i];
      newdepths1[j] = arc1->depths[dindex];
    }
    if (arc1n != arc2)
    {
      size2 = dindex + arc1n->cusps + 2;
      newdepths2 = (int *) malloc (size2 * sizeof (int));
      for (i = 0, j = 0; i <= dindex; i++)
        newdepths2[j++] = arc2->depths[i];
      for (i = 0; i <= arc1n->cusps; i++)
        newdepths2[j++] = arc1n->depths[i];
      assert (abs(newdepths2[dindex + 1] - newdepths2[dindex]) == 1);
    } else {
      isclosed2 = 1;
      size2 = dindex + 2;
      newdepths2 = (int *) malloc (size2 * sizeof (int));
      for (i = 0, j = 0; i <= dindex; i++)
        newdepths2[j++] = arc1n->depths[i];
      newdepths2[j] = arc1n->depths[0];
    }
    free (arc1->depths);
    arc1->depths = 0;
    if (arc1n->depths) free (arc1n->depths);
    arc1n->depths = 0;
    if (arc2->depths) free (arc2->depths);
    arc2->depths = 0;
    arc1->depths = newdepths1;
    arc1->depthsdim = arc1->dvalues = arc1->cusps = size1;
    arc1->cusps--;
    arc1->endpoints = 1;
    if (b1->next != b1) arc1->endpoints = 2;
    if (isclosed1) {arc1->dvalues--; arc1->endpoints = 0;}
    if (arc1n->depths != 0) arc1n = newarc (sketch);
    arc1n->depths = newdepths2;
    arc1n->depthsdim = arc1n->dvalues = arc1n->cusps = size2;
    arc1n->cusps--;
    arc1n->endpoints = 1;
    if (b1n->next != b2) arc1n->endpoints = 2;
    if (isclosed2) {arc1n->dvalues--; arc1n->endpoints = 0;}
  }
  if (b1 == b2)
  {
    b1 = newborder (b2->border);
    b1->orientation = b2->orientation;
    b1->next = b2->next;
    b2->next = b1;
  }
  arc1->regionleft = b1;
  arc1->regionright = b2t;
  b1->info = b2t->info = arc1;
  if (b1nt == b2t)
  {
    b1nt = newborder (b2t->border);
    b1nt->orientation = b2t->orientation;
    b1nt->next = b2t->next;
    b2t->next = b1nt;
  }
  if (arc1n)
  {
    arc1n->regionleft = b2;
    arc1n->regionright = b1nt;
    b2->info = b1nt->info = arc1n;
  }
  if (arc2->depths == 0)
  {
    arc2->regionleft = arc2->regionright = 0;
    removearc (arc2, sketch);
  }
  if (b1t != b2t) b1t->info = 0; /* perche' ora e' b2t che punta all'arco */
  assert (b1 != b2);
  topo_change_g (b1, b2, TC_EROSION, sketch);
  assert (b1tp != b1nt);
  topo_change_g (b1tp, b1nt, TC_DILATION, sketch);
  assert (b1nt != b2t);
  topo_change_g (b1nt, b2t, TC_DILATION, sketch);

  /* rimozioni */
  if (b2 != b2->next) {b2->next->info = 0; removeborder (b2->next);}
  if (b2t != b2t->next) removeborder (b2t->next);
  if (b1tp != b1tp->next) {b1tp->next->info = 0; removeborder (b1tp->next);}
  if (b1ntnt != b1ntnt->next) removeborder (b1ntnt->next);
  if (arc1n == 0)
  {
    b2->info = 0; removeborder (b2);
    b1nt->info = 0; removeborder (b1nt);
  }
  if (debug) checkconsistency (sketch);
  if (debug) printsketch (sketch);
  count = adjust_isexternalinfo (sketch);
  if (debug) printf ("isexternalinfo, count = %d\n", count);
}
