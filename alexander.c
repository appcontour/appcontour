/*
 * computation of Alexander polinomial of
 * a finitely presented group with infinite ciclic
 * commutator quotient
 */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "contour.h"
#include "fundamental.h"
#include "laurent.h"
#include "alexander.h"
#include "fox.h"
#include "parser.h"
#include "groebner.h"

static int int_overflow_encountered = 0;

/* local prototypes */
void print_matrix1_n_m (struct laurentpoly ***matrix, int nrows, int ncols);
void print_matrix_n_m (struct laurentpoly ***matrix, int nrows, int ncols, int indets);
int laurentpoly_linf (struct laurentpoly *l);

int
alexander (struct presentation *p)
{
  extern int verbose, quiet, foxd, shuffle, outformat;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  struct laurentpoly *determinant2;
  struct laurentpoly **extradeterminants;
  struct alexanderideal *ai;
  int j, rank, matrixrank, deficiency;
  int numcols = 0;
  int printdvalue = 0;
  int foxdtoolarge = 0;
  int extradets, gconj, gconj2;

  topreabelian (p);

  if (verbose) print_presentation (p);

  for (r = p->rules; r; r = r->next) numcols++;
  /* foxd is the index of the Alexander ideal (dth Alexander ideal) */
  deficiency = p->gennum - numcols;
  if (deficiency < p->espected_deficiency && numcols > 0)
  {
    start_comment ();
    printf ("WARNING: deficiency: %d, espected value: %d.", deficiency, p->espected_deficiency);
    if (verbose)
    {
      printf ("\n\nYou can increase the deficiency by applying suitable punchhole surgeries that do not\n");
      printf ("affect the fundamental group.  The --autosurgery option automates this process.\n");
      printf ("Alternatively, use \"contour suggest_p_surgery\" on the apparent contour to obtain the correct surgery.\n");
      printf ("Prefix with \"contour wrap\" filter if you are computing the fundamental group of the outside\n\n");
    } else printf ("  Use --verbose option to get more details.\n");
  }
  if (foxd != FOXD_UNDEF && outformat == OUTFORMAT_APPCONTOUR) printf ("#\n# --foxd %d\n#\n", foxd);
  if (foxd < 0 && foxd >= FOXD_MINVALID)
    return (printout_constant_ideal ("Trivial zero ideal for negative d:\n", 0));
  if (foxd == 0)
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {}\n");
    printf ("Computation of order ideal not implemented\n");
    return (0);
  }
  if (foxd == FOXD_MAXINTERESTING)
  {
    printdvalue = 1;
    foxd = numcols;
    if (foxd < 1) foxd = 1;
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("#\n# --foxd %d (generators: %d)\n#\n", foxd, p->gennum);
  }
  if (foxd == FOXD_UNDEF)
  {
    printdvalue = 1;
    foxd = deficiency;
    if (foxd < 1) foxd = 1;
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("#\n# --foxd %d\n#\n", foxd);
  }
  if (printdvalue && outformat != OUTFORMAT_APPCONTOUR)
  {
    if (quiet) printf ("d = %d\n", foxd);
    else {
      printf ("Computing ");
      switch (foxd)
      {
        case 1: printf ("first"); break;
        case 2: printf ("second"); break;
        case 3: printf ("third"); break;
        default: printf ("%d-th", foxd); break;
      }
      printf (" Alexander ideal (--foxd = %d)\n", foxd);
    }
  }
  if (foxd >= p->gennum)
    return (printout_constant_ideal ("Trivial whole ring ideal:\n", 1));
  if (foxd < deficiency)
    /* Jacobian matrix has too low rank */
    return (printout_constant_ideal ("Alexander polynomial (special large deficiency case):\n", 0));
  rank = compute_fg_rank (p);
  if (rank < 0)
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {\n");
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
    printf ("Cannot compute Alexander polynomial for groups with torsion\n");
    return (0);
  }
  matrixrank = p->gennum - rank;  /* this is probably wrong for rank >= 3! */
  if (rank > 0) gconj = p->gennum;
  if (rank > 26)
  {
    assert (matrixrank > 0 || numcols > 0);
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {}\n");
    printf ("Cannot compute Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  if (rank < 1)
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {}\n");
    printf ("Cannot compute Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  gconj = p->gennum;
  if (rank >= 2) gconj2 = p->gennum - 1;

  if (verbose && rank < 3) printf ("Matrix a-la-Magnus-Karrass-Solitar has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols && deficiency != 1 && rank < 3)  /* deficiency == 1: link */
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {}\n");
    printf ("Cannot compute Alexander polynomial because there are too many relators\n");
    return (0);
  }

  if (matrixrank > numcols && rank < 3)
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {}\n");
    printf ("Cannot compute Alexander polynomial because there are too few relators\n");
    return (0);
  }

  switch (rank)
  {
    case 1:  /* case of knot groups */
    switch (foxd)
    {
      case 1:  /* this is the most interesting */
      determinant = laurent_eliminate_one_indeterminate (p, gconj);
      laurent_canonify1 (determinant);
      printout_ideal1 (0, determinant);
      break;

      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      ai = laurent_notfirst_elementary_ideal (p, gconj, foxd - 1);
      if (ai == 0) { foxdtoolarge++; break; }
      if (ai->l1num > 1)
      {
        start_comment ();
        printf ("*** Warning: result can be noncanonical ***\n");
      }
      alexander_fromideal (ai);
      break;

      default:
      foxdtoolarge++;
      break;
    }
    break;

    case 2:  /* two-component links or surface of genus 2 */
    assert (foxd >= deficiency);
    switch (foxd - deficiency)
    {
      case 0: /* 2 components link: foxd = 1; genus 2 surface: foxd = 2 */
      determinant2 = laurent_eliminate_two_indeterminates (p, gconj2, gconj, &extradeterminants);
      extradets = 1;
      if (deficiency == 2) extradets = numcols;
      if (extradeterminants == 0) extradets = 0;
      if (verbose)
      {
        printf ("Alexander ideal before canonization:\n");
        printf ("-------------\n");
        printf ("Main polynomial:\n");
        print_laurentpoly (determinant2, "uv");
        printf ("\n");
        for (j = 0; j < extradets; j++)
        {
          printf ("Fundamental ideal factor %d:\n", j);
          print_laurentpoly (extradeterminants[j], "uv");
          printf ("\n");
        }
        printf ("-------------\n");
      }
      if (shuffle)
      {
        shuffle_poly2 (&determinant2, extradeterminants, extradets);
      }
      if (canonify_ideal2 (&determinant2, extradeterminants, extradets) == 0)
      {
        start_comment ();
        printf ("*** Warning: result can be noncanonical ***\n");
      }
      if (extradeterminants)
      {
        if (deficiency == 2) laurent_canonify_exponents (determinant2);
        for (j = 0; j < extradets; j++) laurent_canonify_exponents (extradeterminants[j]);
      } else laurent_canonify_exponents (determinant2);
      printout_ideal (0, determinant2, extradeterminants, extradets, (deficiency == 2));
      break;

      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      if (deficiency > 1) { foxdtoolarge++; break; }
      /* case of link and request of second elementary ideal */
      ai = laurent_notfirst_elementary_ideal2 (p, gconj2, gconj, foxd - 1);
      if (ai == 0) { foxdtoolarge++; break; }
      /*
       * test used to be: 
       *    if (ai->l2num + ai->fl2num > 1) ...
       * but if ai->l2num == 1 the message would be printed twice
       */
      if (ai->l2num > 1)
      {
        start_comment ();
        printf ("*** Warning: result can be noncanonical ***\n");
      }
      alexander_fromideal (ai);
      break;

      default:
      foxdtoolarge++;

      break;
    }
    break;

    case 3:
    assert (foxd >= 1);
    switch (foxd)
    {
      case 1:
      assert (deficiency == 1);
      ai = three_components_link (p);
      if (ai == 0) { foxdtoolarge++; break; }
      //if (ai->l2num + ai->fl2num > 1) printf ("# *** Warning: result can be noncanonical ***\n");
      start_comment ();
      printf ("*** Warning: result can be noncanonical ***\n");
      alexander_fromideal (ai);
      break;

      default:
      ai = generic_ideal_computation (p, rank, p->gennum - foxd);
      if (ai == 0) return (0);
      start_comment ();
      printf ("*** Warning: result can be noncanonical ***\n");
      alexander_fromideal (ai);
      break;
    }
    break;

    default:
    ai = generic_ideal_computation (p, rank, p->gennum - foxd);
    if (ai == 0) return (0);
    start_comment ();
    printf ("*** Warning: result can be noncanonical ***\n");
    alexander_fromideal (ai);
    break;
  }

  if (foxdtoolarge)
  {
    printf ("Unable to compute higer d = %d Alexander ideal.\n", foxd);
    return (0);
  }
  return (1);
}

/*
 * We start from an ideal...
 */

int
alexander_fromideal (struct alexanderideal *ai)
{
  int i, extradets;
  extern int verbose, quiet;
  struct laurentpoly *determinant2, **extradeterminants;

  if (ai == 0)
  {
    if (!quiet) printf ("Error: not able to compute the required ideal.\n");
    return (0);
  }
  switch (ai->indets)
  {
    case 0:
      ai->val = abs (ai->val);
      printout_constant_ideal (0, ai->val);
      break;

    case 1:
      if (verbose)
      {
        printf ("Alexander ideal before canonization:\n{\n");
	printout_ideal1 (ai, 0);
      }
      if (ai->l1num == 1) laurent_canonify1 (ai->l[0]);
       else {
        for (i = 0; i < ai->l1num; i++)
        {
          ai->l[i]->minexpon = 0;
          laurent_canonifysign1 (ai->l[i]);
        }
      }
      printout_ideal1 (ai, 0);
      break;

    case 2:
      if (verbose)
      {
        printf ("Alexander ideal before canonization:\n{\n");
        for (i = 0; i < ai->l2num; i++)
        {
          print_laurentpoly (ai->l[i], "uv");
          printf (";\n");
        }
        for (i = 0; i < ai->fl2num; i++)
        {
          printf ("F: ");
          print_laurentpoly (ai->l[i + ai->fl2offset], "uv");
          printf (";\n");
        }
        printf ("}\n");
      }
      if (ai->l2num > 1) /* cannot use canonify_ideal2 */
      {
        for (i = 0; i < ai->l2num; i++) laurent_canonify_exponents (ai->l[i]);
        for (i = 0; i < ai->fl2num; i++) laurent_canonify_exponents (ai->l[i + ai->fl2offset]);
      } else {
        determinant2 = 0;
        assert (ai->l2num <= 1);
        if (ai->l2num == 1) determinant2 = ai->l[0];
        extradets = ai->fl2num;
        extradeterminants = 0;
        if (extradets > 0)
        {
          extradeterminants = (struct laurentpoly **) malloc (extradets*sizeof (struct laurentpoly *));
          for (i = 0; i < extradets; i++) extradeterminants[i] = ai->l[i + ai->fl2offset];
        }
        if (canonify_ideal2 (&determinant2, extradeterminants, extradets) == 0)
        {
          start_comment ();
          printf ("*** Warning: result can be noncanonical ***\n");
          for (i = 0; i < ai->l2num; i++) laurent_canonify_exponents (ai->l[i]);
          for (i = 0; i < ai->fl2num; i++) laurent_canonify_exponents (ai->l[i + ai->fl2offset]);
        }
        if (ai->l2num == 1) ai->l[0] = determinant2;
        if (extradets > 0)
        {
          for (i = 0; i < extradets; i++) ai->l[i + ai->fl2offset] = extradeterminants[i];
          free (extradeterminants);
        }
      }
      printout_ideal (ai, 0, 0, 0, 0);
      break;

    default:
      printout_ideal (ai, 0, 0, 0, 0);
      break;
  }
  return (1);
}

/*
 * printout a one-indet ideal
 */

void
printout_ideal1 (struct alexanderideal *ai, struct laurentpoly *principal) 
{
  extern int outformat, quiet;
  int j;

  if (ai) assert (ai->indets == 1);
  if (ai) assert (principal == 0);

  if (ai && ai->l1num > 1)
  {
    switch (outformat)
    {
      case OUTFORMAT_APPCONTOUR:
      printf ("ideal(t) {\n");
      break;

      case OUTFORMAT_MACAULAY2:
      start_comment ();
      printf ("M2 -- Macaulay2 input commands:\n");
      printf ("S = ZZ[t,u]\n");
      printf ("groebnerBasis ideal (t*u-1,\n");
      break;

      default:
      if (quiet) printf ("Ideal [\n");
       else printf ("Alexander ideal generated by:\n");
      break;
    }
    switch (outformat)
    {
      case OUTFORMAT_MACAULAY2:
      if (principal) {printf ("  "); print_laurentpoly (principal, "t"); printf (")\n");}
      if (ai)
      {
        for (j = 0; j < ai->l1num; j++)
        {
          print_laurentpoly (ai->l[j], "t");
          printf ("%c\n", (j < ai->l1num-1)?',':')');
        }
      }
      break;

      case OUTFORMAT_APPCONTOUR:
      default:
      if (principal) {print_laurentpoly (principal, "t"); printf (";\n");}
      if (ai)
      {
        for (j = 0; j < ai->l1num; j++)
        {
          print_laurentpoly (ai->l[j], "t");
          printf (";\n");
        }
      }
      break;
    }
    switch (outformat)
    {
      case OUTFORMAT_APPCONTOUR:
      case OUTFORMAT_MACAULAY2:
      break;

      default:
      if (quiet) printf ("]\n");
      break;
    }
  } else {
    if (ai) assert (ai->spread >= 1);
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander(t) {\n");
    if (quiet)
    {
      if (ai && ai->spread > 1) printf ("spread = %d\n", ai->spread);
    } else {
      if (ai && ai->spread > 1) printf ("Fuzzy ");
      printf ("Alexander polynomial (up to t -> 1/t)");
      if (ai && ai->spread > 1) printf (" with spread factor %d", ai->spread);
      printf (":\n");
    }
    if (ai == 0) print_laurentpoly (principal, "t");
     else {
      assert (ai->l1num == 1);
      print_laurentpoly (ai->l[0], "t");
    }
    printf (";\n");
  }
  if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
}

/*
 * printout a two-indets ideal
 */

void
printout_ideal (struct alexanderideal *ai, struct laurentpoly *principal, 
                     struct laurentpoly **fundamentals, int fnum, int printprincipal)
{
  extern int outformat, quiet;
  char extraindets[]=",w,x,y,z,...";
  char extraindets2[]=",w,ww,x,xx,y,yy,z,zz,...";
  char extraindets3[]=",w*ww-1,x*xx-1,y*yy-1,z*zz-1,...";
  int j;

  if (fundamentals == 0) assert (fnum == 0);
  if (ai)
  {
    assert (ai->indets >= 2);
    if (ai->indets <= 6)
    {
      extraindets[2*(ai->indets - 2)] = 0;
      extraindets2[5*(ai->indets - 2)] = 0;
      extraindets3[7*(ai->indets - 2)] = 0;
    }
    assert (principal == 0 && fundamentals == 0);
  } else extraindets[0] = 0;

  if (fnum > 0 || (ai && ai->fl2num > 0) || (ai && ai->l2num >= 2))
  {
    switch (outformat)
    {
      case OUTFORMAT_APPCONTOUR:
      printf ("ideal(u,v%s) {\n", extraindets);
      break;

      case OUTFORMAT_MACAULAY2:
      start_comment ();
      printf ("M2 -- Macaulay2 input commands:\n");
      printf ("S = ZZ[u,uu,v,vv%s]\n", extraindets2);
      printf ("groebnerBasis ideal (u*uu-1,v*vv-1%s,\n", extraindets3);
      break;

      default:
      if (quiet) printf ("Ideal [\n");
       else printf ("Alexander ideal generated by:\n");
      break;
    }

    if (printprincipal || principal || (ai && ai->l2num > 0))
    {
      if (principal || printprincipal)
      {
        print_laurentpoly (principal, "uvwxyz");
        printf ("%c\n", (outformat == OUTFORMAT_MACAULAY2)?',':';');
      }
      if (ai)
      {
        for (j = 0; j < ai->l2num; j++)
        {
          print_laurentpoly (ai->l[j], "uvwxyz");
          printf ("%c\n", (outformat == OUTFORMAT_MACAULAY2)?',':';');
        }
      }
    }
    for (j = 0; j < fnum; j++)
    {
      if (fundamentals[j] == 0) continue;
      switch (outformat)
      {
        case OUTFORMAT_APPCONTOUR:
        printf ("F: ");
        print_laurentpoly (fundamentals[j], "uvwxyz");
        printf (";\n");
        break;

        case OUTFORMAT_MACAULAY2:
        printf ("(");
        print_laurentpoly (fundamentals[j], "uvwxyz");
        printf (")*(u - 1),\n");
        printf ("(");
        print_laurentpoly (fundamentals[j], "uvwxyz");
        printf (")*(v - 1),\n");
        if (ai && ai->indets > 2) printf ("[...],\n");
        break;

        default:
        printf ("(");
        print_laurentpoly (fundamentals[j], "uvwxyz");
        printf (") (u - 1);\n");
        printf ("(");
        print_laurentpoly (fundamentals[j], "uvwxyz");
        printf (") (v - 1);\n");
        if (ai && ai->indets > 2) printf ("[...];\n");
        break;
      }
    }
    if (ai)
    {
      for (j = 0; j < ai->fl2num; j++)
      {
        if (ai->l[j + ai->fl2offset] == 0) continue;
        switch (outformat)
        {
          case OUTFORMAT_APPCONTOUR:
          printf ("F: ");
          print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
          printf (";\n");
          break;

          case OUTFORMAT_MACAULAY2:
          printf ("(");
          print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
          printf (")*(u - 1),\n");
          printf ("(");
          print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
          printf (")*(v - 1),\n");
          if (ai->indets == 3)
          {
            printf ("(");
            print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
            printf (")*(w - 1),\n");
          }
          if (ai->indets > 3) printf ("[...],\n");
          break;

          default:
          printf ("(");
          print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
          printf (") (u - 1);\n");
          printf ("(");
          print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
          printf (") (v - 1);\n");
          if (ai->indets == 3)
          {
            printf ("(");
            print_laurentpoly (ai->l[j + ai->fl2offset], "uvwxyz");
            printf (") (w - 1);\n");
          }
          if (ai->indets > 3) printf ("[...];\n");
          break;
        }
      }
    }
    switch (outformat)
    {
      case OUTFORMAT_APPCONTOUR:
      break;

      case OUTFORMAT_MACAULAY2:
      printf ("0)\n");
      break;

      default:
      if (quiet) printf ("]\n");
    }
  } else {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander(u,v%s) {\n", extraindets);
    if (!quiet) printf ("Alexander polynomial:\n");
    if (ai == 0) print_laurentpoly (principal, "uvwxyz");
     else {
      //if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander(u,v%s) {\n", extraindets);
      assert (ai->l2num == 1);
      print_laurentpoly (ai->l[0], "uvwxyz");
    }
    printf (";\n");
  }
  if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
}

/*
 * printout of a constant ideal
 */

int
printout_constant_ideal (char *msg, int val)
{
  extern int quiet, outformat;

  if (!quiet && msg) printf ("%s", msg);
  if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {\n%d;\n}\n", val);
    else printf ("%d;\n", val);
  return (1);
}

/*
 * compute the linking number of a two-component link
 */

int
linkingnumber (struct presentation *p)
{
  extern int verbose, quiet;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  int i, rank, matrixrank;
  int deriv, val, e;
  int numcols = 0;
  int gconj;

  topreabelian (p);

  if (verbose) print_presentation (p);

  for (r = p->rules; r; r = r->next) numcols++;
  rank = compute_fg_rank (p);
  if (rank < 0)
  {
    printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
    return (0);
  }
  matrixrank = p->gennum - rank;
  //for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  //{
  //  sum = get_exp_sum1 (r, i);
  //  assert (sum >= 0);
  //  if (sum) matrixrank = i;
  //  if (sum && sum != 1)
  //  {
  //    printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
  //    return (0);
  //  }
  //}
  //rank = p->gennum - matrixrank;
  if (rank <= 1)
  {
    printf ("Cannot compute linking number for groups with rank %d\n", rank);
    return (0);
  }
  if (rank > 2)
  {
    printf ("Cannot compute linking number for groups with rank %d\n", rank);
    return (0);
  }
  rank = 1;
  matrixrank = p->gennum - rank;
  gconj = p->gennum;

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols)
  {
    printf ("Cannot compute linking number because there are too many relators\n");
    return (0);
  }

  determinant = laurent_eliminate_one_indeterminate (p, gconj);

  //laurent_canonify1 (determinant);

  if (verbose)
  {
    printf ("Determinant of matrix: ");
    print_laurentpoly (determinant, "t");
    printf (";\n");
  }

  /* linking number is the derivative in t=1 */

  deriv = 0;
  val = 0;
  if (determinant != 0)
  {
    for (i = 0; i <= determinant->stemdegree; i++)
    {
      e = i + determinant->minexpon;
      deriv += e*determinant->stem[i].l0;
      val += determinant->stem[i].l0;
    }
  }
  if (val != 0) printf ("Warning: Value in 1 is %d, should be zero.\n", val);
  deriv = abs(deriv);
  if (quiet) printf ("%d\n", deriv);
   else printf ("Linking number is %d\n", deriv);
  return (1);
}

/*
 * linking number computed from Alexander ideal
 */

int
linkingnumber_fromideal (struct alexanderideal *ai)
{
  printf ("Not implemented!\n");
  return (0);
}

/*
 * corank one (experimental)
 */

int
corank_one_alexander (struct presentation *p)
{
  extern int verbose, quiet;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  int i, rank, matrixrank;
  int deriv, val, e;
  int numcols = 0;
  int gconj;

  topreabelian (p);

  if (verbose) print_presentation (p);

  for (r = p->rules; r; r = r->next) numcols++;
  rank = compute_fg_rank (p);
  if (rank < 0)
  {
    printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
    return (0);
  }
  matrixrank = p->gennum - rank;
  //matrixrank = 0;
  //for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  //{
  //  sum = get_exp_sum1 (r, i);
  //  assert (sum >= 0);
  //  if (sum) matrixrank = i;
  //  if (sum && sum != 1)
  //  {
  //    printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
  //    return (0);
  //  }
  //}
  //rank = p->gennum - matrixrank;
  if (rank <= 1)
  {
    printf ("Cannot compute corank one Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  if (rank > 2)
  {
    printf ("Cannot compute corank one Alexander polynomial for groups with rank %d\n", rank);
    return (0);
  }
  rank = 1;
  matrixrank = p->gennum - rank;
  gconj = p->gennum;

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols)
  {
    printf ("Cannot compute Alexander polynomial because there are too many relators\n");
    return (0);
  }

  determinant = laurent_eliminate_one_indeterminate (p, gconj);

  laurent_canonify1 (determinant);

  if (verbose)
  {
    printf ("Determinant of matrix (before canonization): ");
    print_laurentpoly (determinant, "t");
    printf (";\n");
  }

  /* calcolo della derivata in t=1 */

  deriv = 0;
  val = 0;
  if (determinant != 0)
  {
    for (i = 0; i <= determinant->stemdegree; i++)
    {
      e = i + determinant->minexpon;
      deriv += e*determinant->stem[i].l0;
      val += determinant->stem[i].l0;
    }
  }
  printf ("Value in 1: %d, deriv = %d\n", val, deriv);
  if (!quiet) printf ("Alexander polynomial (up to t -> 1/t):\n");
  print_laurentpoly (determinant, "t");
  printf (";\n");
  return (1);
}

/*
 * eliminate one indeterminate by conjugacy
 */

struct laurentpoly *
laurent_eliminate_one_indeterminate (struct presentation *p, int eliminate)
{
  extern int verbose;
  struct laurentpoly *determinant;
  struct laurentmatrix *matrix;

  matrix = laurent_build_matrix1 (p, eliminate);

  assert (matrix->numcols <= matrix->numrows);
  if (matrix->numcols < matrix->numrows) return (0);

  determinant = laurent_compute_determinant1 (matrix->columns, matrix->numcols);

  laurent_free_matrix (matrix);

  return (determinant);
}

/*
 * second elementary ideal for knots (one indet)
 */

struct alexanderideal *
laurent_notfirst_elementary_ideal (struct presentation *p, int eliminate, int corank)
{
  extern int verbose, simplifyideal;
  struct alexanderideal *ai;
  struct laurentmatrix *matrix, *minor;
  struct laurentpoly *l, **matrixcolumn, **columnj, **columnjj;
  int rank, i, ii, j, jj;

  assert (corank >= 1);
  matrix = laurent_build_matrix1 (p, eliminate);

  assert (matrix->numcols <= matrix->numrows);
  if (matrix->numcols < matrix->numrows - corank) return (0); /* TODO: this should be the trivial "0" ideal, instead */

  rank = matrix->numrows - corank;
  ai = (struct alexanderideal *) malloc (AI_DIM(IDEAL_DEF_GENERATORS_NUM));
  ai->max_generators_num = IDEAL_DEF_GENERATORS_NUM;
  ai->indets = 1;
  ai->spread = 1;
  if (rank <= corank)
  {
    switch (rank)
    {
      case 0:
      ai->indets = 0;
      ai->l1num = 1;
      ai->val = 1;
      laurent_free_matrix (matrix);
      return (ai);
      break;

      case 1:
      ai->l1num = 0;
      for (i = 0; i < matrix->numrows; i++)
      {
        for (j = 0; j < matrix->numcols; j++)
        {
          matrixcolumn = matrix->columns[j];
          l = matrixcolumn[i];
          if (l == 0) continue;
          if (ai->l1num >= ai->max_generators_num) ai = ai_increase_size (ai);
          if (ai->l1num >= ai->max_generators_num)
          {
            printf ("Fatal: too many generators (%d) for the ideal\n", ai->l1num);
            laurent_free_matrix (matrix);
            free (ai);
            return (0);
          } else {
            // assert (idx < ai->l1num);
            ai->l[ai->l1num++] = laurent_dup(l);
          }
        }
      }
      break;

      case 2:
      ai->l1num = 0;
      minor = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
      minor->numcols = 2;
      minor->numrows = 2;
      minor->columns = (struct laurentpoly ***) malloc (2*sizeof(struct laurentpoly **));
      minor->columns[0] = (struct laurentpoly **) malloc (2*sizeof(struct laurentpoly *));
      minor->columns[1] = (struct laurentpoly **) malloc (2*sizeof(struct laurentpoly *));
      for (j = 0; j < matrix->numcols; j++)
      {
        columnj = matrix->columns[j];
        for (jj = j+1; jj < matrix->numcols; jj++)
        {
          columnjj = matrix->columns[jj];
          for (i = 0; i < matrix->numrows; i++)
          {
            for (ii = i+1; ii < matrix->numrows; ii++)
            {
              minor->columns[0][0] = columnj[i];
              minor->columns[0][1] = columnj[ii];
              minor->columns[1][0] = columnjj[i];
              minor->columns[1][1] = columnjj[ii];
              if (ai->l1num >= ai->max_generators_num) ai = ai_increase_size (ai);
              if (ai->l1num >= ai->max_generators_num)
              {
                printf ("Fatal: too many generators (%d) for the ideal\n", ai->l1num);
                free (minor->columns[0]);
                free (minor->columns[1]);
                free (minor->columns);
                free (minor);
                laurent_free_matrix (matrix);
                free (ai);
                return (0);
              }
              ai->l[ai->l1num++] = laurent_compute_determinant1 (minor->columns, minor->numcols);
              minor->columns[0][0] = 0;
              minor->columns[0][1] = 0;
              minor->columns[1][0] = 0;
              minor->columns[1][1] = 0;
            }
          }
        }
      }
      free (minor->columns[0]);
      free (minor->columns[1]);
      free (minor->columns);
      free (minor);
      break;

      default:
      printf ("Rank larger than 1 is not yet implemented\n");
      free (ai);
      laurent_free_matrix (matrix);
      return (0);
      break;
    }
  } else { /* corank less than rank */
    if (corank > 1)
    {
      printf ("Corank larger than 1 is not yet implemented\n");
      free (ai);
      laurent_free_matrix (matrix);
      return (0);
    }
    /* this is the corank 1 case, with rank larger than 1 */
    assert (matrix->numrows == matrix->numcols);
    ai->l1num = 0;
    for (i = 0; i < matrix->numrows; i++)
    {
      for (j = 0; j < matrix->numcols; j++)
      {
        if (ai->l1num >= ai->max_generators_num) ai = ai_increase_size (ai);
        if (ai->l1num >= ai->max_generators_num)
        {
          printf ("Fatal: too many generators (%d) for the ideal\n", ai->l1num);
          laurent_free_matrix (matrix);
          free (ai);
          return (0);
        }
        assert (ai->l1num < ai->max_generators_num);
        minor = minor_matrix1_corank1 (matrix, i, j);
        ai->l[ai->l1num++] = laurent_compute_determinant1 (minor->columns, minor->numcols);
        for (jj = 0; jj < minor->numcols; jj++) free (minor->columns[jj]);
        free (minor->columns);
        free (minor);
      }
    }
    if (verbose) printout_ideal1 (ai, 0);
  }

  laurent_free_matrix (matrix);
  if (simplifyideal) ai = laurent_simplify_ideal (ai);
  return (ai);
}

/*
 * second elementary ideal for two-components links (two indets)
 */

struct alexanderideal *
laurent_notfirst_elementary_ideal2 (struct presentation *p, int e1, int e2, int corank)
{
  extern int verbose;
  struct laurentmatrix *matrix, *minor;
  struct laurentpoly *l, **matrixcolumn;
  struct alexanderideal *ai;
  int i, j, jj, rank, lastrow;

  matrix = laurent_build_matrix2 (p, e1, e2);

  assert (matrix->numcols <= matrix->numrows);
  rank = matrix->numrows - corank;

  ai = (struct alexanderideal *) malloc (AI_DIM(IDEAL_DEF_GENERATORS_NUM));
  ai->max_generators_num = IDEAL_DEF_GENERATORS_NUM;
  ai->indets = 2;
  ai->fl2offset = ai->max_generators_num/2;
  ai->spread = 1;
  if (rank > matrix->numcols)
  {
    ai->indets = 0;
    ai->val = 0;
    return (ai);
  }

  if (rank <= corank)
  {
    switch (rank)
    {
      case 0:
      ai->indets = 0;
      ai->l2num = 1;
      ai->val = 1;
      laurent_free_matrix (matrix);
      return (ai);
      break;

      case 1:
      /* the main part of the ideal comes from a double loop: last row is excluded */
      ai->l2num = ai->fl2num = 0;
      for (j = 0; j < matrix->numcols; j++)
      {
        matrixcolumn = matrix->columns[j];
        for (i = 0; i < matrix->numrows - 1; i++)
        {
          l = matrixcolumn[i];
          if (l == 0) continue;
          if (ai->l2num >= ai->fl2offset) ai = ai_increase_size (ai);
          if (ai->l2num >= ai->fl2offset)
          {
            printf ("Fatal: too many generators in two indets (%d) for the ideal\n", ai->l2num);
            laurent_free_matrix (matrix);
            free (ai);
            return (0);
          }
          ai->l[ai->l2num++] = laurent_dup(l);
        }
      }
      /* the "fundamental" part of the ideal comes from a single loop: last row is included */
      lastrow = matrix->numrows - 1;
      for (j = 0; j < matrix->numcols; j++)
      {
        matrixcolumn = matrix->columns[j];
        l = matrixcolumn[lastrow];
        if (l == 0) continue;
        if (ai->fl2num + ai->fl2offset >= ai->max_generators_num) ai = ai_increase_size (ai);
        if (ai->fl2num + ai->fl2offset >= ai->max_generators_num)
        { 
          printf ("Fatal: too many f-generators in two indets (%d) for the ideal\n", ai->fl2num);
          laurent_free_matrix (matrix);
          free (ai);
          return (0);
        }
        ai->l[ai->fl2num++ + ai->fl2offset] = laurent_dup(l);
      }
      if (verbose) printout_ideal (ai, 0, 0, 0, 0);
      break;

      default:
      printf ("Rank larger than 0 is not yet implemented\n");
      free (ai);
      laurent_free_matrix (matrix);
      return (0);
      break;
    }
  } else { /* corank less than rank */
    if (corank > 1)
    {
      printf ("Corank larger than 1 is not yet implemented\n");
      free (ai);
      laurent_free_matrix (matrix);
      return (0);
    }
    assert (matrix->numrows == matrix->numcols);  /* the genus-2 case should be treated elsewhere */
    /* however it would be feasible to treat here the genus-2 case (foxd = 2) also */

    /* TODO:
     * jacobian matrix is (p+1)x(p+2) if p+2 is the number of generators
     * the last two columns are a multiple of (v-1) and (u-1) respectively
     * We need to compute determinants for two sets of pxp minors:
     * the (p+1) minors with the first p columns contribute to the main part of the ideal
     * the (p+1)x(p) minors with as last column the common factor of the last two columns
     * of the jacobian matrix: they contribute in the part that multiplies the fundamental
     * ideal
     *
     * in the special case p = 1 (three generators) we have 2 generators in the main part
     * and 2 generators in the fundamental part
     *
     * here the matrix is transposed and the last two columns of the jacobian is substituted
     * by one last row containing the common factors
     */

    /* this is the corank 1 case, with rank larger than 1 */
    /* the main part of the ideal comes from a single loop: last row is excluded */
    ai->l2num = ai->fl2num = 0;
    lastrow = matrix->numrows - 1;
    for (j = 0; j < matrix->numcols; j++)
    {
      if (ai->l2num >= ai->fl2offset) ai = ai_increase_size (ai);
      if (ai->l2num >= ai->fl2offset)
      {
        printf ("Fatal: too many generators in two indets (%d) for the ideal\n", ai->l2num);
        laurent_free_matrix (matrix);
        free (ai);
        return (0);
      }
      minor = minor_matrix2_corank1 (matrix, lastrow, j);
      ai->l[ai->l2num++] = laurent_compute_determinant (minor->columns, minor->numcols, 2);
      for (jj = 0; jj < minor->numcols; jj++) free (minor->columns[jj]);
      free (minor->columns);
      free (minor);
    }
    /* the "fundamental" part of the ideal comes from a couple of loops: last row is included */
    for (j = 0; j < matrix->numcols; j++)
    {
      for (i = 0; i < matrix->numrows - 1; i++)
      {
        if (ai->fl2num + ai->fl2offset >= ai->max_generators_num) ai = ai_increase_size (ai);
        if (ai->fl2num + ai->fl2offset >= ai->max_generators_num)
        {
          printf ("Fatal: too many f-generators in two indets (%d) for the ideal\n", ai->fl2num);
          laurent_free_matrix (matrix);
          free (ai);
          return (0);
        }
        minor = minor_matrix2_corank1 (matrix, i, j);
        ai->l[ai->fl2num++ + ai->fl2offset] = laurent_compute_determinant (minor->columns, minor->numcols, 2);
        for (jj = 0; jj < minor->numcols; jj++) free (minor->columns[jj]);
        free (minor->columns);
        free (minor);
      }
    }
    if (verbose) printout_ideal (ai, 0, 0, 0, 0);
  }

  laurent_free_matrix (matrix);
  //ai = laurent_simplify_ideal (ai);
  return (ai);
}

/*
 * build Alexander matrix (rank 1: e.g. knots)
 */

struct laurentmatrix *
laurent_build_matrix1 (struct presentation *p, int eliminate)
{
  struct presentationrule *r;
  struct laurentmatrix *matrix;
  struct laurentpoly **matrixcolumn;
  int numcols = 0;
  int i, ii, j, numrows;
  extern int verbose, debug;

  assert (eliminate >= 1);
  assert (eliminate <= p->gennum);

  for (r = p->rules; r; r = r->next) numcols++;
  numrows = p->gennum - 1;

  matrix = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
  matrix->numcols = numcols;
  matrix->numrows = numrows;

  matrix->columns = (struct laurentpoly ***) malloc (numcols*sizeof(struct laurentpoly **));
  for (j = 0; j < numcols; j++)
  {
    matrix->columns[j] = (struct laurentpoly **) malloc (numrows*sizeof(struct laurentpoly *));
  }

  /*
   * fill matrix
   */

  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix->columns[j];
    for (i = 1, ii = 0; i <= p->gennum; i++)
    {
      if (i == eliminate) continue;
      matrixcolumn[ii++] = laurent_get_exp_sum1 (r, i, eliminate);
    }
  }

  if (debug)
  {
    for (j = 0, r = p->rules; r; j++, r = r->next)
    {
      printf ("extra row, entry %d: ", j+1);
      print_laurentpoly (laurent_get_exp_sum1 (r, eliminate, eliminate), "t");
      printf (";\n");
    }
  }

  /*
   * stampa della matrice
   */

  if (verbose) print_matrix1 (matrix);
  return (matrix);
}

/*
 * build Alexander matrix (rank 2: e.g. links or genus-2 surfaces)
 */

struct laurentmatrix *
laurent_build_matrix2 (struct presentation *p, int e1, int e2)
{
  struct presentationrule *r;
  struct laurentmatrix *matrix;
  struct laurentpoly **matrixcolumn, *l;
  int numcols = 0;
  int i, ii, j, numrows;
  extern int verbose, debug;

  assert (e1 >= 1 && e2 >= 1);
  assert (e1 <= p->gennum);
  assert (e2 <= p->gennum);
  assert (e1 != e2);

  for (r = p->rules; r; r = r->next) numcols++;

  numrows = p->gennum - 1;

  matrix = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
  matrix->numcols = numcols;
  matrix->numrows = numrows;

  /*
   * last row will contain the common extra factor
   */

  matrix->columns = (struct laurentpoly ***) malloc (numcols*sizeof(struct laurentpoly **));
  for (j = 0; j < numcols; j++)
  {
    matrix->columns[j] = (struct laurentpoly **) malloc (numrows*sizeof(struct laurentpoly *));
  }

  /*
   * fill matrix
   */

  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix->columns[j];
    for (i = 1, ii = 0; i <= p->gennum; i++)
    {
      if (i == e1 || i == e2) continue;
      matrixcolumn[ii++] = laurent_get_exp_sum (r, i, e1, e2);
    }
    assert (ii == numrows - 1);
    matrixcolumn[ii] = laurent_common_factor (r, e1, e2);
  }

  if (debug)
  {
    for (j = 1, r = p->rules; r; j++, r = r->next)
    {
      printf ("Extra row 1, entry %d: ", j);
      print_laurentpoly (laurent_get_exp_sum (r, e1, e1, e2), "uv");
      printf (";\n");
      printf ("Extra row 2, entry %d: ", j);
      print_laurentpoly (laurent_get_exp_sum (r, e2, e1, e2), "uv");
      printf (";\n");
      printf ("Fox mixed derivative, entry %d: ", j);
      print_laurentpoly (laurent_mixed_derivative2 (r, e1, e2), "uv");
      printf (";\n");
      printf ("Common factor, entry %d: ", j);
      print_laurentpoly (matrix->columns[j-1][numrows-1], "uv");
      printf (";\n");
    }
  }

  /*
   * stampa della matrice
   */

  if (verbose)
  {
    printf ("Matrix entries:\n");
    for (i = 0; i < numrows; i++)
    {
      for (j = 0; j < numcols; j++)
      {
        matrixcolumn = matrix->columns[j];
        l = matrixcolumn[i];
        print_laurentpoly (l, "uv");
        printf ("; \t");
      }
      printf ("\n");
    }
    printf ("The last row contains the commutator factors.\n");
  }

  return (matrix);
}

void
print_matrix1 (struct laurentmatrix *matrix)
{
  print_matrix1_n_m (matrix->columns, matrix->numrows, matrix->numcols);
}

void
print_matrix (struct laurentmatrix *matrix, int indets)
{
  print_matrix_n_m (matrix->columns, matrix->numrows, matrix->numcols, indets);
}

void
print_matrix1_n_m (struct laurentpoly ***matrixcolumns, int numrows, int numcols)
{
  struct laurentpoly *l, **matrixcolumn;
  int i, j;

  printf ("Matrix entries:\n");
  for (i = 0; i < numrows; i++)
  {
    for (j = 0; j < numcols; j++)
    {
      matrixcolumn = matrixcolumns[j];
      l = matrixcolumn[i];
      print_laurentpoly (l, "t");
      printf ("; \t");
    }
    printf ("\n");
  }
}

void
print_matrix_n_m (struct laurentpoly ***matrixcolumns, int numrows, int numcols, int indets)
{
  struct laurentpoly *l, **matrixcolumn;
  int i, j;

  printf ("Matrix entries:\n");
  for (i = 0; i < numrows; i++)
  {
    for (j = 0; j < numcols; j++)
    {
      matrixcolumn = matrixcolumns[j];
      l = matrixcolumn[i];
      print_laurentpoly (l, "uvwxyzabcdefghijklmnopqrst");
      printf ("; \t");
    }
    printf ("\n");
  }
}

/* one indet */

struct laurentmatrix *
minor_matrix1_corank1 (struct laurentmatrix *matrix, int dropi, int dropj)
{
  struct laurentmatrix *minor;
  struct laurentpoly **matrixcolumn, **minorcolumn;
  int i, j, ii, jj;

  minor = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
  minor->numcols = matrix->numcols - 1;
  minor->numrows = matrix->numrows - 1;

  minor->columns = (struct laurentpoly ***) malloc (minor->numcols*sizeof(struct laurentpoly **));
  for (j = 0; j < minor->numcols; j++)
  {
    minor->columns[j] = (struct laurentpoly **) malloc (minor->numrows*sizeof(struct laurentpoly *));
  }

  for (j = 0, jj = 0; j < matrix->numcols; j++)
  {
    if (j == dropj) continue;
    matrixcolumn = matrix->columns[j];
    minorcolumn = minor->columns[jj++];
    for (i = 0, ii = 0; i < matrix->numrows; i++)
    {
      if (i == dropi) continue;
      minorcolumn[ii++] = matrixcolumn[i];
    }
  }
  return (minor);
}

/* two indets */

struct laurentmatrix *
minor_matrix2_corank1 (struct laurentmatrix *matrix, int dropi, int dropj)
{
  struct laurentmatrix *minor;
  struct laurentpoly **matrixcolumn, **minorcolumn;
  int i, j, ii, jj;

  minor = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
  minor->numcols = matrix->numcols - 1;
  minor->numrows = matrix->numrows - 1;

  minor->columns = (struct laurentpoly ***) malloc (minor->numcols*sizeof(struct laurentpoly **));
  for (j = 0; j < minor->numcols; j++)
  {
    minor->columns[j] = (struct laurentpoly **) malloc (minor->numrows*sizeof(struct laurentpoly *));
  }

  for (j = 0, jj = 0; j < matrix->numcols; j++)
  {
    if (j == dropj) continue;
    matrixcolumn = matrix->columns[j];
    minorcolumn = minor->columns[jj++];
    for (i = 0, ii = 0; i < matrix->numrows; i++)
    {
      if (i == dropi) continue;
      minorcolumn[ii++] = matrixcolumn[i];
    }
  }
  return (minor);
}

/*
 * free allocated space for Alexander matrix (k indets)
 */

void
laurent_free_matrix (struct laurentmatrix *matrix)
{
  struct laurentpoly *l, **matrixcolumn;
  int i, j;

  for (j = 0; j < matrix->numcols; j++)
  {
    matrixcolumn = matrix->columns[j];
    for (i = 0; i < matrix->numrows; i++)
    {
      l = matrixcolumn[i];
      if (l) free_laurentpoly (l);
    }
    free (matrixcolumn);
  }
  free (matrix->columns);
  free (matrix);
}

/*
 * compute all the determinants of kxk minors in matrix
 */

int next_subset (int *vec, int k, int n);

struct alexanderideal *
compute_invariant_factor (struct laurentpoly ***columns, int numrows, int numcols, int minordim, int indets)
{
  struct alexanderideal *ai = 0;
  struct laurentmatrix *minor;
  struct laurentpoly **column, *determinant;
  int i, j, k, numminors;
  int *veci, *vecj, s, ss;

  assert (indets >= 2);
  switch (minordim)
  {
    case 1:
      ai = (struct alexanderideal *) malloc (AI_DIM(numrows*numcols));
      ai->indets = indets;
      ai->fl2num = 0;
      ai->fl2offset = numrows*numcols;
      ai->max_generators_num = numrows*numcols;
      k = 0;
      for (j = 0; j < numcols; j++)
      {
        column = columns[j];
        for (i = 0; i < numrows; i++)
        {
          if (column[i]) ai->l[k++] = laurent_dup (column[i]);
        }
      }
      ai->l2num = k;
    break;

    default:
      numminors = binomial (numrows,minordim)*binomial(numcols,minordim);
      ai = (struct alexanderideal *) malloc (AI_DIM(numminors));
      ai->indets = indets;
      ai->fl2num = 0;
      ai->fl2offset = numminors;
      ai->max_generators_num = numminors;

      minor = (struct laurentmatrix *) malloc (sizeof (struct laurentmatrix));
      minor->numcols = minordim;
      minor->numrows = minordim;
      minor->columns = (struct laurentpoly ***) malloc (minordim*sizeof(struct laurentpoly **));
      for (j = 0; j < minordim; j++)
        minor->columns[j] = (struct laurentpoly **) malloc (minordim*sizeof(struct laurentpoly *));

      k = 0;
      vecj = (int *) malloc (minordim*sizeof(int));
      veci = (int *) malloc (minordim*sizeof(int));
      for (s = 0; s < minordim; s++) vecj[s] = s;
      do
      {
        for (s = 0; s < minordim; s++) veci[s] = s;
        do
        {
          for (s = 0; s < minordim; s++)
          {
            for (ss = 0; ss  < minordim; ss++)
            {
              minor->columns[s][ss] = columns[vecj[s]][veci[ss]];
            }
          }
          determinant = laurent_compute_determinant (minor->columns, minor->numcols, indets);
          if (determinant) ai->l[k++] = determinant;
        } while (next_subset (veci, minordim, numrows));
      } while (next_subset (vecj, minordim, numcols));
      free (vecj);
      free (veci);
      for (s = 0; s < minordim; s++) free (minor->columns[s]);
      free (minor->columns);
      free (minor);
      ai->l2num = k;
    break;
  }
  return (ai);
}

int
next_subset (int *vec, int k, int n)
{
  if (vec[k-1] < n - 1)
  {
    vec[k-1]++;
    return (1);
  }
  if (k == 0) return (0);

  if (next_subset (vec, k-1, n-1))
  {
    vec[k-1] = vec[k-2] + 1;
    return (1);
  }
  return (0);
}

void laurent_sort_entries (int num, struct laurentpoly *l[]);

struct alexanderideal *
laurent_simplify_ideal (struct alexanderideal *ai)
{
  struct laurentpoly *oldgcd, *newgcd;
  extern int principal, verbose, experimental;
  int last, i, spread, lspread;
  int linf, maxcoef, loop = 1;

  if (ai->indets == 1 && (experimental || principal == 0))
  {
    ai = groebner1 (ai);
    /* BIG WARNING:  it is not yet proved that this is really a UNIQUE Groebner basis */
    /* moreover, canonization should be performed with respect to t -> 1/t also */
    /* (simultaneously on all generators) this is not done yet */
    laurent_sort_entries (ai->l1num, ai->l);
    for (i = 0; i < ai->l1num; i++)
      laurent_canonifysign1 (ai->l[i]);
    return (ai);
  }

  /* principal = 1 in test.28 and test.29 */
  /* ai->indets > 1 in test.24 and test.25 */
  /* it seems that indets > 1 only if called directly from contour.c! */
  if (ai->indets > 1 && principal)
  {
    fprintf (stderr, "Warning: computation of smallest principal ideal is not implemented for polynomials\n");
    fprintf (stderr, "in two or more indeterminates.\n");
  }
  while (loop)
  {
    loop = 0;
    while (laurent_try_simplify_ideal (ai));

    if (ai->l1num > 1)
    {
      spread = 1;
      last = ai->l1num - 1;
      if (int_overflow_encountered)
      {
        /* chek if size of coefficients is too large */
        maxcoef = 0;
        for (i = 0; i < ai->l1num; i++)
        {
          if ((linf = laurentpoly_linf (ai->l[i])) > maxcoef) maxcoef = linf;
        }
        if (maxcoef > (INT_MAX >> (4*sizeof(int))))
        {
          start_comment ();
          printf ("WARNING: inhibiting gcd computation due to integer size\n");
          break;
        }
      }
      newgcd = laurent_gcd (spread, ai->l[last], ai->l[last-1], &lspread);
      spread = lspread;
      for (i = last-2; i >= 0; i--)
      {
        oldgcd = newgcd;
        newgcd = laurent_gcd (spread, oldgcd, ai->l[i], &lspread);
        spread = lspread;
        free (oldgcd);
      }
      if (newgcd->stemdegree > ai->l[0]->stemdegree ||
          (newgcd->stemdegree == ai->l[0]->stemdegree && abs(spread*newgcd->stem[0].l0) >= abs(ai->l[0]->stem[0].l0)) )
      {
        free (newgcd);
      } else {
        if (verbose) printf ("Information: degree reduction in ideal generators!\n");
        for (i = 0; i <= newgcd->stemdegree; i++) newgcd->stem[i].l0 *= spread;
        if (ai->l1num >= ai->max_generators_num) ai = ai_increase_size (ai);
        if (ai->l1num < ai->max_generators_num)
        {
          for (i = ai->l1num; i > 0; i--) ai->l[i] = ai->l[i-1];
          ai->l1num++;
          ai->l[0] = newgcd;
          loop = 1;
        } else {
          printf ("Warning: no space left to add new generator!\n");
          free (newgcd);
        }
      }
    }
  }

  for (i = 0; i < ai->l1num; i++)
    laurent_canonifysign1 (ai->l[i]);

  if (principal && ai->l1num > 1)
  {
    spread = 1;
    for (i = 1; i < ai->l1num; i++)
    {
      newgcd = laurent_gcd (spread, ai->l[0], ai->l[i], &lspread);
      spread = lspread;
      free (ai->l[0]);
      free (ai->l[i]);
      ai->l[0] = newgcd;
    }
    ai->l1num = 1;
    ai->spread = spread;
    return (ai);
  }

  return (ai);  /* TODO: not implemented at the moment! */
}

//#define USE_EUCLID_GCD 1

/*
 * try to simplify the ideal by removal of some generator or
 * reduction in complexity of one generator
 * no increase is allowed
 * return 1 if at least one simplification was performed
 */

int laurent_try_reduce_pair (struct laurentpoly *l1, struct laurentpoly *l2);

int
laurent_try_simplify_ideal (struct alexanderideal *ai)
{
  struct laurentpoly *l, *l1, *l2;
#ifdef USE_EUCLID_GCD
  struct laurentpoly *newgcd;
#endif
  int i, j, spread, numreductions;
  extern int verbose;

  if (verbose >= 2)
  {
    printf ("==============================\n");
    printf ("Try simplify ideal: \n");
    printout_ideal1 (ai, 0);
    printf ("==============================\n");
  }
  numreductions = spread = 0;
  if (ai->indets != 1 ) return (0);

  if (ai->l1num <= 1) return (0); // at least two generators...

  laurent_sort_entries (ai->l1num, &ai->l[0]);
  for (i = 0; i < ai->l1num; i++)
  {
    if (ai->l[i] == 0)
    {
      for (j = i; j < ai->l1num - 1; j++) ai->l[j] = ai->l[j+1];
      ai->l1num--;
      if (verbose) printf ("Ideal simplification: zero polynomial removed\n");
      return (1);
    }
  }
  for (i = 0; i < ai->l1num; i++)
  {
    l = ai->l[i];
    if (l->stemdegree == 0 && abs(l->stem[0].l0) == 1)
    {
      for (j = 0; j < ai->l1num; j++)
      {
        if (j != i && ai->l[j]) free (ai->l[j]);
        ai->l[j] = 0;
      }
      ai->l[0] = l;
      ai->l1num = 1;
      if (verbose) printf ("Ideal simplification: unit polynomial generates everything\n");
      return (1);
    }
#ifdef USE_EUCLID_GCD
    for (j = i+1; j < ai->l1num; j++)
    {
      /*
       * check whether the two polynomials generate a principal ideal
       */
      newgcd = laurent_gcd (1, ai->lx.l1[i], ai->lx.l1[j], &spread);
      if (spread == 1)
      {
        if (ai->lx.l1[i]) free (ai->lx.l1[i]);
        if (ai->lx.l1[j]) free (ai->lx.l1[j]);
        ai->lx.l1[i] = newgcd;
        ai->lx.l1[j] = 0;
        if (verbose) printf ("Ideal simplification: two pols generate a principal ideal\n");
        assert (newgcd);
        assert (newgcd->stem[0].l0);
        assert (newgcd->stem[newgcd->stemdegree].l0);
        return (1);
      } else free (newgcd);
    }
#endif
  }
  for (i = 0; i < ai->l1num; i++)
  {
    l1 = ai->l[i];
    if (l1 == 0) continue;
    assert (l1->stem[0].l0);
    assert (l1->stem[l1->stemdegree].l0);
    for (j = i+1; j < ai->l1num; j++)
    {
      /*
       * kind of euclidean division step
       */
      l2 = ai->l[j];
      if (l2 == 0) continue;
      assert (l2->stem[0].l0);
      assert (l2->stem[l2->stemdegree].l0);
      if (l1->stemdegree <= l2->stemdegree) numreductions += laurent_try_reduce_pair (l2, l1);
      if (l2->stemdegree == 0 && l2->stem[0].l0 == 0)
      {
        free (l2);
        ai->l[j] = 0;
        continue;
      }
      if (l1->stemdegree >= l2->stemdegree) numreductions += laurent_try_reduce_pair (l1, l2);
      if (l1->stemdegree == 0 && l1->stem[0].l0 == 0)
      {
        free (l1);
        ai->l[i] = 0;
        break;
      }
    }
    if (numreductions)
    {
      if (verbose) printf ("Ideal simplification: %d pairwise reductions performed\n", numreductions);
      return (numreductions);
    }
  }
  return (0);
}

/*
 * try division of l1 by a monomial multiple of l2, in case there is a stemdegree reduction
 * note that if l1 is a monomial multiple of l2, the result is a polynomial of stemdegree 0
 * and null coefficient, this should be tested by the caller function and the pointer to
 * the polynomial replaced by the null pointer
 *
 * TODO: the simplification strategy should be modified in order to avoid abnormal
 * size increase of the coefficients in the middle...
 */

int l_safety_check (int f, int maxc);

int
laurent_try_reduce_pair (struct laurentpoly *l1, struct laurentpoly *l2)
{
  int f, i, j, k, c1, c2first, c2last, quotient;
  int deltadegree;
  int s2first = 1;
  int s2last = 1;
  int numreductions = 0;
  int maxc2size;

  assert (l1);
  assert (l2);
  assert (l1->stemdegree >= l2->stemdegree);

  c1 = l1->stem[0].l0;
  c2first = l2->stem[0].l0;
  assert (c2first);
  maxc2size = laurentpoly_linf (l2);

  quotient = c1/c2first;
  if (c1 == c2first*quotient)
  {
    if ((l_safety_check (quotient, maxc2size)) == 0) return (0);
    for (i = 0; i <= l2->stemdegree; i++) l1->stem[i].l0 -= quotient*l2->stem[i].l0;
    if (laurentpoly_linf (l1) >= INT_MAX/2)
    {
      /* backtrack if coefficients grow too much */
      for (i = 0; i <= l2->stemdegree; i++) l1->stem[i].l0 += quotient*l2->stem[i].l0;
      return (0);
    }
    assert (l1->stem[0].l0 == 0);
    while (l1->stem[0].l0 == 0 && l1->stemdegree > 0)
    {
      for (i = 0; i < l1->stemdegree; i++) l1->stem[i].l0 = l1->stem[i+1].l0;
      l1->stemdegree--;
      l1->minexpon++;
    }
    while (l1->stem[l1->stemdegree].l0 == 0 && l1->stemdegree > 0)
      l1->stemdegree--;
    if (l1->stemdegree > 0)
    {
      assert (l1->stem[0].l0);
      assert (l1->stem[l1->stemdegree].l0);
    }
    return (1);
  }
  c1 = l1->stem[l1->stemdegree].l0;
  c2last = l2->stem[l2->stemdegree].l0;
  assert (c2last);
  quotient = c1/c2last;
  if (c1 == c2last*quotient)
  {
    if ((l_safety_check (quotient, maxc2size)) == 0) return (0);
    k = l1->stemdegree - l2->stemdegree;
    for (i = 0; i <= l2->stemdegree; i++) l1->stem[i + k].l0 -= quotient*l2->stem[i].l0;
    if (laurentpoly_linf (l1) >= INT_MAX/2)
    {
      /* backtrack if coefficients grow too much */
      for (i = 0; i <= l2->stemdegree; i++) l1->stem[i + k].l0 += quotient*l2->stem[i].l0;
      return (0);
    }
    assert (l1->stem[l1->stemdegree].l0 == 0);
    while (l1->stem[l1->stemdegree].l0 == 0 && l1->stemdegree > 0)
      l1->stemdegree--;
    while (l1->stem[0].l0 == 0 && l1->stemdegree > 0)
    {
      for (i = 0; i < l1->stemdegree; i++) l1->stem[i].l0 = l1->stem[i+1].l0;
      l1->stemdegree--;
      l1->minexpon++;
    }
    if (l1->stemdegree > 0)
    {
      assert (l1->stem[0].l0);
      assert (l1->stem[l1->stemdegree].l0 != 0);
    }
    return (1);
  }

  /* no reduction in stemdegree is possible, try to reduce coefficients size */

  deltadegree = l1->stemdegree - l2->stemdegree;

  if (c2first < 0) s2first = -1;
  c2first *= s2first;
  if (c2last < 0) s2last = -1;
  c2last *= s2last;

  for (k = 0; k <= deltadegree/2; k++)
  {
    if (l1->stem[k].l0 > c2first/2)
    {
      f = (l1->stem[k].l0 - c2first/2 + c2first - 1)/c2first;
      if ((l_safety_check (f, maxc2size)) == 0) return (0);
      for (i = 0, j = k; i <= l2->stemdegree; i++) l1->stem[j++].l0 -= f*s2first*l2->stem[i].l0;
      if (laurentpoly_linf (l1) >= INT_MAX/2)
      {
        /* backtrack if coefficients grow too much */
        for (i = 0, j = k; i <= l2->stemdegree; i++) l1->stem[j++].l0 += f*s2first*l2->stem[i].l0;
        return (0);
      }
      numreductions++;
    }
    if (l1->stem[k].l0 < -(c2first-1)/2)
    {
      f = (-(c2first-1)/2 - l1->stem[k].l0 + c2first - 1)/c2first;
      if ((l_safety_check (f, maxc2size)) == 0) return (0);
      for (i = 0, j = k; i <= l2->stemdegree; i++) l1->stem[j++].l0 += f*s2first*l2->stem[i].l0;
      if (laurentpoly_linf (l1) >= INT_MAX/2)
      {
        /* backtrack if coefficients grow too much */
        for (i = 0, j = k; i <= l2->stemdegree; i++) l1->stem[j++].l0 -= f*s2first*l2->stem[i].l0;
        return (0);
      }
      numreductions++;
    }
    if (k >= (deltadegree+1)/2) continue;
    if (l1->stem[l1->stemdegree - k].l0 > c2last/2)
    {
      f = (l1->stem[l1->stemdegree - k].l0 - c2last/2 + c2last - 1)/c2last;
      if ((l_safety_check (f, maxc2size)) == 0) return (0);
      for (i = 0, j = deltadegree - k; i <= l2->stemdegree; i++) l1->stem[j++].l0 -= f*s2last*l2->stem[i].l0;
      if (laurentpoly_linf (l1) >= INT_MAX/2)
      {
        /* backtrack if coefficients grow too much */
        for (i = 0, j = deltadegree - k; i <= l2->stemdegree; i++) l1->stem[j++].l0 += f*s2first*l2->stem[i].l0;
        return (0);
      }
      numreductions++;
    }
    if (l1->stem[l1->stemdegree - k].l0 < -(c2last-1)/2)
    {
      f = (-(c2last-1)/2 - l1->stem[l1->stemdegree - k].l0 + c2last - 1)/c2last;
      if ((l_safety_check (f, maxc2size)) == 0) return (0);
      for (i = 0, j = deltadegree - k; i <= l2->stemdegree; i++) l1->stem[j++].l0 += f*s2last*l2->stem[i].l0;
      if (laurentpoly_linf (l1) >= INT_MAX/2)
      {
        /* backtrack if coefficients grow too much */
        for (i = 0, j = deltadegree - k; i <= l2->stemdegree; i++) l1->stem[j++].l0 -= f*s2first*l2->stem[i].l0;
        return (0);
      }
      numreductions++;
    }
  }

  return (numreductions);
}

int
laurentpoly_linf (struct laurentpoly *l)
{
  int i;
  int norm = 0;

  for (i = 0; i <= l->stemdegree; i++)
  {
    if (abs(l->stem[i].l0) > norm) norm = abs(l->stem[i].l0);
  }
  return (norm);
}

int
l_safety_check (int f, int maxc)
{
  int safesize = INT_MAX;

  if (f) safesize /= abs(f);
  safesize /= 4;
  if (maxc > safesize)
  {
    if (int_overflow_encountered++ == 0)
    {
      start_comment ();
      printf ("WARNING: above max allowed int size, cannot complete simplification\n");
    }
    return (0);
  }
  return (1);
}

int
ll_safety_check (long long int f, long long int maxc)
{
  long long int safesize = LLONG_MAX;

  if (f < 0) f = -f;
  if (f) safesize /= f;
  safesize /= 4;
  if (maxc > safesize)
  {
    if (int_overflow_encountered++ == 0)
    {
      start_comment ();
      printf ("WARNING: above max allowed int size, cannot complete simplification\n");
    }
    return (0);
  }
  return (1);
}

void
signal_int_overflow (void)
{
  int_overflow_encountered++;
}

void laurent_sort_entries_buf (int num, struct laurentpoly *l[], struct laurentpoly *buffer[]);

void
laurent_sort_entries (int num, struct laurentpoly *l[])
{
  struct laurentpoly **buffer;

  if (num <= 1) return;

  buffer = (struct laurentpoly **) malloc (num * sizeof (struct laurentpoly *));

  laurent_sort_entries_buf (num, l, buffer);
  free (buffer);
}

void
laurent_sort_entries_buf (int num, struct laurentpoly *l[], struct laurentpoly *buffer[])
{
  struct laurentpoly *li, *lj;
  int i, j, k, s;
  int half, iwins;

  if (num <= 1) return;

  half = num/2;

  /*
   * sort each half
   */
  laurent_sort_entries_buf (half, l, buffer);
  laurent_sort_entries_buf (num - half, l + half, buffer);

  /*
   * merge...
   */

  i = 0; j = half; k = 0;
  while (k < num)
  {
    /* compare i and j entries */
    if (i >= half) {buffer[k++] = l[j++]; continue;}
    if (j >= num) {buffer[k++] = l[i++]; continue;}
    if (l[i] == 0) {buffer[k++] = l[i++]; continue;}
    if (l[j] == 0) {buffer[k++] = l[j++]; continue;}
    if (l[i]->stemdegree < l[j]->stemdegree) {buffer[k++] = l[i++]; continue;}
    if (l[i]->stemdegree > l[j]->stemdegree) {buffer[k++] = l[j++]; continue;}
    li = l[i];
    lj = l[j];
    iwins = 1;
    for (s = 0; s <= li->stemdegree; s++)
    {
      if (abs (li->stem[s].l0) < abs (lj->stem[s].l0)) break;
      if (abs (li->stem[s].l0) > abs (lj->stem[s].l0)) {iwins = 0; break;}
      if (li->stem[s].l0 < lj->stem[s].l0) break;
      if (li->stem[s].l0 > lj->stem[s].l0) {iwins = 0; break;}
    }
    if (iwins) {buffer[k++] = l[i++]; continue;}
    buffer[k++] = l[j++];
  }
  assert (i == half);
  assert (j == num);
  assert (k == num);

  for (k = 0; k < num; k++) l[k] = buffer[k];
}

/*
 * eliminate two indeterminates by conjugacy
 */

struct laurentpoly *
laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2,
             struct laurentpoly ***extradeterminantspt)
{
  extern int verbose, debug, quiet;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  struct laurentpoly **matrixcolumn;
  struct laurentpoly **extradeterminants;
  struct laurentmatrix *matrix;
  int numcols = 0;
  int lastrow, deficiency, j;
  int extrafound = 0;

  matrix = laurent_build_matrix2 (p, e1, e2);

  for (r = p->rules; r; r = r->next) numcols++;
  assert (numcols == matrix->numcols);
  deficiency = p->gennum - numcols;
  assert (deficiency >= 1);
  if (deficiency > 2) return (0);

  /*
   * tricky problem...
   * the fundamental groups of both the splitted two-components link and
   * of the figure-eight are the free group of rank 2, however in the case
   * of the figure-eight it is natural to return +1, since p(1,1) = 1 for
   * all (other) genus 2 complements, whereas in the case of the splitted
   * link Fox (e.g.) gives 0 as the Alexander polynomial, to be multiplied
   * by the fundamental ideal, giving anyway the trivial ideal.
   * It seems to be a matter of deficiency, since the splitted link has
   * deficiency 2, whereas in all other cases the deficiency is 1
   *
   * Indeed the program gives "0" as a result if we force the deficiency to
   * be one, e.g. by giving a presentation such as <a, b, aA>:
   * $ echo "fpgroup{<a,b; aA>}" | contour --nosimplify --out alexander
   * *** Warning: result can be noncanonical ***
   * Alexander polynomial:
   * 0
   *
   * This problem is solved by allowing the user to explicitly require a value for
   * the index "d" used by R.H. Fox to select one of the Alexander ideals.
   *
   * If no indication is given, then "contour" selects a default value and prints the
   * value selected on output, to avoid confusion.
   */

  lastrow = matrix->numrows - 1;
  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix->columns[j];
    if (matrixcolumn[lastrow]) extrafound++;
  }

  /* determinant of the square matrix in case of links, of the principal minor otherwise */
  determinant = laurent_compute_determinant (matrix->columns, numcols, 2);

  switch (deficiency)
  {
    case 1:  /* link with two components */
    if (determinant == 0)
    {
      if (extradeterminantspt) *extradeterminantspt = 0;
    } else {
      if (extradeterminantspt == 0)
      {
        if (!quiet)
        {
          printf ("Result must be multiplied by the fundamental ideal!\n");
        }
      } else {
        *extradeterminantspt = extradeterminants =
             (struct laurentpoly **) malloc (sizeof(struct laurentpoly *));
        extradeterminants[0] = determinant;
        determinant = 0;
      }
    }
    break;

    case 2:  /* bouquet of two loops */
    assert (matrix->numrows == numcols + 1);
    if (extrafound && extradeterminantspt == 0)
    {
      printf ("FATAL: extra commutator factors found but computation forbidden!\n");
      printf ("  printed result will be incorrect\n");
    }
    if (extrafound && extradeterminantspt)
    {
      *extradeterminantspt = extradeterminants =
           (struct laurentpoly **) malloc (numcols*sizeof(struct laurentpoly *));
      for (j = 0; j < numcols; j++)
      {
        extradeterminants[j] = laurent_minor_determinant (matrix->columns, numcols, j, 2);
      }
    }
    if (extrafound == 0 && extradeterminantspt) *extradeterminantspt = 0;
    break;

    default:
    printf ("FATAL: invalid deficiency: %d\n", deficiency);
  }

  laurent_free_matrix (matrix);

  return (determinant);
}

/*
 * Compute matrix entries for Alexander polynomial computation (one indet)
 */

struct laurentpoly *
laurent_get_exp_sum1 (struct presentationrule *r, int n, int gconj)
{
  int k, d;
  int stemdegree;
  int runningexp, minexp, maxexp;
  struct laurentpoly *l;

  runningexp = 0;
  minexp = INT_MAX;
  maxexp = INT_MIN;
  assert (gconj > 0);
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -gconj) runningexp--;
    if (abs(r->var[k]) == n)
    {
      if (runningexp < minexp) minexp = runningexp;
      if (runningexp > maxexp) maxexp = runningexp;
    }
    if (r->var[k] == gconj) runningexp++;
  }

  if (minexp > maxexp) return (0);
  stemdegree = maxexp - minexp;
  l = (struct laurentpoly *) malloc (POLYSIZE(stemdegree + 1));
  l->indets = 1;
  l->stemdegree = stemdegree;
  l->minexpon = minexp;
  for (d = 0; d <= stemdegree; d++) l->stem[d].l0 = 0;

  runningexp = 0;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -gconj) runningexp--;
    if (abs(r->var[k]) == n)
    {
      d = runningexp - minexp;
      if (r->var[k] > 0) l->stem[d].l0++;
        else l->stem[d].l0--;
    }
    if (r->var[k] == gconj) runningexp++;
  }

  assert (l);
  if (laurent_normalize(l) == 0)
  {
    free (l);
    return (0);
  }
  assert (l && l->stem[0].l0);
  return (l);
}

/*
 * Compute matrix entries for Alexander polynomial computation (two indeterminates)
 */

struct laurentpoly *
laurent_get_exp_sum (struct presentationrule *r, int n, int e1, int e2)
{
  int k, d;
  int stemdegree;
  int re1, re2, minexp1, minexp2, maxexp1, maxexp2;
  struct laurentpoly *l;

  assert (e1 > 0 && e2 > 0);
  re1 = re2 = 0;
  minexp1 = minexp2 = INT_MAX;
  maxexp1 = maxexp2 = INT_MIN;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -e1) re1--;
    if (r->var[k] == -e2) re2--;
    if (abs(r->var[k]) == n)
    {
      if (re1 < minexp1) minexp1 = re1;
      if (re2 < minexp2) minexp2 = re2;
      if (re1 > maxexp1) maxexp1 = re1;
      if (re2 > maxexp2) maxexp2 = re2;
    }
    if (r->var[k] == e1) re1++;
    if (r->var[k] == e2) re2++;
  }

  if (minexp1 > maxexp1 || minexp2 > maxexp2) return (0);
  stemdegree = maxexp2 - minexp2;
  l = (struct laurentpoly *) malloc (POLYSIZE(stemdegree + 1));
  l->indets = 2;
  l->stemdegree = stemdegree;
  l->minexpon = minexp2;
  for (d = 0; d <= stemdegree; d++) l->stem[d].lx = 0;

  re1 = re2 = 0;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -e1) re1--;
    if (r->var[k] == -e2) re2--;
    if (abs(r->var[k]) == n)
    {
      d = re2 - minexp2;
      if (r->var[k] > 0) l->stem[d].lx = laurentpoly1_addmonom (l->stem[d].lx, re1, 1);
        else l->stem[d].lx = laurentpoly1_addmonom (l->stem[d].lx, re1, -1);
    }
    if (r->var[k] == e1) re1++;
    if (r->var[k] == e2) re2++;
  }

  assert (l);
  if (laurent_normalize(l) == 0)
  {
    free_laurentpoly (l);
    return (0);
  }
  assert (l && l->stem[0].lx);
  return (l);
}

/*
 * Compute matrix entries on pairs of indeterminates using Fox
 * second derivatives, assuming expsum is zero for both
 * It is \partial_{x2}\partial_{x1} w (1,u,v)
 * where u is substituted in place of x1 and v in place of x2,
 * 1 in all others
 *
 * polynomial in two indets
 */

struct laurentpoly *
laurent_mixed_derivative2 (struct presentationrule *r, int x1, int x2)
{
  int k, d;
  int stemdegree;
  int re1, re2, minexp1, minexp2, maxexp1, maxexp2;
  struct laurentpoly *l;

  assert (x1 > 0 && x2 > 0);
  re1 = re2 = 0;
  minexp1 = minexp2 = INT_MAX;
  maxexp1 = maxexp2 = INT_MIN;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -x1) re1--;
    if (r->var[k] == -x2) re2--;
    if (abs(r->var[k]) == x2)
    {
      if (re1 < minexp1) minexp1 = re1;
      if (re2 < minexp2) minexp2 = re2;
      if (re1 > maxexp1) maxexp1 = re1;
      if (re2 > maxexp2) maxexp2 = re2;
    }
    if (r->var[k] == x1) re1++;
    if (r->var[k] == x2) re2++;
  }
  assert (re1 == 0);
  assert (re2 == 0);

  if (minexp1 > maxexp1 || minexp2 > maxexp2) return (0);
  stemdegree = maxexp2 - minexp2;
  l = (struct laurentpoly *) malloc (POLYSIZE(stemdegree + 1));
  l->indets = 2;
  l->stemdegree = stemdegree;
  l->minexpon = minexp2;
  for (d = 0; d <= stemdegree; d++) l->stem[d].lx = 0;

  re1 = 0;
  re2 = -minexp2;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -x1) re1--;
    if (r->var[k] == -x2) re2--;
    if (r->var[k] == x2) l->stem[re2].lx = laurentpoly1_addmonom (l->stem[re2].lx, re1, -re1);
    if (r->var[k] == -x2) l->stem[re2].lx = laurentpoly1_addmonom (l->stem[re2].lx, re1, re1);
    //if (abs(r->var[k]) == x2)
    //{
    //  d = re2 - minexp2;
    //  coef = -re1;
    //  if (r->var[k] < 0) coef = re1;
    //  l->stem[d] = laurentpoly1_addmonom (l->stem[d], re1, coef);
    //}
    if (r->var[k] == x1) re1++;
    if (r->var[k] == x2) re2++;
  }

  assert (l);
  if (laurent_normalize(l) == 0)
  {
    free_laurentpoly (l);
    return (0);
  }
  assert (l && l->stem[0].lx);
  return (l);
}

/*
 * Compute common factor on pairs of indeterminates
 * expsum is zero for both
 */

struct laurentpoly *
laurent_common_factor (struct presentationrule *r, int x1, int x2)
{
  int k;
  struct laurentpoly *l, *poly1, *poly2;
  struct laurentpoly *res, *addres, *mulres, *uminus1; /* these must be polynomials in one indet */

  poly1 = laurent_get_exp_sum (r, x1, x1, x2);
  if (poly1) assert (poly1->indets == 2);
  /* this should be divisible by (x2 - 1) */
  res = laurent_sum_coefficients2 (poly1);
  assert (res == 0);

  /* now divide by x2-1 */
  if (poly1)
  {
    for (k = poly1->stemdegree - 1; k >= 0; k--)
    {
      addres = laurent_add (poly1->stem[k+1].lx, poly1->stem[k].lx);
      if (poly1->stem[k].lx) free_laurentpoly (poly1->stem[k].lx);
      poly1->stem[k].lx = addres;
    }
    assert (poly1->stem[0].lx == 0);
    poly1->minexpon--;
    assert (laurent_normalize(poly1));
  }

  poly2 = laurent_get_exp_sum (r, x2, x1, x2);
  if (poly2) assert (poly2->indets == 2);
  /* this should be divisible by (x1 - 1) */
  res = laurent_sum_each_coefficient2 (poly2);
  assert (res == 0);

  /* poly2 + (u-1)*poly1 should be zero */
  /* doing it by hands... */
  if (poly2 == 0) assert (poly1 == 0);
  if (poly2)
  {
    assert (poly2->minexpon == poly1->minexpon);
    assert (poly2->stemdegree == poly1->stemdegree);
    uminus1 = (struct laurentpoly *) malloc (POLYSIZE(2));
    uminus1->indets = 1;
    uminus1->minexpon = 0;
    uminus1->stemdegree = 1;
    uminus1->stem[0].l0 = -1;
    uminus1->stem[1].l0 = 1;
    for (k = 0; k <= poly1->stemdegree; k++)
    {
      mulres = laurent_mul1 (poly1->stem[k].lx, uminus1);
      addres = laurent_add (poly2->stem[k].lx, mulres);
      if (mulres) free_laurentpoly (mulres);
      assert (addres == 0);
    }
    free_laurentpoly (uminus1);
  }

  l = poly1;

  free_laurentpoly (poly2);
  return (l);
}

/*
 * compute the determinant of the matrix (one indet polynomials)
 */

struct laurentpoly *
laurent_compute_determinant1 (struct laurentpoly ***matrix, int n)
{
  extern int debug;
  int i, ii, jj, i1;
  int sign;
  struct laurentpoly *determinant = 0, *subdeterminant, *product, *sum;
  struct laurentpoly ***submatrix;
  struct laurentpoly **matrixcol, **submatrixcol, **firstcolumn;

  assert (n >= 0);
  if (n == 0)
  {
    determinant = (struct laurentpoly *) malloc (POLYSIZE(1));
    determinant->indets = 1;
    determinant->minexpon = 0;
    determinant->stemdegree = 0;
    determinant->stem[0].l0 = 1;
    return (determinant);
  }
  if (n == 1) return (laurent_dup(matrix[0][0]));

  /*
   * developping the determinant about the first column
   */

  if (debug >= 2) printf ("Computing determinant of a %dx%d matrix\n", n, n);
  if (debug >= 2) print_matrix1_n_m (matrix, n, n);
  submatrix = (struct laurentpoly ***) malloc ( (n-1)*sizeof (struct laurentpoly **) );
  for (jj = 0; jj < n-1; jj++)
    submatrix[jj] = (struct laurentpoly **) malloc ( (n-1)*sizeof (struct laurentpoly *) );

  firstcolumn = matrix[0];
  sign = 1;
  for (i = 0, sign = 1; i < n; i++, sign = -sign)
  {
    if (firstcolumn[i] == 0) continue;
    for (jj = 0; jj < n - 1; jj++)
    {
      submatrixcol = submatrix[jj];
      matrixcol = matrix[jj + 1];
      for (i1 = 0, ii = 0; i1 < n; i1++)
      {
        if (i1 == i) continue;
        submatrixcol[ii++] = matrixcol[i1];
      }
    }

    subdeterminant = laurent_compute_determinant1 (submatrix, n-1);
    product = laurent_mul1 (subdeterminant, firstcolumn[i]);
    if (debug >= 2)
    {
      printf ("Development row: %d of %d\n", i, n);
      printf ("element: "); print_laurentpoly (firstcolumn[i], "t");
      printf ("\nminor determinant: "); print_laurentpoly (subdeterminant, "t");
      printf ("\n product: "); print_laurentpoly (product, "t");
      printf ("\n");
    }
    if (subdeterminant) free (subdeterminant);
    if (sign < 0) laurent_negate (product);
    sum = laurent_add (determinant, product);
    if (product && product != sum) free (product);
    if (determinant && sum != determinant) free (determinant);
    determinant = sum;
  }

  for (jj = 0; jj < n-1; jj++) free (submatrix[jj]);
  free (submatrix);

  if (debug >= 2)
  {
    printf ("Resulting determinant of %d matrix: ", n);
    print_laurentpoly (determinant, "t");
    printf ("\n");
  }
  return (determinant);
}

/*
 * compute the determinant of the matrix (k indeterminates)
 */

struct laurentpoly *
laurent_compute_determinant (struct laurentpoly ***matrix, int n, int k)
{
  int i, ii, jj, i1;
  int sign;
  struct laurentpoly *determinant = 0, *subdeterminant, *product, *sum;
  struct laurentpoly ***submatrix;
  struct laurentpoly **matrixcol, **submatrixcol, **firstcolumn;

  assert (n >= 0);
  if (n == 0) return (laurent_buildconstant (k, 1));
  if (n == 1) return (laurent_dup(matrix[0][0]));

  /*
   * developping the determinant about the first column
   */

  submatrix = (struct laurentpoly ***) malloc ( (n-1)*sizeof (struct laurentpoly **) );
  for (jj = 0; jj < n-1; jj++)
    submatrix[jj] = (struct laurentpoly **) malloc ( (n-1)*sizeof (struct laurentpoly *) );

  firstcolumn = matrix[0];
  sign = 1;
  for (i = 0, sign = 1; i < n; i++, sign = -sign)
  {
    if (firstcolumn[i] == 0) continue;
    for (jj = 0; jj < n - 1; jj++)
    {
      submatrixcol = submatrix[jj];
      matrixcol = matrix[jj + 1];
      for (i1 = 0, ii = 0; i1 < n; i1++)
      {
        if (i1 == i) continue;
        submatrixcol[ii++] = matrixcol[i1];
      }
    }

    subdeterminant = laurent_compute_determinant (submatrix, n-1, k);
    product = laurent_mul (subdeterminant, firstcolumn[i]);
    if (subdeterminant) free_laurentpoly (subdeterminant);
    if (sign < 0) laurent_negate (product);
    sum = laurent_add (determinant, product);
    if (product && product != sum) free_laurentpoly (product);
    if (determinant && sum != determinant) free_laurentpoly (determinant);
    determinant = sum;
  }

  for (jj = 0; jj < n-1; jj++) free (submatrix[jj]);
  free (submatrix);

  return (determinant);
}

/*
 * compute the determinant of the matrix (k indeterminates)
 * with row 
 */

struct laurentpoly *
laurent_minor_determinant (struct laurentpoly ***matrix, int numcols, int row_to_subst, int k)
{
  int j;
  struct laurentpoly **matrixcolumn;
  struct laurentpoly *saved, *determinant;

  /* exchange rows */
  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    saved = matrixcolumn[row_to_subst];
    matrixcolumn[row_to_subst] = matrixcolumn[numcols];
    matrixcolumn[numcols] = saved;
  }

  determinant = laurent_compute_determinant (matrix, numcols, k);

  /* exchange back */
  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    saved = matrixcolumn[row_to_subst];
    matrixcolumn[row_to_subst] = matrixcolumn[numcols];
    matrixcolumn[numcols] = saved;
  }

  return (determinant);
}

/*
 * canonification procedure for ideals resulting as mixed
 * polinomial + extra times fundamental ideal
 * return 0 if cannot canonify
 */

int
canonify_ideal2 (struct laurentpoly **wpt, struct laurentpoly **wi, int winum)
{
  if ((*wpt != 0) && (wi!= 0)) return (0);  // for now do not try to canonify
  if (winum > 1) return (0);

  if (*wpt) return (base_canonify2 (wpt));
  return (base_canonify2 (wi));
}

/*
 * full canonification with base-equivalence
 * (e.g. for a principal ideal)
 */

int
base_canonify2 (struct laurentpoly **ppt)
{
  extern int nobasecanonify;
  int suppdim;

  if (ppt == 0) return (1);
  suppdim = laurent_suppdimx2 (*ppt);
 
  if (suppdim < 0) return (1);

  if (nobasecanonify) return (0);
  switch (suppdim)
  {
    case 0:
    case 1:
    return (base_canonify2_onedim (ppt));
    break;

    case 2:
    return (base_canonify2_twodim (ppt));
    return (0);

    default:
    printf ("UNKNOWN DIMENSION: %d\n", suppdim);
  }

  return (0);
}

/*
 * canonification in the special one-dim case
 */

int
base_canonify2_onedim (struct laurentpoly **ppt)
{
  extern int verbose;
  struct laurentpoly *p;
  struct laurentpoly *p1, *newp1;  /* one indet polynomials */
  int origy, dx, dy, k, xk;
  int mygcd, xstep, ystep;

  p = *ppt;
  assert (p && p->indets == 2);
  if (p->stemdegree == 0)
  {
    laurent_canonify1 (p->stem[0].lx);
    laurent_canonifysign1 (p->stem[0].lx);
    return (1);  // already canonical!
  }

  p1 = p->stem[0].lx;
  origy = p1->minexpon;  // orig u

  dx = p->stemdegree;    // dv
  assert (dx);
  p1 = p->stem[p->stemdegree].lx;
  assert (p1);
  assert (p1->stemdegree == 0);
  dy = p1->minexpon - origy;  // du

  mygcd = gcd (dx, dy);  /* this is the degree of the resulting polynomial */
  xstep = dx/mygcd;
  ystep = dy/mygcd;
  assert (xstep > 0);

  assert (mygcd >= 1);

  newp1 = (struct laurentpoly *) malloc (POLYSIZE(mygcd+1));
  newp1->indets = 1;
  newp1->stemdegree = mygcd;

  for (k = 0; k <= mygcd; k++) newp1->stem[k].l0 = 0;
  for (k = 0, xk = 0; k <= mygcd; k++, xk += xstep)
  {
    /* fill up new polynomial */
    p1 = p->stem[xk].lx;
    if (p1 == 0) continue;
    assert (p1->minexpon == origy + k*ystep);
    assert (p1->stemdegree == 0);
    newp1->stem[k].l0 = p1->stem[0].l0;
  }
  newp1->minexpon = 0;

  /* free mem space of old poly */
  for (xk = 0; xk <= p->stemdegree; xk++)
  {
    if (xk % xstep) assert (p->stem[xk].lx == 0);
    if (p->stem[xk].lx) free_laurentpoly (p->stem[xk].lx);
  }
  p->minexpon = 0;
  p->stemdegree = 0;
  p->stem[0].lx = newp1;
  laurent_canonifysign1 (p->stem[0].lx);

  return (1);
}

/*
 * canonification in the two-dim case
 */

int
base_canonify2_twodim (struct laurentpoly **ppt)
{
  struct laurentpoly *origp, *newp, *optp;
  struct laurentpoly *tempp;
  extern int verbose;
  int x1[2], x2[2], x3[2], du2, dv2, du3, dv3;
  int origtotdegree, newtotdegree, optdegree;
  int i, k;
  int xiu[3], xiv[3], sumxiu, sumxiv, amin, amax, bmin, bmax;
  int area2;  // area of the triangle with those vertices
  int a, b, c, d, bm[2][2];

  origp = *ppt;
  assert (origp);
  assert (origp->stem[0].lx);
  laurent_canonify_exponents (origp);
  optp = laurent_dup (origp);
  origtotdegree = laurentx_totdegree (origp);

  laurent_getthree (*ppt, x1, x2, x3);

  /* translate x1 to the origin */
  du2 = x2[0] - x1[0];
  dv2 = x2[1] - x1[1];
  du3 = x3[0] - x1[0];
  dv3 = x3[1] - x1[1];
  area2 = du2*dv3 - du3*dv2;
  if (area2 < 0)
  {
    du2 = x3[0] - x1[0];
    dv2 = x3[1] - x1[1];
    du3 = x2[0] - x1[0];
    dv3 = x2[1] - x1[1];
    area2 = du2*dv3 - du3*dv2;
  }

  if (verbose > 1)
  {
    printf ("total degree of original translated polynomial: %d\n", origtotdegree);
    printf ("three non-aligned points in the support:\n");
    printf ("  [u,v]=[%d,%d],[%d,%d],[%d,%d]\n", x1[0], x1[1], x2[0], x2[1], x3[0], x3[1]);
    printf ("twice area of triangle: %d\n", area2);
  }

  /* computing xi[0], xi[1], xi[2]: edges of the "octant" cone */
  xiu[0] = dv2 - dv3;
  xiv[0] = du3 - du2;
  xiu[1] = dv3;
  xiv[1] = -du3;
  xiu[2] = -dv2;
  xiv[2] = du2;
  sumxiu = xiu[0] + xiu[1] + xiu[2];
  sumxiv = xiv[0] + xiv[1] + xiv[2];

  /*
   * computing bounds for matrix B entries
   * f*amin <= a,c <= f*amax
   * f*bmin <= b,d <= f*bmax
   *
   * where f = origtotdegree/area2
   */
  amin = amax = bmin = bmax = 0;
  if (sumxiu < amin) amin = sumxiu;
  if (sumxiu > amax) amax = sumxiu;
  if (sumxiv < bmin) bmin = sumxiv;
  if (sumxiv > bmax) bmax = sumxiv;

  for (i = 0; i < 3; i++)
  {
    if (xiu[i] < amin) amin = xiu[i];
    if (xiu[i] > amax) amax = xiu[i];
    if (xiv[i] < bmin) bmin = xiv[i];
    if (xiv[i] > bmax) bmax = xiv[i];
    if (sumxiu - xiu[i] < amin) amin = sumxiu - xiu[i];
    if (sumxiu - xiu[i] > amax) amax = sumxiu - xiu[i];
    if (sumxiv - xiv[i] < bmin) bmin = sumxiv - xiv[i];
    if (sumxiv - xiv[i] > bmax) bmax = sumxiv - xiv[i];
  }

  if (verbose)
  {
    printf ("Bound on matrix entry 'a' (and 'c'): ");
    printf ("%d*%d <= %d*a <= %d*%d\n", origtotdegree, amin, area2, origtotdegree, amax);
    printf ("Bound on matrix entry 'b' (and 'd'): ");
    printf ("%d*%d <= %d*b <= %d*%d\n", origtotdegree, bmin, area2, origtotdegree, bmax);
  }

  /* walk through feasible base-change matrices */

  optdegree = origtotdegree;
  for (k = 0; k <= optdegree; k++)
  {
    for (a = k*amin/area2 ; a <= k*amax/area2; a++)
    {
      for (b = k*bmin/area2; b <= k*bmax/area2; b++)
      {
        for (c = k*amin/area2; c <= k*amax/area2; c++)
        {
          for (d = k*bmin/area2; d <= k*bmax/area2; d++)
          {
            bm[0][0] = a;
            bm[0][1] = b;
            bm[1][0] = c;
            bm[1][1] = d;
            if (isinvertible_base(bm) == 0) continue;
            newp = laurent_dup (origp);
            newp = base_change2 (newp, bm);
            newtotdegree = laurentx_totdegree (newp);
            laurent_canonify_exponents (newp);
            if (newtotdegree > optdegree)
            {
              free_laurentpoly (newp);
              continue;
            }
            if (newtotdegree < optdegree)
            {
              tempp = newp;
              newp = optp;
              optp = tempp;
              if (verbose > 1) printf ("Decreasing optdegree from %d to %d\n", optdegree, newtotdegree);
              optdegree = newtotdegree;
            } else {
              if (laurentx_lexicocompare (newp, optp) < 0)
              {
                if (verbose > 1) printf ("same totdegree, but better comparison\n");
                tempp = newp;
                newp = optp;
                optp = tempp;
              }
            }
            free_laurentpoly (newp);
            if (verbose > 1)
            {
              printf ("Found feasible matrix: [%d %d; %d %d]\n", a, b, c, d);
              printf ("New polynomial: ");
              print_laurentpoly (optp, "uv");
              printf (";\n");
            }
          }
        }
      }
    }
  }

  free_laurentpoly (origp);
  laurent_canonifysign2 (optp);
  *ppt = optp;
  return (1);
}

/*
 * compute dimension of support
 */

int
laurent_suppdimx2 (struct laurentpoly *p)
{
  struct laurentpoly *l1;   /* one indet polynomial */
  int dirfound = 0;
  int k, origx, origy, dx, dy, ddx, ddy;

  if (p == 0) return (-1);  //empty support
  assert (p->indets == 2);

  if (p->stemdegree == 0)
  {
    l1 = p->stem[0].lx;
    if (l1->stemdegree == 0) return (0);
  }

  if (p->stemdegree == 0) return (1);

  /* let's see if dimension is one... */
  l1 = p->stem[0].lx;
  assert (l1);
  if (l1->stemdegree > 0) return (2);
  origx = 0;
  origy = l1->minexpon;

  for (k = 1; k <= p->stemdegree; k++)
  {
    if (p->stem[k].lx == 0) continue;
    l1 = p->stem[k].lx;
    if (l1->stemdegree > 0) return (2);
    if (dirfound == 0)
    {
      dx = k - origx;
      dy = l1->minexpon - origy;
      dirfound = 1;
      continue;
    }
    ddx = k - origx;
    ddy = l1->minexpon - origy;
    /* check if (dx,dy) and (ddx,ddy) have same direction */
    if (dx*ddy != dy*ddx) return (2);
  }
  return (1);
}

/*
 * get three non-aligned points in the support
 */

void
laurent_getthree (struct laurentpoly *p, int *x1, int *x2, int *x3)
{
  struct laurentpoly *l1, *l2; /* one indet polynomials */
  int dirfound = 0;
  int k, origv, origu, dv, du, ddv, ddu;

  assert (p != 0);
  assert (p->indets == 2);
  assert (p->stemdegree > 0); // two different exponents for v

  l1 = p->stem[0].lx;
  assert (l1);
  if (l1->stemdegree > 0)     // two different exponents for u
  {
    assert (l1->stem[0].l0 && l1->stem[l1->stemdegree].l0);
    x1[0] = l1->minexpon;     // u exponent of first point
    x2[0] = l1->minexpon + l1->stemdegree;  // u exponent of second point
    x1[1] = x2[1] = p->minexpon;  // v exponent of first and second point
    l2 = p->stem[p->stemdegree].lx;
    assert (l2);
    assert (l2->stem[0].l0);
    x3[0] = l2->minexpon;
    x3[1] = p->minexpon + p->stemdegree;
    return;
  }
  origv = 0;  // must add offset p->minexpon
  origu = l1->minexpon;
  x1[0] = origu;
  x1[1] = origv + p->minexpon;

  for (k = 1; k <= p->stemdegree; k++)
  {
    if (p->stem[k].lx == 0) continue;
    l1 = p->stem[k].lx;
    if (l1->stemdegree > 0)
    {
      assert (l1->stem[0].l0 && l1->stem[l1->stemdegree].l0);
      x1[0] = l1->minexpon;     // u exponent of first point
      x2[0] = l1->minexpon + l1->stemdegree;  // u exponent of second point
      x1[1] = x2[1] = p->minexpon + k;  // v exponent of first and second point
      l2 = p->stem[0].lx;
      assert (l2->stem[0].l0);
      x3[0] = l2->minexpon;
      x3[1] = p->minexpon;
      return;
    }
    if (dirfound == 0)
    {
      x2[0] = l1->minexpon;
      x2[1] = k + p->minexpon;
      dv = k - origv;
      du = l1->minexpon - origu;
      dirfound = 1;
      continue;
    }
    ddv = k - origv;
    ddu = l1->minexpon - origu;
    /* check if (dv,du) and (ddv,ddu) have same direction */
    if (dv*ddu != du*ddv)
    {
      x3[0] = l1->minexpon;
      x3[1] = k + p->minexpon;
      return;
    }
  }
  assert (0);
}

/*
 * apply a change of base to an alexander ideal
 */

struct laurentpoly *
base_change2 (struct laurentpoly *l, int matrixb[2][2])
{
  struct laurentpoly *newl = 0;
  struct laurentpoly *lu;  /* one indet poly */
  int a, b, c, d, iu, iv, degu, degv, newdegu, newdegv;
  int coef;

  a = matrixb[0][0];
  b = matrixb[0][1];
  c = matrixb[1][0];
  d = matrixb[1][1];

  assert (isinvertible_base(matrixb));

  if (l == 0) return (0);
  assert (l->indets == 2);

  /* loop through monomials of l */
  for (iv = 0; iv <= l->stemdegree; iv++)
  {
    degv = l->minexpon + iv;
    lu = l->stem[iv].lx;
    if (lu == 0) continue;
    for (iu = 0; iu <= lu->stemdegree; iu++)
    {
      degu = lu->minexpon + iu;
      coef = lu->stem[iu].l0;
      if (coef == 0) continue;
      newdegu = a*degu + b*degv;
      newdegv = c*degu + d*degv;
      newl = laurentpoly2_addmonom (newl, newdegu, newdegv, coef);
    }
  }

  free_laurentpoly (l);

  return (newl);
}

/*
 * randomly change base of Z^2 for the alexander ideal
 */

void
shuffle_poly2 (struct laurentpoly **lpt, struct laurentpoly **extradets, int extranum)
{
  extern int verbose;
  int b[2][2];
  int i, j;

  base_random (b);
  if (verbose)
  {
    printf ("random matrix B:\n");
    for (i = 0; i < 2; i++)
    {
      printf("[");
      for (j = 0; j < 2; j++)
      {
        if (j > 0) printf ("\t");
        printf ("%d", b[i][j]);
      }
      printf ("]\n");
    }
  }

  *lpt = base_change2 (*lpt, b);
  for (j = 0; j < extranum; j++)
  {
    if ((extradets == 0) || (extradets[j] == 0)) continue;
    extradets[j] = base_change2(extradets[j], b);
  }

  if (verbose)
  {
    printf ("Alexander ideal after shuffling:\n");
    printf ("-------------\n");
    printf ("Main polynomial:\n");
    print_laurentpoly (*lpt, "uv");
    printf ("\n");
    for (j = 0; j < extranum; j++)
    {
      printf ("Fundamental ideal factor %d:\n", j);
      print_laurentpoly (extradets[j], "uv");
      printf ("\n");
    }
    printf ("-------------\n");
  }
}

/*
 * generate a random element in GL(2,Z)
 */

#define B_RAND_MAX 25

void
base_random (int b[2][2])
{
  int i, j;

  while (1)
  {
    for (i = 0; i < 2; i++)
    {
      for (j = 0; j < 2; j++)
      {
        b[i][j] = (rand() % (2*B_RAND_MAX)) - B_RAND_MAX;
      }
    }
    if (isinvertible_base (b)) break;
  }
  return;
}

/*
 * check if determinant is +- 1
 */

int
isinvertible_base (int b[2][2])
{
  int det;

  det = b[0][0]*b[1][1] - b[0][1]*b[1][0];

  if (abs(det) == 1) return (1);
  return (0);
}

/*
 * read an alexander ideal from file
 * return number of indeterminates
 */

struct alexanderideal *
read_alexander_ideal (FILE *file)
{
  char ch, indet_names[2];
  struct alexanderideal *ai;
  struct laurentpoly *l2;
  struct laurentpoly *l1;
  int sign, tok, ffound;

  tok = gettoken (file);

  assert (tok == TOK_ALEXANDER || tok == TOK_IDEAL);

  tok = gettoken (file);
  assert (tok == TOK_LPAREN);
  ai = (struct alexanderideal *) malloc (AI_DIM(IDEAL_DEF_GENERATORS_NUM));
  ai->max_generators_num = IDEAL_DEF_GENERATORS_NUM;
  ai->l1num = ai->l2num = ai->fl2num = 0;
  ai->fl2offset = IDEAL_DEF_GENERATORS_NUM/3;  /* assume there are more 'F' polynomials */
  ai->spread = 1;
  ai->val = 0;
  ai->indets = read_generators_list (file, indet_names, 2);

  tok = gettoken (file);
  assert (tok == TOK_RPAREN);

  tok = gettoken (file);
  assert (tok == TOK_LBRACE);

  switch (ai->indets)
  {
    case 0:
      sign = 1;
      ch = mygetchar (file);
      if (ch == '-') sign = -1;
      if (ch != '+' && ch != '-') ungetc (ch, file);
      ai->val = sign*get_unsignednum (file);
      ch = mygetchar (file);
      if (ch == '}') break;
      assert (ch == ';');
      ch = mygetchar (file);
      assert (ch == '}');
      break;

    case 1:
      while (1)
      {
        ch = mygetchar (file);
        if (ch == '}') break;
        ungetc (ch, file);
        l1 = read_laurentpoly1 (file, indet_names);
        ai->l[ai->l1num++] = l1;
      }
    break;

    case 2:
      while (1)
      {
        ffound = 0;
        ch = mygetchar (file);
        if (ch == '}') break;
        if (ch == 'F')
        {
          ch = mygetchar (file);
          assert (ch == ':');
          ffound = 1;
        } else ungetc (ch, file);
        l2 = read_laurentpoly2 (file, indet_names);
        if (ffound)
        {
          if (ai->fl2offset + ai->fl2num >= ai->max_generators_num) ai = ai_increase_size (ai);
          assert (ai->fl2offset + ai->fl2num < ai->max_generators_num);
          ai->l[ai->fl2offset + ai->fl2num++] = l2;
        } else {
          assert (ai->l2num < ai->fl2offset);
          ai->l[ai->l2num++] = l2;
        }
      }
    break;

    default:
      printf ("Reading an Alexander ideal from file with %d indeterminates is not implemented!\n", ai->indets);
      free (ai);
      return (0);
  }

  return (ai);
}

/*
 * reallocate ideal to accomodate for more generators
 */

struct alexanderideal *
ai_increase_size (struct alexanderideal *ai)
{
  struct alexanderideal *newai;
  extern int debug;
  int i, fincrease;
  int oldsize = ai->max_generators_num;
  int newsize = 2 + oldsize + oldsize/2;  /* circa 50% increase in size */

  if (debug) printf ("IDEAL: size increase from %d to %d\n", oldsize, newsize);
  assert (ai->indets > 0);
  assert (ai->indets <= 2);
  newai = (struct alexanderideal *) malloc (AI_DIM(newsize));
  newai->indets = ai->indets;
  newai->max_generators_num = newsize;
  newai->spread = ai->spread;
  newai->val = ai->val;
  switch (ai->indets)
  {
    case 1:
      newai->l1num = ai->l1num;
      for (i = 0; i < ai->l1num; i++) newai->l[i] = ai->l[i];
      break;

    case 2:
      assert (newsize - oldsize >= 2);
      fincrease = (newsize - oldsize)/3;
      if (fincrease <= 0) fincrease = 1;
      assert (newsize - oldsize - fincrease > 0);
      newai->fl2offset = ai->fl2offset + fincrease;
      newai->l2num = ai->l2num;
      newai->fl2num = ai->fl2num;
      for (i = 0; i < ai->l2num; i++) newai->l[i] = ai->l[i];
      for (i = 0; i < ai->fl2num; i++) newai->l[i+newai->fl2offset] = ai->l[i+ai->fl2offset];
      break;
  }
  free (ai);
  return (newai);
}

/* one-line comment */

void
start_comment (void)
{
  extern int outformat;

  switch (outformat)
  {
    case OUTFORMAT_MACAULAY2:
    printf ("-- ");
    break;

    case OUTFORMAT_APPCONTOUR:
    default:
    printf ("# ");
    break;
  }
}
