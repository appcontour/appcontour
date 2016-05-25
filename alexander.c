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
#include "alexander.h"
#include "laurent.h"
#include "parser.h"

static int int_overflow_encountered = 0;

int
alexander (struct presentation *p)
{
  extern int verbose, quiet, foxd, shuffle, outformat;
  struct presentationrule *r;
  struct laurentpoly *determinant;
  struct laurentpoly2 *determinant2;
  struct laurentpoly2 **extradeterminants;
  struct alexanderideal *ai;
  int i, j, sum, rank, matrixrank, deficiency;
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
  rank = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    assert (sum >= 0);
    if (sum) matrixrank = i;
      else gconj = i;
    if (sum && sum != 1)
    {
      if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {\n");
      if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
      printf ("Cannot compute Alexander polynomial for groups with torsion\n");
      return (0);
    }
  }
  rank = p->gennum - matrixrank;
  if (rank > 2)
  {
    assert (matrixrank > 0 || numcols > 0);
    //if (matrixrank == 0 && numcols == 0)
    //{
    //  /* special case, trivial polynomial */
    //  if (!quiet) printf ("Alexander polynomial (special rank=%d case):\n", rank);
    //  printf ("1\n");
    //  return (1);
    //}
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

  if (verbose) printf ("Matrix has %d rows and %d columns\n", matrixrank, numcols);

  if (matrixrank < numcols && deficiency != 1)  /* deficiency == 1: link */
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander() {}\n");
    printf ("Cannot compute Alexander polynomial because there are too many relators\n");
    return (0);
  }

  if (matrixrank > numcols)
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
      laurent_canonify (determinant);
      printout_ideal1 (0, determinant);
      break;

      case 2:
      case 3:
      case 4:
      ai = laurent_notfirst_elementary_ideal (p, gconj, foxd - 1);
      if (ai->l1num > 1) printf ("# *** Warning: result can be noncanonical ***\n");
      alexander_fromideal (ai);
      break;

      default:
      foxdtoolarge++;
      break;
    }
    break;

    case 2:
    if (foxd > deficiency)
    {
      foxdtoolarge++;
      break;
    }
    assert (foxd == deficiency);
    determinant2 = laurent_eliminate_two_indeterminates (p, gconj2, gconj, &extradeterminants);
    extradets = 1;
    if (deficiency == 2) extradets = numcols;
    if (extradeterminants == 0) extradets = 0;
    if (verbose)
    {
      printf ("Alexander ideal before canonization:\n");
      printf ("-------------\n");
      printf ("Main polynomial:\n");
      print_laurentpoly2 (determinant2, 'u', 'v');
      printf ("\n");
      for (j = 0; j < extradets; j++)
      {
        printf ("Fundamental ideal factor %d:\n", j);
        print_laurentpoly2 (extradeterminants[j], 'u', 'v');
        printf ("\n");
      }
      printf ("-------------\n");
    }
    if (shuffle)
    {
      shuffle_poly2 (&determinant2, extradeterminants, extradets);
    }
    if (canonify_ideal2 (&determinant2, extradeterminants, extradets) == 0)
      printf ("# *** Warning: result can be noncanonical ***\n");
    if (extradeterminants)
    {
      if (deficiency == 2) laurent_canonify2 (determinant2);
      for (j = 0; j < extradets; j++) laurent_canonify2 (extradeterminants[j]);
    } else laurent_canonify2 (determinant2);
    printout_ideal2 (0, determinant2, extradeterminants, extradets, (deficiency == 2));
    break;

    default:
    printf ("Software unable to compute Alexander polynomial with these many indeterminates.\n");
    return (0);
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
  struct laurentpoly2 *determinant2, **extradeterminants;

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
        for (i = 0; i < ai->l1num; i++)
        {
          print_laurentpoly (ai->l1[i], 't');
          printf (";\n");
        }
        printf ("}\n");
      }
      for (i = 0; i < ai->l1num; i++) laurent_canonify (ai->l1[i]);
      printout_ideal1 (ai, 0);
    break;

    case 2:
      if (verbose)
      {
        printf ("Alexander ideal before canonization:\n{\n");
        for (i = 0; i < ai->l2num; i++)
        {
          print_laurentpoly2 (ai->l2[i], 'u', 'v');
          printf (";\n");
        }
        for (i = 0; i < ai->fl2num; i++)
        {
          printf ("F: ");
          print_laurentpoly2 (ai->fl2[i], 'u', 'v');
          printf (";\n");
        }
        printf ("}\n");
      }
      if (ai->l2num > 1) /* cannot use canonify_ideal2 */
      {
        for (i = 0; i < ai->l2num; i++) laurent_canonify2 (ai->l2[i]);
        for (i = 0; i < ai->fl2num; i++) laurent_canonify2 (ai->fl2[i]);
      } else {
        determinant2 = 0;
        assert (ai->l2num <= 1);
        if (ai->l2num == 1) determinant2 = ai->l2[0];
        extradets = ai->fl2num;
        extradeterminants = 0;
        if (extradets > 0)
        {
          extradeterminants = (struct laurentpoly2 **) malloc (extradets*sizeof (struct laurentpoly2 *));
          for (i = 0; i < extradets; i++) extradeterminants[i] = ai->fl2[i];
        }
        if (canonify_ideal2 (&determinant2, extradeterminants, extradets) == 0)
          printf ("# *** Warning: result can be noncanonical ***\n");
        if (ai->l2num == 1) ai->l2[0] = determinant2;
        if (extradets > 0)
        {
          for (i = 0; i < extradets; i++) ai->fl2[i] = extradeterminants[i];
          free (extradeterminants);
        }
        printout_ideal2 (ai, 0, 0, 0, 0);
      }
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
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("ideal(u,v) {\n");
    if (!quiet) printf ("Alexander ideal generated by:\n");
    if (quiet && outformat != OUTFORMAT_APPCONTOUR) printf ("Ideal [\n");
    if (principal) {print_laurentpoly (principal, 't'); printf (";\n");}
    if (ai)
    {
      for (j = 0; j < ai->l1num; j++)
      {
        print_laurentpoly (ai->l1[j], 't');
        printf (";\n");
      }
    }
    if (quiet && outformat != OUTFORMAT_APPCONTOUR) printf ("]\n");
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
    if (ai == 0) print_laurentpoly (principal, 't');
     else {
      assert (ai->l1num == 1);
      print_laurentpoly (ai->l1[0], 't');
    }
    printf (";\n");
  }
  if (outformat == OUTFORMAT_APPCONTOUR) printf ("}\n");
}

/*
 * printout a two-indets ideal
 */

void
printout_ideal2 (struct alexanderideal *ai, struct laurentpoly2 *principal, 
                     struct laurentpoly2 **fundamentals, int fnum, int printprincipal)
{
  extern int outformat, quiet;
  int j;

  if (fundamentals == 0) assert (fnum == 0);
  if (ai) assert (ai->indets == 2);
  if (ai) assert (principal == 0 && fundamentals == 0);

  if (fnum > 0 || (ai && ai->fl2num > 0) || (ai && ai->l2num >= 2))
  {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("ideal(u,v) {\n");
    if (!quiet) printf ("Alexander ideal generated by:\n");
    if (quiet && outformat != OUTFORMAT_APPCONTOUR) printf ("Ideal [\n");
    if (printprincipal || principal || (ai && ai->l2num > 0))
    {
      if (principal || printprincipal) print_laurentpoly2 (principal, 'u', 'v');
      printf (";\n");
      if (ai)
      {
        for (j = 0; j < ai->l2num; j++)
        {
          print_laurentpoly2 (ai->l2[j], 'u', 'v');
          printf (";\n");
        }
      }
    }
    for (j = 0; j < fnum; j++)
    {
      if (fundamentals[j] == 0) continue;
      if (outformat == OUTFORMAT_APPCONTOUR)
      {
        printf ("F: ");
        print_laurentpoly2 (fundamentals[j], 'u', 'v');
        printf (";\n");
      } else {
        printf ("(");
        print_laurentpoly2 (fundamentals[j], 'u', 'v');
        printf (") (u - 1);\n");
        printf ("(");
        print_laurentpoly2 (fundamentals[j], 'u', 'v');
        printf (") (v - 1);\n");
      }
    }
    if (ai)
    {
      for (j = 0; j < ai->fl2num; j++)
      {
        if (ai->fl2[j] == 0) continue;
        if (outformat == OUTFORMAT_APPCONTOUR)
        {
          printf ("F: ");
          print_laurentpoly2 (ai->fl2[j], 'u', 'v');
          printf (";\n");
        } else {
          printf ("(");
          print_laurentpoly2 (ai->fl2[j], 'u', 'v');
          printf (") (u - 1);\n");
          printf ("(");
          print_laurentpoly2 (ai->fl2[j], 'u', 'v');
          printf (") (v - 1);\n");
        }
      }
    }
    if (quiet && outformat != OUTFORMAT_APPCONTOUR) printf ("]\n");
  } else {
    if (outformat == OUTFORMAT_APPCONTOUR) printf ("alexander(u,v) {\n");
    if (!quiet) printf ("Alexander polynomial:\n");
    if (ai == 0) print_laurentpoly2 (principal, 'u', 'v');
     else {
      assert (ai->l2num == 1);
      print_laurentpoly2 (ai->l2[0], 'u', 'v');
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
  int i, sum, rank, matrixrank;
  int deriv, val, e;
  int numcols = 0;
  int gconj;

  topreabelian (p);

  if (verbose) print_presentation (p);

  for (r = p->rules; r; r = r->next) numcols++;
  rank = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    assert (sum >= 0);
    if (sum) matrixrank = i;
    if (sum && sum != 1)
    {
      printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
      return (0);
    }
  }
  rank = p->gennum - matrixrank;
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

  //laurent_canonify (determinant);

  if (verbose)
  {
    printf ("Determinant of matrix: ");
    print_laurentpoly (determinant, 't');
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
      deriv += e*determinant->stem[i];
      val += determinant->stem[i];
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
  int i, sum, rank, matrixrank;
  int deriv, val, e;
  int numcols = 0;
  int gconj;

  topreabelian (p);

  if (verbose) print_presentation (p);

  for (r = p->rules; r; r = r->next) numcols++;
  rank = 0;
  matrixrank = 0;
  for (i = 1, r = p->rules; r && i <= p->gennum; i++, r = r->next)
  {
    sum = get_exp_sum (r, i);
    assert (sum >= 0);
    if (sum) matrixrank = i;
    if (sum && sum != 1)
    {
      printf ("Cannot compute corank one Alexander polynomial for groups with torsion\n");
      return (0);
    }
  }
  rank = p->gennum - matrixrank;
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

  laurent_canonify (determinant);

  if (verbose)
  {
    printf ("Determinant of matrix (before canonization): ");
    print_laurentpoly (determinant, 't');
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
      deriv += e*determinant->stem[i];
      val += determinant->stem[i];
    }
  }
  printf ("Value in 1: %d, deriv = %d\n", val, deriv);
  if (!quiet) printf ("Alexander polynomial (up to t -> 1/t):\n");
  print_laurentpoly (determinant, 't');
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

  matrix = laurent_build_matrix (p, eliminate);

  assert (matrix->numcols <= matrix->numrows);
  if (matrix->numcols < matrix->numrows) return (0);

  determinant = laurent_compute_determinant (matrix->columns, matrix->numcols);

  laurent_free_matrix (matrix);

  return (determinant);
}

/*
 * second elementary ideal
 */

struct alexanderideal *
laurent_notfirst_elementary_ideal (struct presentation *p, int eliminate, int corank)
{
  extern int verbose;
  struct alexanderideal *ai;
  struct laurentmatrix *matrix, *minor;
  struct laurentpoly *l, **matrixcolumn;
  int rank, i, j, jj, idx;

  matrix = laurent_build_matrix (p, eliminate);

  assert (matrix->numcols <= matrix->numrows);
  if (matrix->numcols < matrix->numrows - corank) return (0);

  rank = matrix->numrows - corank;
  ai = (struct alexanderideal *) malloc (sizeof (struct alexanderideal));
  ai->indets = 1;
  ai->spread = 1;
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
    ai->l1num = matrix->numrows;
    if (matrix->numcols == matrix->numrows) ai->l1num *= ai->l1num;
    if (ai->l1num > IDEAL_MAX_GENERATORS_NUM)
    {
      printf ("Fatal: too many generators (%d) for the ideal\n", ai->l1num);
      laurent_free_matrix (matrix);
      free (ai);
      return (0);
    }
    idx = 0;
    for (i = 0; i < matrix->numrows; i++)
    {
      for (j = 0; j < matrix->numcols; j++)
      {
        matrixcolumn = matrix->columns[j];
        l = matrixcolumn[i];
        assert (idx < ai->l1num);
        ai->l1[idx++] = laurent_dup(l);
      }
    }
    while (idx < ai->l1num) ai->l1[idx++] = 0;
    break;

    default:
    assert (corank >= 1);
    if (corank > 1)
    {
      printf ("Corank larger than 1 is not yet implemented\n");
      free (ai);
      laurent_free_matrix (matrix);
      return (0);
    }
    /* this is the corank 1 case, with rank larger than 1 */
    assert (matrix->numrows == matrix->numcols);
    idx = 0;
    for (i = 0; i < matrix->numrows; i++)
    {
      for (j = 0; j < matrix->numcols; j++)
      {
        if (idx >= IDEAL_MAX_GENERATORS_NUM)
        {
          printf ("Fatal: too many generators (%d) for the ideal\n", idx+1);
          laurent_free_matrix (matrix);
          free (ai);
          return (0);
        }
        assert (idx < IDEAL_MAX_GENERATORS_NUM);
        minor = minor_matrix_corank1 (matrix, i, j);
        ai->l1[idx++] = laurent_compute_determinant (minor->columns, minor->numcols);
        for (jj = 0; jj < minor->numcols; jj++) free (minor->columns[jj]);
        free (minor->columns);
        free (minor);
        ai->l1num = idx;
      }
    }
    break;
  }

  laurent_free_matrix (matrix);
  laurent_simplify_ideal (ai);
  return (ai);
}

/*
 * build Alexander matrix (rank 1: e.g. knots)
 */

struct laurentmatrix *
laurent_build_matrix (struct presentation *p, int eliminate)
{
  struct presentationrule *r;
  struct laurentmatrix *matrix;
  struct laurentpoly *l, **matrixcolumn;
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
      matrixcolumn[ii++] = laurent_get_exp_sum (r, i, eliminate);
    }
  }

  if (debug)
  {
    for (j = 0, r = p->rules; r; j++, r = r->next)
    {
      printf ("extra row, entry %d: ", j+1);
      print_laurentpoly (laurent_get_exp_sum (r, eliminate, eliminate), 't');
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
        print_laurentpoly (l, 't');
        printf ("; \t");
      }
      printf ("\n");
    }
  }
  return (matrix);
}

struct laurentmatrix *
minor_matrix_corank1 (struct laurentmatrix *matrix, int dropi, int dropj)
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
 * free allocated space for Alexander matrix
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
      if (l) free (l);
    }
    free (matrixcolumn);
  }
  free (matrix->columns);
  free (matrix);
}

void
laurent_simplify_ideal (struct alexanderideal *ai)
{
  struct laurentpoly *oldgcd, *newgcd;
  extern int principal, verbose;
  int last, i, spread, lspread;
  int loop = 1;

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
        printf ("WARNING: inhibiting gcd computation due to integer size\n");
        break;
      }
      newgcd = laurent_gcd (spread, ai->l1[last], ai->l1[last-1], &lspread);
      spread = lspread;
      for (i = last-2; i >= 0; i--)
      {
        oldgcd = newgcd;
        newgcd = laurent_gcd (spread, oldgcd, ai->l1[i], &lspread);
        spread = lspread;
        free (oldgcd);
      }
      if (newgcd->stemdegree > ai->l1[0]->stemdegree ||
          (newgcd->stemdegree == ai->l1[0]->stemdegree && abs(spread*newgcd->stem[0]) >= abs(ai->l1[0]->stem[0])) )
      {
        free (newgcd);
      } else {
        if (verbose) printf ("Information: degree reduction in ideal generators!\n");
        for (i = 0; i <= newgcd->stemdegree; i++) newgcd->stem[i] *= spread;
        if (ai->l1num < IDEAL_MAX_GENERATORS_NUM)
        {
          for (i = ai->l1num; i > 0; i--) ai->l1[i] = ai->l1[i-1];
          ai->l1num++;
          ai->l1[0] = newgcd;
          loop = 1;
        } else {
          printf ("Warning: no space left to add new generator!\n");
          free (newgcd);
        }
      }
    }
  }

  for (i = 0; i < ai->l1num; i++)
    laurent_canonifysign (ai->l1[i]);

  if (principal && ai->l1num > 1)
  {
    spread = 1;
    for (i = 1; i < ai->l1num; i++)
    {
      newgcd = laurent_gcd (spread, ai->l1[0], ai->l1[i], &lspread);
      spread = lspread;
      free (ai->l1[0]);
      free (ai->l1[i]);
      ai->l1[0] = newgcd;
    }
    ai->l1num = 1;
    ai->spread = spread;
    return;
  }

  return;  /* TODO: not implemented at the moment! */
}

//#define USE_EUCLID_GCD 1

/*
 * try to simplify the ideal by removal of some generator or
 * reduction in complexity of one generator
 * no increase is allowed
 * return 1 if at least one simplification was performed
 */

int laurent_try_reduce_pair (struct laurentpoly *l1, struct laurentpoly *l2);
void laurent_sort_entries (int num, struct laurentpoly *l[]);

int
laurent_try_simplify_ideal (struct alexanderideal *ai)
{
  struct laurentpoly *l, *l1, *l2;
#ifdef USE_EUCLID_GCD
  struct laurentpoly *newgcd;
#endif
  int i, j, spread, numreductions;
  extern int verbose;

  numreductions = spread = 0;
  if (ai->indets != 1 ) return (0);

  if (ai->l1num <= 1) return (0); // at least two generators...

  laurent_sort_entries (ai->l1num, ai->l1);
  for (i = 0; i < ai->l1num; i++)
  {
    if (ai->l1[i] == 0)
    {
      for (j = i; j < ai->l1num - 1; j++) ai->l1[j] = ai->l1[j+1];
      ai->l1num--;
      if (verbose) printf ("Ideal simplification: zero polynomial removed\n");
      return (1);
    }
  }
  for (i = 0; i < ai->l1num; i++)
  {
    l = ai->l1[i];
    if (l->stemdegree == 0 && abs(l->stem[0]) == 1)
    {
      for (j = 0; j < ai->l1num; j++)
      {
        if (j != i && ai->l1[j]) free (ai->l1[j]);
        ai->l1[j] = 0;
      }
      ai->l1[0] = l;
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
      newgcd = laurent_gcd (1, ai->l1[i], ai->l1[j], &spread);
      if (spread == 1)
      {
        if (ai->l1[i]) free (ai->l1[i]);
        if (ai->l1[j]) free (ai->l1[j]);
        ai->l1[i] = newgcd;
        ai->l1[j] = 0;
        if (verbose) printf ("Ideal simplification: two pols generate a principal ideal\n");
        assert (newgcd);
        assert (newgcd->stem[0]);
        assert (newgcd->stem[newgcd->stemdegree]);
        return (1);
      } else free (newgcd);
    }
#endif
  }
  for (i = 0; i < ai->l1num; i++)
  {
    l1 = ai->l1[i];
    if (l1 == 0) continue;
    assert (l1->stem[0]);
    assert (l1->stem[l1->stemdegree]);
    for (j = i+1; j < ai->l1num; j++)
    {
      /*
       * kind of euclidean division step
       */
      l2 = ai->l1[j];
      if (l2 == 0) continue;
      assert (l2->stem[0]);
      assert (l2->stem[l2->stemdegree]);
      if (l1->stemdegree <= l2->stemdegree) numreductions += laurent_try_reduce_pair (l2, l1);
      if (l2->stemdegree == 0 && l2->stem[0] == 0)
      {
        free (l2);
        ai->l1[j] = 0;
        continue;
      }
      if (l1->stemdegree >= l2->stemdegree) numreductions += laurent_try_reduce_pair (l1, l2);
      if (l1->stemdegree == 0 && l1->stem[0] == 0)
      {
        free (l1);
        ai->l1[i] = 0;
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
 */

int
laurent_try_reduce_pair (struct laurentpoly *l1, struct laurentpoly *l2)
{
  int f, i, j, k, c1, c2first, c2last, quotient;
  int deltadegree;
  int s2first = 1;
  int s2last = 1;
  int numreductions = 0;
  int maxc2size, safesize;

  assert (l1);
  assert (l2);
  assert (l1->stemdegree >= l2->stemdegree);

  c1 = l1->stem[0];
  c2first = l2->stem[0];
  assert (c2first);
  maxc2size = 0;
  for (i = 0; i <= l2->stemdegree; i++) if (abs(l2->stem[i]) > maxc2size) maxc2size = abs(l2->stem[i]);

  quotient = c1/c2first;
  safesize = INT_MAX;
  if (quotient) safesize /= abs(quotient);
  safesize /= 4;
  if (maxc2size > safesize)
  {
    if (int_overflow_encountered++ == 0)
      printf ("WARNING: above max allowed int size, cannot complete simplification\n");
    return (0);
  }

  if (c1 == c2first*quotient)
  {
    for (i = 0; i <= l2->stemdegree; i++)
    {
      l1->stem[i] -= quotient*l2->stem[i];
    }
    assert (l1->stem[0] == 0);
    while (l1->stem[0] == 0 && l1->stemdegree > 0)
    {
      for (i = 0; i < l1->stemdegree; i++) l1->stem[i] = l1->stem[i+1];
      l1->stemdegree--;
      l1->minexpon++;
    }
    while (l1->stem[l1->stemdegree] == 0 && l1->stemdegree > 0)
      l1->stemdegree--;
    if (l1->stemdegree > 0)
    {
      assert (l1->stem[0]);
      assert (l1->stem[l1->stemdegree]);
    }
    return (1);
  }
  c1 = l1->stem[l1->stemdegree];
  c2last = l2->stem[l2->stemdegree];
  assert (c2last);
  quotient = c1/c2last;
  if (c1 == c2last*quotient)
  {
    k = l1->stemdegree - l2->stemdegree;
    for (i = 0; i <= l2->stemdegree; i++)
    {
      l1->stem[i + k] -= quotient*l2->stem[i];
    }
    assert (l1->stem[l1->stemdegree] == 0);
    while (l1->stem[l1->stemdegree] == 0 && l1->stemdegree > 0)
      l1->stemdegree--;
    while (l1->stem[0] == 0 && l1->stemdegree > 0)
    {
      for (i = 0; i < l1->stemdegree; i++) l1->stem[i] = l1->stem[i+1];
      l1->stemdegree--;
      l1->minexpon++;
    }
    if (l1->stemdegree > 0)
    {
      assert (l1->stem[0]);
      assert (l1->stem[l1->stemdegree] != 0);
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
    if (l1->stem[k] > c2first/2)
    {
      numreductions++;
      f = (l1->stem[k] - c2first/2 + c2first - 1)/c2first;
      for (i = 0, j = k; i <= l2->stemdegree; i++)
      {
        l1->stem[j++] -= f*s2first*l2->stem[i];
      }
    }
    if (l1->stem[k] < -(c2first-1)/2)
    {
      numreductions++;
      f = (-(c2first-1)/2 - l1->stem[k] + c2first - 1)/c2first;
      for (i = 0, j = k; i <= l2->stemdegree; i++)
      {
        l1->stem[j++] += f*s2first*l2->stem[i];
      }
    }
    if (k >= (deltadegree+1)/2) continue;
    if (l1->stem[l1->stemdegree - k] > c2last/2)
    {
      numreductions++;
      f = (l1->stem[l1->stemdegree - k] - c2last/2 + c2last - 1)/c2last;
      for (i = 0, j = deltadegree - k; i <= l2->stemdegree; i++)
      {
        l1->stem[j++] -= f*s2last*l2->stem[i];
      }
    }
    if (l1->stem[l1->stemdegree - k] < -(c2last-1)/2)
    {
      numreductions++;
      f = (-(c2last-1)/2 - l1->stem[l1->stemdegree - k] + c2last - 1)/c2last;
      for (i = 0, j = deltadegree - k; i <= l2->stemdegree; i++)
      {
        l1->stem[j++] += f*s2last*l2->stem[i];
      }
    }
  }

  return (numreductions);
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
      if (abs (li->stem[s]) < abs (lj->stem[s])) break;
      if (abs (li->stem[s]) > abs (lj->stem[s])) {iwins = 0; break;}
      if (li->stem[s] < lj->stem[s]) break;
      if (li->stem[s] > lj->stem[s]) {iwins = 0; break;}
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

struct laurentpoly2 *
laurent_eliminate_two_indeterminates (struct presentation *p, int e1, int e2,
             struct laurentpoly2 ***extradeterminantspt)
{
  extern int verbose, debug, quiet;
  struct presentationrule *r;
  struct laurentpoly2 *l, *determinant;
  struct laurentpoly2 **matrixcolumn;
  struct laurentpoly2 ***matrix;
  struct laurentpoly2 **extradeterminants;
  int numcols = 0;
  int numrows, deficiency;
  int i, ii, j;
  int extrafound = 0;

  assert (e1 >= 1 && e2 >= 1);
  assert (e1 <= p->gennum);
  assert (e2 <= p->gennum);
  assert (e1 != e2);

  for (r = p->rules; r; r = r->next) numcols++;
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
   * $ echo "fpgroup{<a,b; aA>}" | contour --out alexander
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

  numrows = p->gennum - 1;
  /*
   * last row will contain the common extra factor
   */

  matrix = (struct laurentpoly2 ***) malloc (numcols*sizeof(struct laurentpoly2 **));
  for (j = 0; j < numcols; j++)
  {
    matrix[j] = (struct laurentpoly2 **) malloc (numrows*sizeof(struct laurentpoly2 *));
  }

  /*
   * fill matrix
   */

  for (j = 0, r = p->rules; r; j++, r = r->next)
  {
    matrixcolumn = matrix[j];
    for (i = 1, ii = 0; i <= p->gennum; i++)
    {
      if (i == e1 || i == e2) continue;
      matrixcolumn[ii++] = laurent_get_exp_sum2 (r, i, e1, e2);
    }
    assert (ii == numrows - 1);
    matrixcolumn[ii] = laurent_common_factor2 (r, e1, e2);
    if (matrixcolumn[ii]) extrafound++;
  }

  if (debug)
  {
    for (j = 1, r = p->rules; r; j++, r = r->next)
    {
      printf ("Extra row 1, entry %d: ", j);
      print_laurentpoly2 (laurent_get_exp_sum2 (r, e1, e1, e2), 'u', 'v');
      printf (";\n");
      printf ("Extra row 2, entry %d: ", j);
      print_laurentpoly2 (laurent_get_exp_sum2 (r, e2, e1, e2), 'u', 'v');
      printf (";\n");
      printf ("Fox mixed derivative, entry %d: ", j);
      print_laurentpoly2 (laurent_mixed_derivative2 (r, e1, e2), 'u', 'v');
      printf (";\n");
      printf ("Common factor, entry %d: ", j);
      print_laurentpoly2 (matrix[j-1][numrows-1], 'u', 'v');
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
        matrixcolumn = matrix[j];
        l = matrixcolumn[i];
        print_laurentpoly2 (l, 'u', 'v');
        printf ("; \t");
      }
      printf ("\n");
    }
    printf ("The last row contains the commutator factors.\n");
  }

  /* determinant of the square matrix in case of links, of the principal minor otherwise */
  determinant = laurent_compute_determinant2 (matrix, numcols);

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
             (struct laurentpoly2 **) malloc (sizeof(struct laurentpoly2 *));
        extradeterminants[0] = determinant;
        determinant = 0;
      }
    }
    break;

    case 2:  /* bouquet of two loops */
    assert (numrows == numcols + 1);
    if (extrafound && extradeterminantspt == 0)
    {
      printf ("FATAL: extra commutator factors found but computation forbidden!\n");
      printf ("  printed result will be incorrect\n");
    }
    if (extrafound && extradeterminantspt)
    {
      *extradeterminantspt = extradeterminants =
           (struct laurentpoly2 **) malloc (numcols*sizeof(struct laurentpoly2 *));
      for (j = 0; j < numcols; j++)
      {
        extradeterminants[j] = laurent_minor_determinant2 (matrix, numcols, j);
      }
    }
    if (extrafound == 0 && extradeterminantspt) *extradeterminantspt = 0;
    break;

    default:
    printf ("FATAL: invalid deficiency: %d\n", deficiency);
  }

  /*
   * free allocated space
   */

  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    for (i = 0; i < numrows; i++)
    {
      l = matrixcolumn[i];
      if (l) free_laurentpoly2 (l);
    }
    free (matrixcolumn);
  }
  free (matrix);

  return (determinant);
}

/*
 * Compute matrix entries for Alexander polynomial computation
 */

struct laurentpoly *
laurent_get_exp_sum (struct presentationrule *r, int n, int gconj)
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
  l = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + (stemdegree + 1)*sizeof(int));
  l->stemdegree = stemdegree;
  l->minexpon = minexp;
  for (d = 0; d <= stemdegree; d++) l->stem[d] = 0;

  runningexp = 0;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -gconj) runningexp--;
    if (abs(r->var[k]) == n)
    {
      d = runningexp - minexp;
      if (r->var[k] > 0) l->stem[d]++;
        else l->stem[d]--;
    }
    if (r->var[k] == gconj) runningexp++;
  }

  assert (l);
  if (laurent_normalize(l) == 0)
  {
    free (l);
    return (0);
  }
  assert (l && l->stem[0]);
  return (l);
}

/*
 * Compute matrix entries for Alexander polynomial computation (two indeterminates)
 */

struct laurentpoly2 *
laurent_get_exp_sum2 (struct presentationrule *r, int n, int e1, int e2)
{
  int k, d;
  int stemdegree;
  int re1, re2, minexp1, minexp2, maxexp1, maxexp2;
  struct laurentpoly2 *l;

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
  l = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2)
            + (stemdegree + 1)*sizeof(struct laurentpoly *));
  l->stemdegree = stemdegree;
  l->minexpon = minexp2;
  for (d = 0; d <= stemdegree; d++) l->stem[d] = 0;

  re1 = re2 = 0;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -e1) re1--;
    if (r->var[k] == -e2) re2--;
    if (abs(r->var[k]) == n)
    {
      d = re2 - minexp2;
      if (r->var[k] > 0) l->stem[d] = laurentpoly_addmonom (l->stem[d], re1, 1);
        else l->stem[d] = laurentpoly_addmonom (l->stem[d], re1, -1);
    }
    if (r->var[k] == e1) re1++;
    if (r->var[k] == e2) re2++;
  }

  assert (l);
  if (laurent_normalize2(l) == 0)
  {
    free_laurentpoly2 (l);
    return (0);
  }
  assert (l && l->stem[0]);
  return (l);
}

/*
 * Compute matrix entries on pairs of indeterminates using Fox
 * second derivatives, assuming expsum is zero for both
 * It is \partial_{x2}\partial_{x1} w (1,u,v)
 * where u is substituted in place of x1 and v in place of x2,
 * 1 in all others
 */

struct laurentpoly2 *
laurent_mixed_derivative2 (struct presentationrule *r, int x1, int x2)
{
  int k, d;
  int stemdegree;
  int re1, re2, minexp1, minexp2, maxexp1, maxexp2;
  struct laurentpoly2 *l;

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
  l = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2)
            + (stemdegree + 1)*sizeof(struct laurentpoly *));
  l->stemdegree = stemdegree;
  l->minexpon = minexp2;
  for (d = 0; d <= stemdegree; d++) l->stem[d] = 0;

  re1 = 0;
  re2 = -minexp2;
  for (k = 0; k < r->length; k++)
  {
    if (r->var[k] == -x1) re1--;
    if (r->var[k] == -x2) re2--;
    if (r->var[k] == x2) l->stem[re2] = laurentpoly_addmonom (l->stem[re2], re1, -re1);
    if (r->var[k] == -x2) l->stem[re2] = laurentpoly_addmonom (l->stem[re2], re1, re1);
    //if (abs(r->var[k]) == x2)
    //{
    //  d = re2 - minexp2;
    //  coef = -re1;
    //  if (r->var[k] < 0) coef = re1;
    //  l->stem[d] = laurentpoly_addmonom (l->stem[d], re1, coef);
    //}
    if (r->var[k] == x1) re1++;
    if (r->var[k] == x2) re2++;
  }

  assert (l);
  if (laurent_normalize2(l) == 0)
  {
    free_laurentpoly2 (l);
    return (0);
  }
  assert (l && l->stem[0]);
  return (l);
}

/*
 * Compute common factor on pairs of indeterminates
 * expsum is zero for both
 */

struct laurentpoly2 *
laurent_common_factor2 (struct presentationrule *r, int x1, int x2)
{
  int k;
  struct laurentpoly2 *l, *poly1, *poly2;
  struct laurentpoly *res, *addres, *mulres, *uminus1;

  poly1 = laurent_get_exp_sum2 (r, x1, x1, x2);
  /* this should be divisible by (x2 - 1) */
  res = laurent_sum_coefficients2 (poly1);
  assert (res == 0);

  /* now divide by x2-1 */
  if (poly1)
  {
    for (k = poly1->stemdegree - 1; k >= 0; k--)
    {
      addres = laurent_add (poly1->stem[k+1], poly1->stem[k]);
      if (poly1->stem[k]) free(poly1->stem[k]);
      poly1->stem[k] = addres;
    }
    assert (poly1->stem[0] == 0);
    poly1->minexpon--;
    assert (laurent_normalize2(poly1));
  }

  poly2 = laurent_get_exp_sum2 (r, x2, x1, x2);
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
    uminus1 = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + 2*sizeof(int));
    uminus1->minexpon = 0;
    uminus1->stemdegree = 1;
    uminus1->stem[0] = -1;
    uminus1->stem[1] = 1;
    for (k = 0; k <= poly1->stemdegree; k++)
    {
      mulres = laurent_mul (poly1->stem[k], uminus1);
      addres = laurent_add (poly2->stem[k], mulres);
      if (mulres) free (mulres);
      assert (addres == 0);
    }
    free (uminus1);
  }

  l = poly1;

  free_laurentpoly2 (poly2);
  return (l);
}

/*
 * compute the determinant of the matrix
 */

struct laurentpoly *
laurent_compute_determinant (struct laurentpoly ***matrix, int n)
{
  int i, ii, jj, i1;
  int sign;
  struct laurentpoly *determinant = 0, *subdeterminant, *product, *sum;
  struct laurentpoly ***submatrix;
  struct laurentpoly **matrixcol, **submatrixcol, **firstcolumn;

  assert (n >= 0);
  if (n == 0)
  {
    determinant = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + sizeof (int));
    determinant->minexpon = 0;
    determinant->stemdegree = 0;
    determinant->stem[0] = 1;
    return (determinant);
  }
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

    subdeterminant = laurent_compute_determinant (submatrix, n-1);
    product = laurent_mul (subdeterminant, firstcolumn[i]);
    if (subdeterminant) free (subdeterminant);
    if (sign < 0) laurent_negate (product);
    sum = laurent_add (determinant, product);
    if (product && product != sum) free (product);
    if (determinant && sum != determinant) free (determinant);
    determinant = sum;
  }

  for (jj = 0; jj < n-1; jj++) free (submatrix[jj]);
  free (submatrix);

  return (determinant);
}

/*
 * compute the determinant of the matrix (two indeterminates)
 */

struct laurentpoly2 *
laurent_compute_determinant2 (struct laurentpoly2 ***matrix, int n)
{
  int i, ii, jj, i1;
  int sign;
  struct laurentpoly2 *determinant = 0, *subdeterminant, *product, *sum;
  struct laurentpoly2 ***submatrix;
  struct laurentpoly2 **matrixcol, **submatrixcol, **firstcolumn;
  struct laurentpoly *lp1;

  assert (n >= 0);
  if (n == 0)
  {
    determinant = (struct laurentpoly2 *) malloc (sizeof (struct laurentpoly2) + sizeof (struct laurentpoly *));
    determinant->minexpon = 0;
    determinant->stemdegree = 0;
    lp1 = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + sizeof (int));
    determinant->stem[0] = lp1;
    lp1->minexpon = 0;
    lp1->stemdegree = 0;
    lp1->stem[0] = 1;
    return (determinant);
  }
  if (n == 1) return (laurent_dup2(matrix[0][0]));

  /*
   * developping the determinant about the first column
   */

  submatrix = (struct laurentpoly2 ***) malloc ( (n-1)*sizeof (struct laurentpoly2 **) );
  for (jj = 0; jj < n-1; jj++)
    submatrix[jj] = (struct laurentpoly2 **) malloc ( (n-1)*sizeof (struct laurentpoly2 *) );

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

    subdeterminant = laurent_compute_determinant2 (submatrix, n-1);
    product = laurent_mul2 (subdeterminant, firstcolumn[i]);
    if (subdeterminant) free_laurentpoly2 (subdeterminant);
    if (sign < 0) laurent_negate2 (product);
    sum = laurent_add2 (determinant, product);
    if (product && product != sum) free_laurentpoly2 (product);
    if (determinant && sum != determinant) free_laurentpoly2 (determinant);
    determinant = sum;
  }

  for (jj = 0; jj < n-1; jj++) free (submatrix[jj]);
  free (submatrix);

  return (determinant);
}

/*
 * compute the determinant of the matrix (two indeterminates)
 * with row 
 */

struct laurentpoly2 *
laurent_minor_determinant2 (struct laurentpoly2 ***matrix, int numcols, int row_to_subst)
{
  int j;
  struct laurentpoly2 **matrixcolumn;
  struct laurentpoly2 *saved, *determinant;

  /* exchange rows */
  for (j = 0; j < numcols; j++)
  {
    matrixcolumn = matrix[j];
    saved = matrixcolumn[row_to_subst];
    matrixcolumn[row_to_subst] = matrixcolumn[numcols];
    matrixcolumn[numcols] = saved;
  }

  determinant = laurent_compute_determinant2 (matrix, numcols);

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
canonify_ideal2 (struct laurentpoly2 **wpt, struct laurentpoly2 **wi, int winum)
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
base_canonify2 (struct laurentpoly2 **ppt)
{
  extern int nobasecanonify;
  int suppdim;

  suppdim = laurent_suppdim2 (*ppt);
 
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
base_canonify2_onedim (struct laurentpoly2 **ppt)
{
  extern int verbose;
  struct laurentpoly2 *p;
  struct laurentpoly *p1;
  struct laurentpoly *newp1;
  int origy, dx, dy, k, xk;
  int mygcd, xstep, ystep;

  p = *ppt;
  assert (p);
  if (p->stemdegree == 0)
  {
    laurent_canonify (p->stem[0]);
    laurent_canonifysign (p->stem[0]);
    return (1);  // already canonical!
  }

  p1 = p->stem[0];
  origy = p1->minexpon;  // orig u

  dx = p->stemdegree;    // dv
  assert (dx);
  p1 = p->stem[p->stemdegree];
  assert (p1);
  assert (p1->stemdegree == 0);
  dy = p1->minexpon - origy;  // du

  mygcd = gcd (dx, dy);  /* this is the degree of the resulting polynomial */
  xstep = dx/mygcd;
  ystep = dy/mygcd;
  assert (xstep > 0);

  assert (mygcd >= 1);

  newp1 = (struct laurentpoly *) malloc (sizeof (struct laurentpoly) + (mygcd+1)*sizeof(int));
  newp1->stemdegree = mygcd;

  for (k = 0; k <= mygcd; k++) newp1->stem[k] = 0;
  for (k = 0, xk = 0; k <= mygcd; k++, xk += xstep)
  {
    /* fill up new polynomial */
    p1 = p->stem[xk];
    if (p1 == 0) continue;
    assert (p1->minexpon == origy + k*ystep);
    assert (p1->stemdegree == 0);
    newp1->stem[k] = p1->stem[0];
  }
  newp1->minexpon = 0;

  /* free mem space of old poly */
  for (xk = 0; xk <= p->stemdegree; xk++)
  {
    if (xk % xstep) assert (p->stem[xk] == 0);
    if (p->stem[xk]) free (p->stem[xk]);
  }
  p->minexpon = 0;
  p->stemdegree = 0;
  p->stem[0] = newp1;
  laurent_canonifysign (p->stem[0]);

  return (1);
}

/*
 * canonification in the two-dim case
 */

int
base_canonify2_twodim (struct laurentpoly2 **ppt)
{
  struct laurentpoly2 *origp, *newp, *optp;
  struct laurentpoly2 *tempp;
  extern int verbose;
  int x1[2], x2[2], x3[2], du2, dv2, du3, dv3;
  int origtotdegree, newtotdegree, optdegree;
  int i, k;
  int xiu[3], xiv[3], sumxiu, sumxiv, amin, amax, bmin, bmax;
  int area2;  // area of the triangle with those vertices
  int a, b, c, d, bm[2][2];

  origp = *ppt;
  assert (origp);
  assert (origp->stem[0]);
  laurent_canonify2 (origp);
  optp = laurent_dup2 (origp);
  origtotdegree = laurent2_totdegree (origp);

  laurent_getthree2 (*ppt, x1, x2, x3);

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
            newp = laurent_dup2 (origp);
            newp = base_change2 (newp, bm);
            newtotdegree = laurent2_totdegree (newp);
            laurent_canonify2 (newp);
            if (newtotdegree > optdegree)
            {
              free_laurentpoly2 (newp);
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
              if (laurent2_lexicocompare (newp, optp) < 0)
              {
                if (verbose > 1) printf ("same totdegree, but better comparison\n");
                tempp = newp;
                newp = optp;
                optp = tempp;
              }
            }
            free_laurentpoly2 (newp);
            if (verbose > 1)
            {
              printf ("Found feasible matrix: [%d %d; %d %d]\n", a, b, c, d);
              printf ("New polynomial: ");
              print_laurentpoly2 (optp, 'u', 'v');
              printf (";\n");
            }
          }
        }
      }
    }
  }

  free_laurentpoly2 (origp);
  laurent2_canonifysign (optp);
  *ppt = optp;
  return (1);
}

/*
 * compute dimension of support
 */

int
laurent_suppdim2 (struct laurentpoly2 *p)
{
  struct laurentpoly *l1;
  int dirfound = 0;
  int k, origx, origy, dx, dy, ddx, ddy;

  if (p == 0) return (-1);  //empty support

  if (p->stemdegree == 0)
  {
    l1 = p->stem[0];
    if (l1->stemdegree == 0) return (0);
  }

  if (p->stemdegree == 0) return (1);

  /* let's see if dimension is one... */
  l1 = p->stem[0];
  assert (l1);
  if (l1->stemdegree > 0) return (2);
  origx = 0;
  origy = l1->minexpon;

  for (k = 1; k <= p->stemdegree; k++)
  {
    if (p->stem[k] == 0) continue;
    l1 = p->stem[k];
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
laurent_getthree2 (struct laurentpoly2 *p, int *x1, int *x2, int *x3)
{
  struct laurentpoly *l1, *l2;
  int dirfound = 0;
  int k, origv, origu, dv, du, ddv, ddu;

  assert (p != 0);
  assert (p->stemdegree > 0); // two different exponents for v

  l1 = p->stem[0];
  assert (l1);
  if (l1->stemdegree > 0)     // two different exponents for u
  {
    assert (l1->stem[0] && l1->stem[l1->stemdegree]);
    x1[0] = l1->minexpon;     // u exponent of first point
    x2[0] = l1->minexpon + l1->stemdegree;  // u exponent of second point
    x1[1] = x2[1] = p->minexpon;  // v exponent of first and second point
    l2 = p->stem[p->stemdegree];
    assert (l2);
    assert (l2->stem[0]);
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
    if (p->stem[k] == 0) continue;
    l1 = p->stem[k];
    if (l1->stemdegree > 0)
    {
      assert (l1->stem[0] && l1->stem[l1->stemdegree]);
      x1[0] = l1->minexpon;     // u exponent of first point
      x2[0] = l1->minexpon + l1->stemdegree;  // u exponent of second point
      x1[1] = x2[1] = p->minexpon + k;  // v exponent of first and second point
      l2 = p->stem[0];
      assert (l2->stem[0]);
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

struct laurentpoly2 *
base_change2 (struct laurentpoly2 *l, int matrixb[2][2])
{
  struct laurentpoly2 *newl = 0;
  struct laurentpoly *lu;
  int a, b, c, d, iu, iv, degu, degv, newdegu, newdegv;
  int coef;

  a = matrixb[0][0];
  b = matrixb[0][1];
  c = matrixb[1][0];
  d = matrixb[1][1];

  assert (isinvertible_base(matrixb));

  if (l == 0) return (0);

  /* loop through monomials of l */
  for (iv = 0; iv <= l->stemdegree; iv++)
  {
    degv = l->minexpon + iv;
    lu = l->stem[iv];
    if (lu == 0) continue;
    for (iu = 0; iu <= lu->stemdegree; iu++)
    {
      degu = lu->minexpon + iu;
      coef = lu->stem[iu];
      if (coef == 0) continue;
      newdegu = a*degu + b*degv;
      newdegv = c*degu + d*degv;
      newl = laurentpoly2_addmonom (newl, newdegu, newdegv, coef);
    }
  }

  free_laurentpoly2 (l);

  return (newl);
}

/*
 * randomly change base of Z^2 for the alexander ideal
 */

void
shuffle_poly2 (struct laurentpoly2 **lpt, struct laurentpoly2 **extradets, int extranum)
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
    print_laurentpoly2 (*lpt, 'u', 'v');
    printf ("\n");
    for (j = 0; j < extranum; j++)
    {
      printf ("Fundamental ideal factor %d:\n", j);
      print_laurentpoly2 (extradets[j], 'u', 'v');
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
  struct laurentpoly2 *l2;
  struct laurentpoly *l1;
  int sign, tok, ffound;

  tok = gettoken (file);

  assert (tok == TOK_ALEXANDER || tok == TOK_IDEAL);

  tok = gettoken (file);
  assert (tok == TOK_LPAREN);
  ai = (struct alexanderideal *) malloc (sizeof (struct alexanderideal));
  ai->l1num = ai->l2num = ai->fl2num = 0;
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
        l1 = read_laurentpoly (file, indet_names);
        ai->l1[ai->l1num++] = l1;
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
          ai->fl2[ai->fl2num++] = l2;
        } else {
          ai->l2[ai->l2num++] = l2;
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

