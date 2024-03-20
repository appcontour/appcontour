#include <assert.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#ifdef HAVE_UNISTD_H
  #include <unistd.h>
  #include <sys/types.h>
  #include <sys/wait.h>
#endif

#include "contour.h"
#include "parser.h"
#include "mendes.h"
#include "fundamental.h"
#include "fox.h"
#include "alexander.h"
#include "representations.h"
#include "giovecanonify.h"
#include "readdtcode.h"
#include "wirtinger.h"

#ifndef EXAMPLES_DIR
  #define EXAMPLES_DIR ""
#endif
#define MAXFILELENGTH 2000
#define NCCMAX 100

#define KITANO_SUZUKI_PREFIX "ks_"

int debug = 0;
int quiet = 0;
int outformat = 0;
int verbose = 0;
int interactive = 0;
int heisemberg = -1;
int docanonify = 1;
int useoldcanonify = 0;
static int renumber = 1;

struct global_data globals;
struct tagged_data user_data;

/* prototipi locali */
FILE *open_description_file (char *arg, int argnum);
FILE *new_file (FILE *oldfile, FILE *newfile);
static char examplesfilename[MAXFILELENGTH+1];
int readintlistsub1 (int nccmax, int *ccids, char *s);

int
main (int argc, char *argv[])
{
  struct sketch *sketch, *s1, *s2;
  int tok, ori = 0, i, count, action=ACTION_NONE, res, ncc = 1, ccids[NCCMAX];
  unsigned int rndseed = 0;
  int newextregion = 0;
  int fg_type = FG_SURFACE;
  char rule[10], *endch, *envvar;
  FILE *infile = 0;
  FILE *infile2 = 0;
  struct mendesgraph *mendes;
  struct region *r;
  struct borderlist *bl;
  struct border *bp;
  struct arc *a, *a2;
  struct ccomplex *ccomplex;
  int infiles = 1;
  struct presentation *p;
  struct alexanderideal *ai;
  struct vecofintlist *loiv, *newloiv;
  int ccnum;
  struct ccomplexcc *cccc;

  ccids[0] = 0;
  user_data.mrnum = user_data.manum = 0;
  globals.dtcode_realize = 0;
  globals.dorecomputef = globals.doretagregions = globals.simplifypresentation = 1;
  globals.simplifycomplex = globals.simplifyideal = 1;
  globals.finfinity = globals.preabelian = globals.nobasecanonify = globals.viacc = 0;
  globals.mendesge = HGE_TEXT;
  globals.foxd = FOXD_UNDEF;
  globals.rulenames = RULENAMES_NEW;
  globals.shuffle = globals.autosurgery = 0;
  globals.onlyeven = globals.dontidentify = globals.insist = globals.outer = 0;
  globals.focus_on_fundamental = globals.principal = globals.factorideal = globals.internalcheck = 0;
  globals.abelianize = globals.experimental = 0;
  globals.knotname_fallback = 1;
  globals.summand1cc = globals.summand2cc = -1;
  if ((envvar = getenv ("APPCONTOUR_AUTOSURGERY")) && *envvar) 
    globals.autosurgery++;
  if ((envvar = getenv ("APPCONTOUR_OLDNAMES")) && *envvar) 
    globals.rulenames = RULENAMES_OLD;
  if ((envvar = getenv ("APPCONTOUR_DISABLEKNFALLBACK")) && *envvar) 
    globals.knotname_fallback = 0;
  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"--newnames") == 0)
    {
      globals.rulenames = RULENAMES_NEW;
      continue;
    }
    if (strcmp(argv[i],"--oldnames") == 0)
    {
      globals.rulenames = RULENAMES_OLD;
      continue;
    }
    if (strcmp(argv[i],"--summand1cc") == 0)
    {
      globals.summand1cc = strtol (argv[++i], &endch, 10) - 1;
      continue;
    }
    if (strcmp(argv[i],"--summand2cc") == 0)
    {
      globals.summand2cc = strtol (argv[++i], &endch, 10) - 1;
      continue;
    }
    if (strcmp(argv[i],"-r") == 0 || strcmp(argv[i],"--region") == 0)
    {
      user_data.region[user_data.mrnum] = strtol (argv[++i], &endch, 10);
      if (*endch == ':' && *(++endch)) 
	user_data.stratum[user_data.stnum++] = strtol (endch, &endch, 10);
      if (*endch != '\0') {
        fprintf (stderr, "Conversion error in %s, %x\n", argv[i], (int)*endch);
	exit (1);
      }
      user_data.mrnum++;
      if (user_data.mrnum > MRPTMAX || user_data.stnum > MSPTMAX)
      {
        fprintf (stderr, "Too many tagged regions (%d)\n", user_data.mrnum);
        exit (1);
      }
      continue;
    }
    if (strcmp(argv[i],"-a") == 0 || strcmp(argv[i],"--arc") == 0)
    {
      user_data.arc[user_data.manum] = strtol (argv[++i], &endch, 10);
      user_data.arcl[user_data.manum] = 0;
      if (*endch == ':' && *(++endch)) 
	user_data.arcl[user_data.manum] = strtol (endch, &endch, 10);
      if (*endch != '\0') {
        fprintf (stderr, "Conversion error in %s, %x\n", argv[i], (int)*endch);
	exit (1);
      }
      user_data.manum++;
      if (user_data.manum > MAPTMAX)
      {
        fprintf (stderr, "Too many tagged arcs (%d)\n", user_data.manum);
        exit (1);
      }
      continue;
    }
    if (strcmp(argv[i],"--stratum") == 0 || strcmp(argv[i],"-s") == 0)
    {
      user_data.stratum[user_data.stnum++] = atoi (argv[++i]);
      if (user_data.stnum > MSPTMAX)
      {
        fprintf (stderr, "Too many tagged strata (%d)\n", user_data.stnum);
        exit (1);
      }
      continue;
    }
    if (strcmp(argv[i],"--version") == 0)
    {
      printf ("%s version %s\n", argv[0], PACKAGE_VERSION);
      return (0);
    }
    if (strcmp(argv[i],"-Q") == 0)
    {
      quiet = 1;
      outformat = OUTFORMAT_APPCONTOUR;
      continue;
    }
    if (strcmp(argv[i],"--M2") == 0)
    {
      quiet = 1;
      outformat = OUTFORMAT_MACAULAY2;
      continue;
    }
    if (strcmp(argv[i],"-q") == 0)
    {
      quiet = 1;
      continue;
    }
    if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0)
    {
      verbose++;
      continue;
    }
    if (strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--interactive") == 0)
    {
      interactive++;
      continue;
    }
    if (strcmp(argv[i],"--right") == 0)
    {
      globals.dtcode_realize = 1;
      continue;
    }
    if (strcmp(argv[i],"--left") == 0)
    {
      globals.dtcode_realize = -1;
      continue;
    }
    if (strcmp(argv[i],"--autosurgery") == 0)
    {
      globals.autosurgery++;
      continue;
    }
    if (strcmp(argv[i],"--noautosurgery") == 0)
    {
      globals.autosurgery = 0;
      continue;
    }
    if (strcmp(argv[i],"--principal") == 0 || strcmp(argv[i],"--gcd") == 0)
    {
      globals.principal++;
      continue;
    }
    if (strcmp(argv[i],"--factorideal") == 0 || strcmp(argv[i],"--factor") == 0)
    {
      globals.factorideal++;
      continue;
    }
    if (strcmp(argv[i],"--internalcheck") == 0)
    {
      globals.internalcheck++;
      continue;
    }
    if (strcmp(argv[i],"--abelianize") == 0)
    {
      globals.abelianize = 1;
      continue;
    }
    if (strcmp(argv[i],"--trivialize") == 0)
    {
      globals.abelianize = 2;
      continue;
    }
    if (strcmp(argv[i],"--maxd") == 0)
    {
      globals.foxd = FOXD_MAXINTERESTING;
      continue;
    }
    if (strcmp(argv[i],"--dontidentify") == 0)
    {
      globals.dontidentify++;
      continue;
    }
    if (strcmp(argv[i],"--outer") == 0)
    {
      globals.outer++;
      continue;
    }
    if (strcmp(argv[i],"--inner") == 0)
    {
      fprintf (stderr, "Warning: option --inner is now implied by default.  Use --outer to obtain the previous default behaviour\n");
      globals.outer = 0;
      continue;
    }
    if (strcmp(argv[i],"--insist") == 0)
    {
      globals.insist++;
      continue;
    }
    if (strcmp(argv[i],"--experimental") == 0)
    {
      globals.experimental++;
      continue;
    }
    if (strcmp(argv[i],"--debug") == 0)
    {
      debug++;
      continue;
    }
    if (strcmp(argv[i],"--nocanonify") == 0)
    {
      docanonify = 0;
      continue;
    }
    if (strcmp(argv[i],"--nosimplify") == 0)
    {
      globals.simplifypresentation = globals.simplifycomplex
                                   = globals.simplifyideal = 0;
      continue;
    }
    if (strcmp(argv[i],"--nopresentationsimplify") == 0)
    {
      globals.simplifypresentation = 0;
      continue;
    }
    if (strcmp(argv[i],"--nocomplexsimplify") == 0)
    {
      globals.simplifycomplex = 0;
      continue;
    }
    if (strcmp(argv[i],"--noidealsimplify") == 0)
    {
      globals.simplifyideal = 0;
      continue;
    }
    if (strcmp(argv[i],"--nobasecanonify") == 0)
    {
      globals.nobasecanonify = 1;
      continue;
    }
    if (strcmp(argv[i],"--shuffle") == 0)
    {
      globals.shuffle = 1;
      continue;
    }
    if (strcmp(argv[i],"--inside") == 0 || strcmp(argv[i],"--in") == 0)
    {
      fg_type = FG_INTERNAL;
      continue;
    }
    if (strcmp(argv[i],"--outside") == 0 || strcmp(argv[i],"--out") == 0)
    {
      fg_type = FG_EXTERNAL;
      continue;
    }
    if (strcmp(argv[i],"--surface") == 0)
    {
      fg_type = FG_SURFACE;
      continue;
    }
    if (strcmp(argv[i],"--oldcanonify") == 0)
    {
      useoldcanonify = 1;
      continue;
    }
    if (strcmp(argv[i],"--dontrenumber") == 0)
    {
      renumber = 0;
      continue;
    }
    if (strcmp(argv[i],"--foxd") == 0)
    {
      globals.foxd = atoi (argv[++i]);
      continue;
    }
    if (strcmp(argv[i],"--seed") == 0)
    {
      rndseed = atoi (argv[++i]);
      continue;
    }
    if (strcmp(argv[i],"--heisemberg") == 0 ||
        strcmp(argv[i],"--ti") == 0 ||
        strcmp(argv[i],"--transfer_islands") == 0)
    {
      heisemberg = atoi (argv[++i]);
      continue;
    }
    if (strcmp(argv[i],"--darkmatter") == 0 ||
        strcmp(argv[i],"--finfinity") == 0 ||
        strcmp(argv[i],"--fi") == 0)
    {
      globals.finfinity = atoi (argv[++i]);
      continue;
    }
    if (strcmp(argv[i],"--preabelian") == 0)
    {
      globals.preabelian++;
      continue;
    }
    if (strcmp(argv[i],"--mendes_ge") == 0 ||
        strcmp(argv[i],"--mge") == 0)
    {
      i++;
      if (strcmp(argv[i],"text") == 0) globals.mendesge = HGE_TEXT;
      else if (strcmp(argv[i],"pykig") == 0) globals.mendesge = HGE_PYKIG;
      else if (strcmp(argv[i],"kig") == 0) globals.mendesge = HGE_KIG;
      else
      {
        printf ("Invalid mendes graphic engine selection: %s\n",
                    argv[i]);
        printf ("Valid choices are: text(default), pykig\n");
        exit (111);
      }
      continue;
    }
    if (strcmp(argv[i],"--help") == 0)
    {
      printf ("Usage: %s [options] action [file] [file2]\n", argv[0]);
      printf (" Possible actions:\n");
      printf ("  rule <rule>: apply indicated rule to contour\n");
      printf ("  mergearcs: apply an inverse rule that merges two arcs\n");
      printf ("  wrinkle: apply inverse L (lip) rule\n");
      printf ("  punchhole/removehole: perform vertical surgery\n");
      printf ("  gluearcs/pinchneck: perform horizontal surgery\n");
      printf ("  addsphere/removesphere: add-remove small sphere\n");
      printf ("  wrap: put the 3D surface into a big sphere\n");
      printf ("  extractcc <int>: extract 3D connected component\n");
      printf ("  removecc <int>: remove 3D connected component from contour\n");
      printf ("  leftright: left-right reflection\n");
      printf ("  frontback: front-back reflection\n");
      printf ("  mendes: compute Mendes graph (see Hacon-Mendes-Romero Fuster)\n");
      printf ("  evert <int>: make region <int> become the unbounded region\n");
      printf ("  union: disjoint union of two apparent contours\n");
      printf ("  sum: connected sum of two apparent contours\n");
      printf ("  3devert <int>: evert 3D surface (an internal component becomes the outside)\n");
      printf ("    <int> indicates a region to be used for the 3d eversion,\n");
      printf ("    use \"-1\" to indicate the first region with positive f\n");
      printf ("  knotsum: connected sum of two knotted tori (as connected sum of knots)\n");
      printf ("\n Possible informational actions:\n");
      printf ("  info, characteristic, rules, iscontour, islabelled, countcc\n");
      printf ("  list[ma|invl|invs|strata]\n");
      printf ("  ccorientation <int>: gives the orientation of a 3D connected component\n");
      printf ("  ccparent <cc>: finds the 3D component that directly contains \"cc\"\n");
      printf ("    0 means that \"cc\" is external\n");
      printf ("  ccordering: show 3D inclusion between the connected components\n");
      printf ("  compare: lexicographic comparison between two contours\n");
      printf ("\n Possible conversion and standardization actions:\n");
      printf ("  print, printmorse, knot2morse, any2morse, canonify, giovecanonify\n");
      printf ("\n Cell complex and fundamental group:\n");
      printf ("  cellcomplex, insidecomplex, outsidecomplex\n");
      printf ("  fundamental, insidefundamental, outsidefundamental\n");
      printf ("   abbreviations: fg, ifg, ofg\n");
      printf ("  alexander: compute Alexander polynomial of fundamental group\n");
      printf ("  linkingnumber: compute linking number from fundamental group\n");
      printf ("  deficiency: print deficiency of computed fundamental group\n");
      printf ("  abelianizedfundamental, insideabelianizedfundamental, outsideabelianizedfundamental\n");
      printf ("   abbreviations: afg, iafg, oafg\n");
      printf ("  scharacteristic, icharacteristic, ocharacteristic\n");
      printf ("   abbreviations: sch, ich, och\n");
      printf ("  suggest_p_surgery: display 'punchhole' surgery that does not affect fundamental group of inside\n");
      printf ("  specific options: --in, --out, --surface[default]\n");
      printf ("   indicate which part of space to consider.  E.g. \"ifg\" is equivalent to \"fg --inside\"\n");
      printf ("   --[no]autosurgery: automatically apply punchhole surgeries to increase the initial presentation deficiency\n");
      printf ("\n File2 can be present only for actions that require two descriptions:\n");
      printf ("  'compare', 'union', 'sum' actions.\n");
      printf ("  alternatively for such actions the two descriptions can be contained\n");
      printf ("  in the same file, typically standard input.\n");
      printf ("\n Invariants based on representations:\n");
      printf ("  ks_<group>: compute Kitano-Suzuki invariant with respect to some group,\n");
      printf ("  [See source code for a list of implemented groups and corresponding options]\n");
      printf ("\n Reading from knot tables by knot name: %s [options] action knotname\n", argv[0]);
      printf ("  Examples of knotnames include:\n");
      printf ("    K<x>_<n> (Rofsen table),\n");
      printf ("    K<x>a_<n>, K<x>n_<n> (Dowker-Thistlethwaite numbering)\n");
      printf ("    L<x>a_<n>, L<x>n_<n> (Thistlethwaite link table)\n");
      printf ("    HK<x>_<n> (Handlebody genus 2 knots up to 6 crossings of Ishii-Kishimoto-Moriuchi)\n");
      printf ("  Possible specific actions include:\n");
      printf ("    dtcode\n");
      printf ("    rdtcode\n");
      printf ("    gausscode\n");
      printf ("\n Possible options are:\n");
      printf ("  --help: this help\n");
      printf ("  --version: print program version\n");
      printf ("  -q: be quiet\n");
      printf ("  -v|--verbose: be more verbose\n");
      printf ("  -Q: use appcontour input syntax for fundamental groups and Alexander polynomials\n");
      printf ("  --M2: use Macaulay2 input syntax for Alexander polynomials/ideals\n");
      printf ("  --nocanonify: do not canonify region description before printing\n");
      printf ("  --oldcanonify: use the old (version <= 1.3.0) canonification procedure\n");
      printf ("  --dontrenumber: do not renumber regions and arcs after giovecanonify\n");
      printf ("  --transfer_islands|--ti <int_coded_flags>: information on island\n");
      printf ("      location in case of ambiguity (e.g. rule C2)\n");
      printf ("  --finfinity|--fi <int>: value of f at infinity (default 0)\n");
      printf ("  --seed <int>: initialize random number generator\n");
      printf ("      e.g. for Mendes graph graphic presentation or for random base change\n");
      printf ("  -r|--region <int>: mark region for specific action\n");
      printf ("  -a|--arc <int>: mark arc for specific action\n");
      printf ("  --oldnames|--newnames: select set of names for rules\n");
      printf ("  --preabelian: compute preabelian presentation of fundamental group\n");
      printf ("  --in|--out: apply command to the inside/outside of surface\n");
      printf ("      works for cell complex and fundamental group computations\n");
      printf ("  --foxd <int>: compute the d-th elementary ideal (d given by <int>) instead of the\n");
      printf ("      default value (d=1 for knots and links, d=2 for genus-2 surfaces, etc)\n");
      printf ("  --maxd: compute the elementary ideal corresponding to 1x1 minors of the presentation\n");
      printf ("      matrix\n");
      printf ("  --nopresentationsimplify: do not simplify the presentation of the fundamental group\n");
      printf ("  --nocomplexsimplify: do not try to collapse the cell complex\n");
      printf ("  --noidealsimplify: do not simplify the alexander ideal\n");
      printf ("  --nosimplify: all three simplifications above\n");
      printf ("  --nobasecanonify: do not base-canonify Alexander polynomial in two indeterminates\n");
      printf ("  --shuffle: random change of base for Alexander ideal in two indeterminates\n");
      printf ("\n If 'file' is not given, description is taken from standard input\n");
      exit (0);
    }
    if (*argv[i] == '-')
    {
      printf ("Invalid option %s\n", argv[i]);
      exit (1);
    }
    if (infile2)
    {
      printf ("Too many arguments: %s\n", argv[i]);
      exit (1);
    }
    if (action != ACTION_NONE)
    {
      if (infile == 0)
      {
        infile = open_description_file (argv[i], 1);
      } else infile2 = open_description_file (argv[i], 2);
      if (action == ACTION_FILEPATH)
      {
         printf ("%s\n", examplesfilename);
      }
    }
    if (strcmp(argv[i],"applyrule") == 0 || strcmp(argv[i],"rule") == 0)
    {
      action = ACTION_APPLYRULE;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a rule\n"); exit (11);}
      fprintf (stderr, "applying rule %s\n", argv[i]);
      strncpy (rule, argv[i], 9);
    }
    if (strcmp(argv[i],"extractcc") == 0)
    {
      action = ACTION_EXTRACTCC;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      ncc = readintlistsub1 (NCCMAX, ccids, argv[i]);
      //*ccids = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"removecc") == 0)
    {
      action = ACTION_REMOVECC;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      *ccids = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"ccorientation") == 0)
    {
      action = ACTION_CCORIENTATION;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      *ccids = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"ccparent") == 0)
    {
      action = ACTION_FINDCCPARENT;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      *ccids = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"ccchilds") == 0)
    {
      action = ACTION_CCCHILDS;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a cc id\n"); exit (11);}
      *ccids = atoi (argv[i]) - 1;
    }
    if (strcmp(argv[i],"listholes") == 0) action = ACTION_LISTHOLES;
    if (strcmp(argv[i],"removehole") == 0) action = ACTION_REMOVEHOLE;
    if (strcmp(argv[i],"listspheres") == 0) action = ACTION_LISTSPHERES;
    if (strcmp(argv[i],"removesphere") == 0) action = ACTION_REMOVESPHERE;
    if (strcmp(argv[i],"liststrata") == 0) action = ACTION_LISTSTRATA;
    if (strcmp(argv[i],"addsphere") == 0) action = ACTION_ADDSPHERE;
    if (strcmp(argv[i],"wrap") == 0) action = ACTION_WRAP;
    if (strcmp(argv[i],"union") == 0) {action = ACTION_UNION; infiles = 2;}
    if (strcmp(argv[i],"connectedsum") == 0) {action = ACTION_CONNECTEDSUM; infiles = 2;}
    if (strcmp(argv[i],"sum") == 0) {action = ACTION_CONNECTEDSUM; infiles = 2;}
    if (strcmp(argv[i],"knotsum") == 0) {action = ACTION_KNOTSUM; infiles = 2;}
    if (strcmp(argv[i],"punchhole") == 0) action = ACTION_PUNCHHOLE;
    if (strcmp(argv[i],"mergearcs") == 0) action = ACTION_MERGEARCS;
    if (strcmp(argv[i],"gluearcs") == 0) action = ACTION_GLUEARCS;
    if (strcmp(argv[i],"pinchneck") == 0) action = ACTION_PINCHNECK;
    if (strcmp(argv[i],"listma") == 0) action = ACTION_LISTMA;
    if (strcmp(argv[i],"listmergearcs") == 0) action = ACTION_LISTMA;
    if (strcmp(argv[i],"wrinkle") == 0) action = ACTION_WRINKLE;
    if (strcmp(argv[i],"listwr") == 0) action = ACTION_LISTWR;
    if (strcmp(argv[i],"listinvc1") == 0) action = ACTION_LISTWR;
    if (strcmp(argv[i],"listinvl") == 0) action = ACTION_LISTWR;
    if (strcmp(argv[i],"swallowtail") == 0) action = ACTION_SWALLOWTAIL;
    if (strcmp(argv[i],"listinvcn1") == 0) action = ACTION_LISTST;
    if (strcmp(argv[i],"listinvs") == 0) action = ACTION_LISTST;
    if (strcmp(argv[i],"listst") == 0) action = ACTION_LISTST;
    if (strcmp(argv[i],"listswallowtails") == 0) action = ACTION_LISTST;
    if (strcmp(argv[i],"invcn3") == 0) action = ACTION_PUNCTURE;
    if (strcmp(argv[i],"invc") == 0) action = ACTION_PUNCTURE;
    if (strcmp(argv[i],"puncture") == 0) action = ACTION_PUNCTURE;
    if (strcmp(argv[i],"listinvcn3") == 0) action = ACTION_LISTPU;
    if (strcmp(argv[i],"listinvc") == 0) action = ACTION_LISTPU;
    if (strcmp(argv[i],"listpunctures") == 0) action = ACTION_LISTPU;
    if (strcmp(argv[i],"testallrules") == 0) action = ACTION_TESTALLRULES;
    if (strcmp(argv[i],"rules") == 0) action = ACTION_TESTALLRULES;
    if (strcmp(argv[i],"countcc") == 0) action = ACTION_COUNTCC;
    if (strcmp(argv[i],"ccordering") == 0) action = ACTION_CCORDERING;
    if (strcmp(argv[i],"print") == 0) action = ACTION_PRINTSKETCH;
    if (strcmp(argv[i],"ishuffman") == 0) action = ACTION_ISHUFFMAN;
    if (strcmp(argv[i],"islabelled") == 0) action = ACTION_ISHUFFMAN;
    if (strcmp(argv[i],"isappcon") == 0) action = ACTION_ISHUFFMAN;
    if (strcmp(argv[i],"iscontour") == 0) action = ACTION_ISCONTOUR;
    if (strcmp(argv[i],"compare") == 0) {action = ACTION_COMPARE; infiles = 2;}
    if (strcmp(argv[i],"canonify") == 0) action = ACTION_CANONIFY;
    if (strcmp(argv[i],"giovecanonify") == 0) action = ACTION_GIOVECANONIFY;
    if (strcmp(argv[i],"knot2morse") == 0) action = ACTION_KNOT2MORSE;
    if (strcmp(argv[i],"code") == 0) action = ACTION_CODE;
    if (strcmp(argv[i],"dtcode") == 0) action = ACTION_DTCODE;
    if (strcmp(argv[i],"rdtcode") == 0) action = ACTION_RDTCODE;
    if (strcmp(argv[i],"gausscode") == 0) action = ACTION_GAUSSCODE;
    if (strcmp(argv[i],"wirtinger") == 0) action = ACTION_WIRTINGER;
    if (strcmp(argv[i],"any2morse") == 0) action = ACTION_ANY2MORSE;
    if (strcmp(argv[i],"printmorse") == 0) action = ACTION_PRINTMORSE;
    if (strcmp(argv[i],"characteristic") == 0) action = ACTION_CHARACTERISTIC;
    if (strcmp(argv[i],"info") == 0) action = ACTION_INFO;
    if (strcmp(argv[i],"frontback") == 0) action = ACTION_FRONTBACK;
    if (strcmp(argv[i],"leftright") == 0) action = ACTION_LEFTRIGHT;
    if (strcmp(argv[i],"mendes") == 0) action = ACTION_MENDES;
    if (strcmp(argv[i],"cellcomplex") == 0) action = ACTION_CELLCOMPLEX;
    if (strcmp(argv[i],"insidecomplex") == 0) {action = ACTION_CELLCOMPLEX; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"outsidecomplex") == 0) {action = ACTION_CELLCOMPLEX; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"fg") == 0) action = ACTION_FUNDAMENTAL;
    if (strcmp(argv[i],"fundamental") == 0) action = ACTION_FUNDAMENTAL;
    if (strcmp(argv[i],"sfundamental") == 0) action = ACTION_FUNDAMENTAL;
    if (strcmp(argv[i],"ifg") == 0) {action = ACTION_FUNDAMENTAL; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"insidefundamental") == 0) {action = ACTION_FUNDAMENTAL; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"ofg") == 0) {action = ACTION_FUNDAMENTAL; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"outsidefundamental") == 0) {action = ACTION_FUNDAMENTAL; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"deficiency") == 0) action = ACTION_DEFICIENCY;
    if (strcmp(argv[i],"afg") == 0) action = ACTION_AFUNDAMENTAL;
    if (strcmp(argv[i],"abelianizedfundamental") == 0) action = ACTION_AFUNDAMENTAL;
    if (strcmp(argv[i],"iafg") == 0) {action = ACTION_AFUNDAMENTAL; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"insideabelianizedfundamental") == 0) {action = ACTION_AFUNDAMENTAL; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"oafg") == 0) {action = ACTION_AFUNDAMENTAL; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"outsideabelianizedfundamental") == 0) {action = ACTION_AFUNDAMENTAL; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"foxjacobian") == 0) action = ACTION_FOXJACOBIAN;
    if (strcmp(argv[i],"alexander") == 0) action = ACTION_ALEXANDER;
    if (strcmp(argv[i],"linkingnumber") == 0) action = ACTION_LINKINGNUMBER;
    if (strcmp(argv[i],"scharacteristic") == 0) {action = ACTION_CHARACTERISTIC; globals.viacc = 1;}
    if (strcmp(argv[i],"sch") == 0) action = ACTION_CHARACTERISTIC;
    if (strcmp(argv[i],"icharacteristic") == 0) {action = ACTION_CHARACTERISTIC; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"ich") == 0) {action = ACTION_CHARACTERISTIC; fg_type=FG_INTERNAL;}
    if (strcmp(argv[i],"ocharacteristic") == 0) {action = ACTION_CHARACTERISTIC; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"och") == 0) {action = ACTION_CHARACTERISTIC; fg_type=FG_EXTERNAL;}
    if (strcmp(argv[i],"suggest_p_surgery") == 0) action = ACTION_SUGGEST_P_SURGERY;
    if (strcmp(argv[i],"filepath") == 0) action = ACTION_FILEPATH;
    if (strcmp(argv[i],"newfeature") == 0) action = ACTION_NEWFEATURE;
    if (strcmp(argv[i],"evert") == 0)
    {
      action = ACTION_EVERT;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a region tag\n"); exit (11);}
      newextregion = atoi (argv[i]);
    }
    if (strcmp(argv[i],"3devert") == 0)
    {
      action = ACTION_3DEVERT;
      i++;
      if (i >= argc) {fprintf (stderr, "specify a region tag\n"); exit (11);}
      newextregion = atoi (argv[i]);
    }
    if (strncasecmp (argv[i], KITANO_SUZUKI_PREFIX, strlen(KITANO_SUZUKI_PREFIX)) == 0) /* one ring to rule them all :-) */
    {
      char *group = argv[i] + strlen(KITANO_SUZUKI_PREFIX);
      if (strncasecmp (group, "sl2,", 4) == 0)
      {
        group += 4;
        globals.p = atoi (group);
        action = ACTION_CCCOUNTSL2ZP;
        if (globals.p < 2) {fprintf (stderr, "value %d of p in SL2,p is invalid\n", globals.p); exit (11);}
      }
      if (strncasecmp (group, "s", 1) == 0 && isdigit (group[1]))
      {
        group += 1;
        globals.n = atoi (group);
        action = ACTION_CCCOUNTSN;
        if (globals.n < 2) {fprintf (stderr, "value %d of n in Sn is invalid\n", globals.n); exit (11);}
      }
      if (strncasecmp (group, "a", 1) == 0 && isdigit (group[1]))
      {
        group += 1;
        globals.n = atoi (group);
        globals.onlyeven = 1;
        action = ACTION_CCCOUNTSN;
        if (globals.n < 2) {fprintf (stderr, "value %d of n in An is invalid\n", globals.n); exit (11);}
      }
      if (strncasecmp (group, "psl2,", 5) == 0)
      {
        group += 5;
        globals.q = atoi (group);
        action = ACTION_CCCOUNTPSL2Q;;
        if (globals.q < 2) {fprintf (stderr, "value %d of q in PSL2,q is invalid\n", globals.q); exit (11);}
      }
    }
    if (strcasecmp(argv[i],"countsn") == 0)
    {
      action = ACTION_CCCOUNTSN;
      i++;
      if (i >= argc) {fprintf (stderr, "specify the value of n\n"); exit (11);}
      globals.n = atoi (argv[i]);
      fprintf (stderr, "Command 'countsn %d' is deprecated, use 'ks_S%d' instead\n", globals.n, globals.n);
    }
    if (strcasecmp(argv[i],"countan") == 0)
    {
      globals.onlyeven = 1;
      action = ACTION_CCCOUNTSN;
      i++;
      if (i >= argc) {fprintf (stderr, "specify the value of n\n"); exit (11);}
      globals.n = atoi (argv[i]);
      fprintf (stderr, "Command 'countan %d' is deprecated, use 'ks_A%d' instead\n", globals.n, globals.n);
    }
    if (strcasecmp(argv[i],"countsl2q") == 0 || strcasecmp(argv[i],"countsl2zp") == 0)
    {
      action = ACTION_CCCOUNTSL2ZP;
      i++;
      if (i >= argc) {fprintf (stderr, "specify the value of the prime p\n"); exit (11);}
      globals.p = atoi (argv[i]);
      fprintf (stderr, "Command 'countsl2zp %d' is deprecated, use 'ks_SL2,%d' instead\n", globals.p, globals.p);
    }
    if (strcasecmp(argv[i],"countpsl2q") == 0 || strcasecmp(argv[i],"countpsl2zp") == 0)
    {
      action = ACTION_CCCOUNTPSL2Q;
      i++;
      if (i >= argc) {fprintf (stderr, "specify the value of the prime p\n"); exit (11);}
      globals.q = atoi (argv[i]);
      fprintf (stderr, "Command 'countpsl2q %d' is deprecated, use 'ks_PSL2,%d' instead\n", globals.q, globals.q);
    }
    if (action == ACTION_NONE)
    {
      fprintf (stderr, "invalid arg[%d] = %s\n", i, argv[i]);
      exit (10);
    }
  }
  if (action == ACTION_NONE) action = ACTION_PRINTSKETCH;

  if (infile == 0) infile = stdin;
  if (action == ACTION_FILEPATH)
  {
    if (infile == stdin) printf ("<stdin>\n");
    exit (0);
  }

  if (infiles < 2 && infile2 != 0)
  {
    fprintf (stderr, "Too many arguments\n");
    exit (1);
  }

  srandom (rndseed);
  switch (action)
  {
    case ACTION_PRINTSKETCH:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    // printsketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_APPLYRULE:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    res = apply_rule (rule, sketch);
    printsketch (sketch);
    if (res == 0)
    {
      fprintf (stderr, "Rule does not match!\n");
      exit(14);
    }
    break;

    case ACTION_LISTSPHERES:
    ori = 1;
    case ACTION_LISTHOLES:
    if (ori != 1) ori = -1;
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    for (r = sketch->regions; r; r = r->next) {
      bl = r->border;
      if (bl->next) continue;
      bp = bl->sponda;
      if (bp == 0) continue;
      if (bp->next != bp) continue;
      if (bp->orientation * ori <= 0) continue;
      if (bp->info->cusps != 0) continue;
      if (quiet) printf ("%d\n", r->tag);
         else printf ("Region %d is a %s\n", r->tag,
                      (ori > 0)?"sphere":"hole");
    }
    break;

    case ACTION_REMOVESPHERE:
    ori = 1;
    case ACTION_REMOVEHOLE:
    if (ori != 1) ori = -1;
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.mrnum < 1) {
      fprintf (stderr, "You must specify a region\n");
      fprintf (stderr, "   using options -r\n");
      exit(15);
    }
    r = findregion (sketch, user_data.region[0]);
    if (r == 0) fprintf (stderr, "Cannot find region %d\n",
                        user_data.region[0]);
    if (! (r)) exit(15);
    bl = r->border;
    if (bl->next) exit (16);
    bp = bl->sponda;
    if (bp == 0) exit (17);
    if (bp->next != bp) exit (18);
    if (bp->orientation * ori <= 0) exit (19);
    if (bp->info->cusps != 0) exit (20);
    remove_s1 (bp, sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_LISTSTRATA:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    for (r = sketch->regions; r; r = r->next) {
      if (quiet) printf ("%d:%d\n", r->tag, r->f);
         else printf ("Region %d has %d strata\n", r->tag, r->f);
    }
    break;

    case ACTION_LISTPU:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (sketch->huffman_labelling) list_punctures (sketch);
    break;

    case ACTION_PUNCTURE:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.manum < 2) {
      fprintf (stderr, "You must specify two arcs\n");
      fprintf (stderr, "   using options -a\n");
      list_punctures (sketch);
      exit(15);
    }
    a = findarc (sketch, user_data.arc[0]);
    a2 = findarc (sketch, user_data.arc[1]);
    if (a == 0) fprintf (stderr, "Cannot find arc %d\n", user_data.arc[0]);
    if (a2 == 0) fprintf (stderr, "Cannot find arc %d\n", user_data.arc[1]);
    if (! (a && a2)) exit(15);
    res = apply_puncture (sketch, a, a2,
        user_data.arcl[0], user_data.arcl[1], 0);
    printsketch (sketch);
    if (res == 0)
    {
      fprintf (stderr, "Cannot create puncture!\n");
      exit(14);
    }
    break;

    case ACTION_LISTST:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (sketch->huffman_labelling) list_swallowtails (sketch);
    break;

    case ACTION_SWALLOWTAIL:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.manum < 1) {
      fprintf (stderr, "You must specify an arc\n");
      fprintf (stderr, "   using options -a n:i\n");
      list_swallowtails (sketch);
      exit(15);
    }
    a = findarc (sketch, user_data.arc[0]);
    if (a == 0) fprintf (stderr, "Cannot find arc %d\n", user_data.arc[0]);
    if (! (a)) exit(15);
    res = apply_createswallowtail (sketch, a,  
	user_data.arcl[0], 1, 0);
    printsketch (sketch);
    if (res == 0)
    {
      fprintf (stderr, "Cannot create swallowtail!\n");
      exit(14);
    }
    break;

    case ACTION_LISTWR:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (sketch->huffman_labelling) list_strata (sketch);
    break;

    case ACTION_WRINKLE:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.mrnum < 1 || user_data.stnum < 1) {
      fprintf (stderr, "You must specify a region and a stratum\n");
      fprintf (stderr, "   using options -r and --stratum\n");
      list_strata (sketch);
      exit(15);
    }
    r = findregion (sketch, user_data.region[0]);
    if (r == 0) fprintf (stderr, "Cannot find region %d\n", 
			user_data.region[0]);
    if (! (r)) exit(15);
    res = apply_createwrinkle (sketch, r,  
	user_data.stratum[0], 0);
    printsketch (sketch);
    if (res == 0)
    {
      fprintf (stderr, "Cannot create wrinkle!\n");
      exit(14);
    }
    break;

    case ACTION_LISTMA:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (sketch->huffman_labelling) list_mergearcs (sketch);
    break;

    case ACTION_MERGEARCS:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.mrnum < 1 || user_data.manum < 2) {
      fprintf (stderr, "You must specify a region and two arcs\n");
      fprintf (stderr, "   using options -r and -a\n");
      list_mergearcs (sketch);
      exit(15);
    }
    r = findregion (sketch, user_data.region[0]);
    a = findarc (sketch, user_data.arc[0]);
    a2 = findarc (sketch, user_data.arc[1]);
    if (r == 0) fprintf (stderr, "Cannot find region %d\n", 
			user_data.region[0]);
    if (a == 0) fprintf (stderr, "Cannot find arc %d\n", user_data.arc[0]);
    if (a2 == 0) fprintf (stderr, "Cannot find arc %d\n", user_data.arc[1]);
    if (! (r && a && a2)) exit(15);
    res = apply_mergearcs (sketch, r, a, a2, 
	user_data.arcl[0], user_data.arcl[1], 0);
    printsketch (sketch);
    if (res == 0)
    {
      fprintf (stderr, "Cannot merge arcs!\n");
      exit(14);
    }
    break;

    case ACTION_TESTALLRULES:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (testallrules (sketch) == 0) exit(15);
    break;

    case ACTION_PINCHNECK:
    ori = 1;
    case ACTION_GLUEARCS:
    if (ori != 1) ori = -1;
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.manum < 2) {
      if (! quiet) printf ("You must specify two arcs using option -a:\n");
      res = list_hor_sur (sketch, ori);
      break;
    }
    a = findarc (sketch, user_data.arc[0]);
    a2 = findarc (sketch, user_data.arc[1]);
    if (a == 0 || a2 == 0) {
      fprintf (stderr, "Cannot find arc %d\n", 
	(a == 0)?user_data.arc[0]:user_data.arc[1]);
      exit (15);
    }
    res = gluearcs_or_pinchneck (sketch, a, a2, user_data.arcl[0],
            user_data.arcl[1], ori);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    if (res == 0) {
      fprintf (stderr, "Cannot %s!\n", (ori > 0)?"pinch neck":"glue arcs");
      exit(14);
    }
    break;

    case ACTION_ADDSPHERE:
    ori = 1;
    case ACTION_PUNCHHOLE:
    if (ori != 1) ori = -1;
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    if (user_data.mrnum < 1 || 
        (user_data.stnum < 1 && sketch->huffman_labelling)) {
      fprintf (stderr, "You must specify a region using option -r\n");
      fprintf (stderr, "   if huffman a stratum is required, use -r x:y\n");
      list_strata (sketch);
      exit(15);
    }
    r = findregion (sketch, user_data.region[0]);
    if (r == 0) {
      fprintf (stderr, "Cannot find region %d\n", user_data.region[0]);
      exit(15);
    }
    res = add_s1 (sketch, r, user_data.stratum[0], ori);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    if (res == 0) {
      fprintf (stderr, "Cannot add s1!\n");
      exit(14);
    }
    break;

    case ACTION_WRAP:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    res = put_in_s1 (sketch);
    make_extregion_first (sketch);
    globals.dorecomputef = 0;
    if (docanonify) postprocesssketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    if (res == 0) {
      fprintf (stderr, "Cannot wrap contour in an s1!\n");
      exit(14);
    }
    break;

    case ACTION_UNION:
    case ACTION_CONNECTEDSUM:
    case ACTION_KNOTSUM:
    if ((s1 = readcontour (infile)) == 0) exit (14);
    canonify (s1);
    infile = new_file (infile, infile2);
    if ((s2 = readcontour (infile)) == 0) exit (14);
    canonify (s2);
    switch (action)
    {
      case ACTION_UNION:
      res = sketch_union (s1, s2);   // this adds s2 to s1
      break;

      case ACTION_CONNECTEDSUM:
      res = sketch_sum (s1, s2);     // this sums s2 to s1
      break;

      case ACTION_KNOTSUM:
      res = sketch_knotsum (s1, s2); // sum as knotted tubes
      break;
    }
    if (docanonify) canonify (s1);
    printsketch (s1);
    if (res == 0) exit (15);
    break;

    case ACTION_EXTRACTCC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    //res = extract_connected_component (*ccids, sketch);
    res = extract_connected_components (ncc, ccids, sketch);
    printsketch (sketch);
    if (res == 0) exit (15);
    break;

    case ACTION_REMOVECC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    count = count_connected_components (sketch);
    ccid_isvalidp (ncc, ccids, count);
    assert (ncc == 1);
    res = remove_connected_component (*ccids, sketch);
    printsketch (sketch);
    if (res == 0) exit (15);
    break;

    case ACTION_CCORIENTATION:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    count = count_connected_components (sketch);
    ccid_isvalidp (ncc, ccids, count);
    assert (ncc == 1);
    res = connected_component_orientation (*ccids, sketch);
    if (res == 0) exit (15);
    if (quiet) printf ("%c\n", (res == 1)?'+':'-');
      else printf ("Component %d is %s oriented\n", *ccids+1,
                    (res == 1)?"positively":"negatively");
    break;

    case ACTION_FINDCCPARENT:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    count = count_connected_components (sketch);
    ccid_isvalidp (ncc, ccids, count);
    assert (ncc == 1);
    res = find_connected_component_parent (*ccids, sketch);
    res++;
    if (quiet) printf ("%d\n", res);
      else
    {
      if (res <= 0) printf ("None\n");
      else printf ("Contained in component %d\n", res);
    }
    break;

    case ACTION_CCCHILDS:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    //count = count_connected_components (sketch);
    //ccid_isvalidp (ncc, ccids, count);
    assert (ncc == 1);
    print_connected_component_childs (*ccids, sketch);
    break;

    case ACTION_CCORDERING:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    print_connected_components_ordering (sketch);
    break;

    case ACTION_COUNTCC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    count = count_connected_components (sketch);
    if (quiet) printf ("%d\n", count);
      else printf ("Connected components: %d\n", count);
    break;

    case ACTION_ISCONTOUR:
    case ACTION_ISHUFFMAN:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    res = appcontourcheck (sketch, 
          (action == ACTION_ISHUFFMAN)?1:0, 
          (quiet)?0:1);
    if (res == 0) exit (15);
    break;

    case ACTION_CANONIFY:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_GIOVECANONIFY:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    giovecanonify (sketch);
    if (renumber) giovepostcanonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_COMPARE:
    s1 = readcontour (infile);
    canonify (s1);
    infile = new_file (infile, infile2);
    s2 = readcontour (infile);
    canonify (s2);
    res = sketchcmp (s1, s2);
    switch (res)
    {
      case 0: printf ("s1 = s2\n");
      break;
      case 1: printf ("s1 > s2\n");
      break;
      case -1: printf ("s1 < s2\n");
      res = 2;
      break;
    }
    exit (res);
    break;

    case ACTION_KNOT2MORSE:
    knot2morse (infile);
    break;

    case ACTION_ANY2MORSE:
    any2morse (infile);
    break;

    case ACTION_CODE:
    case ACTION_DTCODE:
    case ACTION_RDTCODE:
    case ACTION_GAUSSCODE:
    loiv = dtorgausscodefromfile (infile);

    switch (loiv->type)
    {
      case LOIV_ISDTCODE:
      case LOIV_ISRDTCODE:
      if (action == ACTION_GAUSSCODE || action == ACTION_RGAUSSCODE)
      {
        newloiv = dtcode2gausscode (loiv);
        if (newloiv)
        {
          freeloiv (loiv);
          loiv = newloiv;
        }
      }
      break;

      case LOIV_ISGAUSSCODE:
      case LOIV_ISRGAUSSCODE:
      if (action == ACTION_DTCODE || action == ACTION_RDTCODE)
      {
        newloiv = gausscode2dtcode (loiv);
        if (newloiv)
        {
          freeloiv (loiv);
          loiv = newloiv;
        }
      }
      break;
    }
 
    if (action == ACTION_RDTCODE && loiv->type == LOIV_ISDTCODE) realize_loiv (loiv);
    printloiv (loiv);
    break;

    case ACTION_WIRTINGER:
    loiv = dtorgausscodefromfile (infile);
    if (loiv->type == LOIV_ISDTCODE) realize_loiv (loiv);
    assert (loiv->type == LOIV_ISDTCODE || loiv->type == LOIV_ISRDTCODE);
    p = wirtingerfromloiv (loiv);
    if (p) print_presentation (p);
    break;

    case ACTION_PRINTMORSE:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    printmorse (sketch);
    break;

    case ACTION_SUGGEST_P_SURGERY:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    if (suggest_p_surgery (sketch, &r, &i))
    {
      if (!quiet)
      {
        printf ("punchhole surgery on region %d, strata %d-%d does not change fundamental group", r->tag, i, i+1);
        printf ("  of the inside solid set:\n\n");
      }
      printf ("contour punchhole -r %d -s %d\n", r->tag, i);
    } else {
      if (!quiet) printf ("Cannot find any suitable punchhole surgery\n");
      exit (13);
    }

    break;

    case ACTION_CHARACTERISTIC:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (fg_type == FG_SURFACE && globals.viacc == 0)
    {
      if (quiet) printf ("%d\n", euler_characteristic (sketch));
        else printf ("Euler characteristic: %d\n", euler_characteristic (sketch));
      break;
    }
    if (appcontourcheck (sketch, 1, 0) == 0)
    {
      fprintf (stderr, "This sketch is NOT a labelled apparent contour.\n");
      exit (13);
    }
    if (sketch->isempty)
    {
      fprintf (stderr, "This sketch is EMPTY.\n");
      exit (13);
    }
    ccomplex = compute_cellcomplex (sketch, fg_type);
    res = complex_characteristic (ccomplex);
    if (quiet) printf ("%d\n", complex_characteristic (ccomplex));
    else {
      switch (fg_type)
      {
        case FG_SURFACE:
          printf ("Euler characteristic (via cell complex): %d\n", res);
        break;

        case FG_INTERNAL:
          printf ("Euler characteristic of inside (via cell complex): %d\n", res);
        break;

        case FG_EXTERNAL:
          printf ("Euler characteristic of outside (via cell complex): %d\n", res);
        break;
      }
    }
    break;

    case ACTION_INFO:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    showinfo (sketch);
    if (user_data.mrnum > 0) {
      for (i = 0; i < user_data.mrnum; i++) {
	r = findregion (sketch, user_data.region[i]);
	if (r)
          printf ("Marked region: %d\n", r->tag);
	else
          printf ("Cannot find region with tag %d\n",
	      user_data.region[i]);
      }
    }
    if (user_data.manum > 0) {
      for (i = 0; i < user_data.manum; i++) {
	a = findarc (sketch, user_data.arc[i]);
	if (a)
          printf ("Marked arc: %d:%d\n", a->tag, user_data.arcl[i]);
	else
          printf ("Cannot find arc with tag %d\n", 
	      user_data.arc[i]);
      }
    }
    break;

    case ACTION_FRONTBACK:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    frontback (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_LEFTRIGHT:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    leftright (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_EVERT:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    changeextregion (sketch, newextregion);
    globals.dorecomputef = 0;
    if (docanonify) postprocesssketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_3DEVERT:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    evert3d (sketch, newextregion);
    globals.dorecomputef = 0;
    if (docanonify) postprocesssketch (sketch);
    if (docanonify) canonify (sketch);
    printsketch (sketch);
    break;

    case ACTION_MENDES:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    mendes = compute_mendes (sketch);
    print_mendes (mendes);
    break;

    case ACTION_CELLCOMPLEX:
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    if (sketch->isempty)
    {
      fprintf (stderr, "Cannot compute cell complex of empty contour.\n");
      exit (13);
    }
    ccomplex = compute_cellcomplex (sketch, fg_type);

    if (globals.simplifycomplex) {
      count = complex_collapse (ccomplex);
      if (debug) printf ("Collapsed %d cell pairs\n", count);
    }
    if (debug) printf ("Euler characteristic %d = %d nodes - %d arcs + %d faces\n",
                        ccomplex->nodenum - ccomplex->arcnum + ccomplex->facenum,
                        ccomplex->nodenum, ccomplex->arcnum, ccomplex->facenum);
    cellcomplex_print (ccomplex, verbose + (quiet?0:1));
    if (!quiet) {
      count =  find_spanning_tree (ccomplex);
      printf ("Found %d connected components\n", count);
    }
    break;

    case ACTION_FUNDAMENTAL:
    case ACTION_AFUNDAMENTAL:
    case ACTION_FOXJACOBIAN:
    case ACTION_ALEXANDER:
    case ACTION_LINKINGNUMBER:
    case ACTION_DEFICIENCY:
    globals.focus_on_fundamental++;
    tok = gettoken (infile);
    ungettoken (tok);
    if (tok == TOK_FPGROUP || tok == TOK_ALEXANDER || tok == TOK_IDEAL)
    {
      switch (tok)
      {
        case TOK_FPGROUP:
          p = (struct presentation *) malloc (sizeof (struct presentation));
          p->rules = p->elements = 0;
          read_group_presentation (infile, p);
          switch (action)
          {
            case ACTION_FUNDAMENTAL:
              fundamental_group (p);
              break;
            case ACTION_AFUNDAMENTAL:
              abelianized_fundamental_group (p);
              break;
            case ACTION_FOXJACOBIAN:
              if (globals.simplifypresentation) simplify_presentation (p);
              foxjacobian (p);
              break;
            case ACTION_ALEXANDER:
              if (globals.simplifypresentation) simplify_presentation (p);
              res = alexander (p);
              if (res == 0) exit (1);
              break;
            case ACTION_LINKINGNUMBER:
              if (globals.simplifypresentation) simplify_presentation (p);
              linkingnumber (p);
              break;
            case ACTION_DEFICIENCY:
              if (globals.simplifypresentation) simplify_presentation (p);
              print_deficiency (p);
              break;
          }
          break;

        case TOK_ALEXANDER:
        case TOK_IDEAL:
          ai = read_alexander_ideal (infile);
          switch (action)
          {
            case ACTION_ALEXANDER:
              if (globals.simplifyideal)
              {
                ai = laurent_simplify_ideal (ai);
              }
              alexander_fromideal (ai);
            break;

            case ACTION_LINKINGNUMBER:
              linkingnumber_fromideal (ai);
            break;

            default:
              printf ("Invalid action %d for Alexander ideal input file\n", action);
              exit (12);
            break;
          }
          break;

        default:
          printf ("Reading an input file with token (TOK_) %d (see parser.h) is not implemented\n", tok);
          exit (11);
          break;
      }
    } else {
      if ((sketch = readcontour (infile)) == 0) exit (14);
      if (docanonify) canonify (sketch);
      if (sketch->isempty)
      {
        fprintf (stderr, "Cannot compute fundamental group of EMPTY contour\n");
        exit (13);
      }
      ccomplex = compute_cellcomplex (sketch, fg_type);
      if (globals.simplifycomplex) {
        count = complex_collapse (ccomplex);
        if (debug) printf ("%d pairs of cells collapsed\n", count);
      }
      compute_fundamental (ccomplex, action);
    }
    break;

    case ACTION_CCCOUNTSL2ZP:
    case ACTION_CCCOUNTSN:
    case ACTION_CCCOUNTPSL2Q:
    globals.focus_on_fundamental++;
    tok = gettoken (infile);
    ungettoken (tok);
    switch (tok)
    {
      struct presentationlist *pstlist;
      case TOK_FPGROUP:
      p = (struct presentation *) malloc (sizeof (struct presentation));
      p->rules = p->elements = 0;
      read_group_presentation (infile, p);
      if (globals.simplifypresentation) simplify_presentation (p);
      switch (action)
      {
        case ACTION_CCCOUNTSL2ZP:
        cccountsl2zp (p);
        break;

        case ACTION_CCCOUNTPSL2Q:
        cccountpsl2q (p);
        break;

        case ACTION_CCCOUNTSN:
        cccountsn (p);
        break;
      }
      break;

      case TOK_FPGROUPLIST:
      pstlist = read_group_presentation_list (infile);
      switch (action)
      {
        case ACTION_CCCOUNTSL2ZP:
        assert (0);
        break;

        case ACTION_CCCOUNTPSL2Q:
        cccountpsl2q_list (pstlist);

        case ACTION_CCCOUNTSN:
        cccountsn_list (pstlist);
        break;
      }
      break;

      default:
      if ((sketch = readcontour (infile)) == 0) exit (14);
      if (docanonify) canonify (sketch);
      if (sketch->isempty)
      {
        fprintf (stderr, "Cannot compute fundamental group of EMPTY contour\n");
        exit (13);
      }
      ccomplex = compute_cellcomplex (sketch, fg_type);
      if (globals.simplifycomplex) {
        count = complex_collapse (ccomplex);
        if (debug) printf ("%d pairs of cells collapsed\n", count);
      }
      compute_fundamental (ccomplex, action);
      break;
    }
    break;

    case ACTION_NEWFEATURE:
    //printf ("Sorry, there is no new feature to experiment with...\n");
    //exit (14);
    if ((sketch = readcontour (infile)) == 0) exit (14);
    if (docanonify) canonify (sketch);
    if (sketch->isempty)
    {
      fprintf (stderr, "Cannot compute fundamental group of EMPTY contour\n");
      exit (13);
    }
    assert (fg_type == FG_SURFACE);
    verbose++;
    debug++;
    ccomplex = compute_cellcomplex (sketch, fg_type);
    //compute_fundamental (ccomplex, action);
    ccnum = find_spanning_tree (ccomplex);
    assert (ccnum == 1);
    cccc = ccomplex->cc;
    cccc->p = compute_fundamental_single (ccomplex, cccc);
    fundamental_group (cccc->p);

    break;

    default:
    printf ("invalid action %d\n", action);
    exit (1);
    break;
  }
  if (infile) fclose (infile);
  exit (0);
}

struct sketch *
readcontour (FILE *file)
{
  int tok;
  struct sketch *s;
  struct vecofintlist *loiv;
#ifdef HAVE_UNISTD_H
  int retcode, status, cpid;
  int pipedes[2];
  FILE *tomorsefile;
#endif

  tok = gettoken (file);

  if (tok == TOK_MORSE) return (readmorse (file));
  if (tok == TOK_SKETCH) return (readsketch (file));
  if (tok == TOK_DTCODE) return (readdtcode (file));
  if (tok == TOK_GAUSSCODE) return (readgausscode (file));
  if (tok == TOK_KNOTSCAPE)
  {
    loiv = readknotscape (file, &s);
    if (s) return (s);
    switch (loiv->type)
    {
      case LOIV_ISDTCODE:
        realize_loiv (loiv);
      case LOIV_ISRDTCODE:
        return (realize_dtcode (loiv->len, loiv->vec, loiv->handedness));
        break;
      case LOIV_ISGAUSSCODE:
        return (readgausscodeloiv (loiv));
        break;
    }
    printf ("FATAL: invalid loiv type: %d\n", loiv->type);
    exit (1);
  }
#ifdef HAVE_UNISTD_H
  if (tok == TOK_KNOT)
  {
    fprintf (stderr, "reading knot description\n");
    retcode = pipe (pipedes);
    cpid = fork ();
    if (cpid < 0) exit (123);
    if (cpid == 0)  /* this is the child */
    {
      /* the child only writes on pipe */
      close (pipedes[0]);
      /* function knot2morse writes on stdout, changing descriptor */
      retcode = dup2 (pipedes[1], 1);
      if (retcode < 0) {perror ("error in dup2"); exit (222);}
      ungettoken (tok);
      retcode = knot2morse (file);
      close (pipedes[1]);
      if (retcode) exit (0);
      exit (333);
    } else {        /* this is the parent */
      /* the parent reads from the pipe */
      close (pipedes[1]);
      tomorsefile = fdopen (pipedes[0],"r");
      tok = gettoken (tomorsefile);
      if (tok != TOK_MORSE)
      {
        fprintf (stderr, "invalid token %d from knot2morse\n", tok);
        fclose (tomorsefile);
        return (0);
      }
      s = readmorse (tomorsefile);
      fclose (tomorsefile);
      waitpid (cpid, &status, 0);
      retcode = WEXITSTATUS (status);
      return (s);
    }
  }
#endif
  fprintf (stderr, "Only 'morse'/'sketch' formats implemented\n");
  exit (2);
}

struct vecofintlist *
dtorgausscodefromfile (FILE *file)
{
  int tok;
  struct vecofintlist *loiv;
  struct sketch *s;

  tok = gettoken (file);

  if (tok == TOK_DTCODE)
  {
    loiv = readvecofintlist (file, LOIV_ISDTCODE);
    return (loiv);
  }

  if (tok == TOK_GAUSSCODE)
  {
    loiv = readvecofintlist (file, LOIV_ISGAUSSCODE);
    return (loiv);
  }

  if (tok == TOK_KNOTSCAPE)
  {
    loiv = readknotscape (file, &s);
    if (s)
    {
      fprintf (stderr, "Cannot recover a dtcode for this object\n");
      exit (2);
    }
    return (loiv);
  }
  fprintf (stderr, "Only 'dtcode'/'knotscape' formats implemented\n");
  exit (2);
}

struct sketch *
readgausscode (FILE *file)
{
  struct vecofintlist *loiv;
  loiv = readvecofintlist (file, LOIV_ISGAUSSCODE);
  return (readgausscodeloiv (loiv));
}

void
ccid_isvalidp (int ncc, int *ccids, int count)
{
  int i;

  for (i = 0; i < ncc; i++)
  {
    if (ccids[i] >= count || ccids[i] < 0)
    {
      fprintf (stderr, "Invalid cc id: %d, valid range:[1,%d]\n",
        ccids[i] + 1, count);
      exit (15);
    }
  }
}

/*
 * open input description file based on argument on command line
 */

FILE *
open_description_file (char *arg, int argnum)
{
  FILE *infile;
  char *subdirs[]={".", "immersed", "knots", "links", "handlebody_knots", 0};
  char *exts[]={"morse", "sketch", "knot", "dtcode", 0};
  char internalname1[]="shortcut_internal.knotname";
  char internalname2[]="shortcut_internal2.knotname";
  int subdirid, extid, len;

  strncpy (examplesfilename, arg, MAXFILELENGTH);
  infile = fopen (arg, "r");
  if (infile) return (infile);
  if (! isalnum (arg[0]) || EXAMPLES_DIR[0] == 0 )
  {
    perror ("Cannot open input file");
    exit (10);
  }
  for (subdirid = 0; subdirs[subdirid]; subdirid++)
  {
    strncpy (examplesfilename, EXAMPLES_DIR, MAXFILELENGTH);
    strncat (examplesfilename, "/", MAXFILELENGTH);
    if (strcmp(subdirs[subdirid],".") != 0)
    {
      strncat (examplesfilename, subdirs[subdirid], MAXFILELENGTH);
      strncat (examplesfilename, "/", MAXFILELENGTH);
    }
    strncat (examplesfilename, arg, MAXFILELENGTH);
    infile = fopen (examplesfilename, "r");
    if (infile && quiet == 0) fprintf (stderr, "Reading from file %s\n", examplesfilename);
    if (infile) return (infile);

    strncat (examplesfilename, ".", MAXFILELENGTH);
    len = strlen (examplesfilename);
    for (extid = 0; exts[extid]; extid++)
    {
      examplesfilename[len] = 0;
      strncat (examplesfilename, exts[extid], MAXFILELENGTH);
      infile = fopen (examplesfilename, "r");
      if (infile && quiet == 0) fprintf (stderr, "Reading from file %s\n", examplesfilename);
      if (infile) return (infile);
    }
  }
  if (globals.knotname_fallback)
  {
    char knchar;
    if ((knchar = check_basic_knotname (arg)) != 0)
    {
      switch (knchar)
      {
        case 'K':
        case 'L':
          assert (argnum == 1 || argnum == 2);
          if (argnum == 1) strncpy (globals.knotname1, arg, MAXKNOTNAMELENGTH);
          if (argnum == 2) strncpy (globals.knotname2, arg, MAXKNOTNAMELENGTH);
          strncpy (examplesfilename, "examples/", MAXFILELENGTH);
          strncat (examplesfilename, (argnum == 1)?internalname1:internalname2, MAXFILELENGTH);
	  infile = fopen (examplesfilename, "r");
          if (infile) return (infile);
          strncpy (examplesfilename, EXAMPLES_DIR, MAXFILELENGTH);
          strncat (examplesfilename, "/", MAXFILELENGTH);
          strncat (examplesfilename, (argnum == 1)?internalname1:internalname2, MAXFILELENGTH);
	  infile = fopen (examplesfilename, "r");
          if (infile) return (infile);
          perror ("Anyway, cannot internal input file");
          exit (11);
        case 'H':
          strncpy (examplesfilename, EXAMPLES_DIR, MAXFILELENGTH);
          strncat (examplesfilename, "/handlebody_knots/hk", MAXFILELENGTH);
          strncat (examplesfilename, &arg[2], MAXFILELENGTH);
          strncat (examplesfilename, ".knot", MAXFILELENGTH);
          infile = fopen (examplesfilename, "r");
          if (infile && quiet == 0) fprintf (stderr, "Reading from file %s\n", examplesfilename);
          if (infile) return (infile);
          perror ("Cannot open corresponding handlebody-knot file");
          exit (11);
        case 'h':
          strncpy (examplesfilename, EXAMPLES_DIR, MAXFILELENGTH);
          strncat (examplesfilename, "/handlebody_knots/hl", MAXFILELENGTH);
          strncat (examplesfilename, &arg[2], MAXFILELENGTH);
          strncat (examplesfilename, ".knot", MAXFILELENGTH);
          infile = fopen (examplesfilename, "r");
          if (infile && quiet == 0) fprintf (stderr, "Reading from file %s\n", examplesfilename);
          if (infile) return (infile);
          perror ("Cannot open corresponding handlebody-link file");
          exit (11);
      }
    }
  }
  perror ("Cannot open input file");
  exit (10);
}

/*
 * a few commands require two descriptions, in which case
 * the user may indicate two files instead of a single file
 * containing both descriptions.
 */

FILE *
new_file (FILE *oldfile, FILE *newfile)
{
  int tok;

  if (newfile == 0) return (oldfile);

  tok = gettoken (oldfile);

  if (tok != TOK_EOF)
    fprintf (stderr, "Warning: discarding extra information from first description file\n");

  fclose (oldfile);
  return (newfile);
}

/*
 * read a list of comma-separated integers
 */

int
readintlistsub1 (int nccmax, int *ccids, char *s)
{
  int i = 0;
  char *spt;

  spt = s;
  while (1)
  {
    assert (i < nccmax);
    ccids[i++] = (int) strtol (s, &spt, 10) - 1;
    if (*spt == ',')
    {
      s = spt + 1;
      continue;
    }
    if (s == spt) return (i-1);
    return (i);
  }
}

/*
 * knotname given as argument, check if name is sintactically correct
 */

char
check_basic_knotname (char *arg)
{
  switch (*arg)
  {
    case 'K':
    case 'L':
      if (isdigit (arg[1])) return (*arg);
      return (0);

    case 'H':
      switch (arg[1])
      {
        case 'K':
          if (isdigit (arg[2])) return ('H');
          return (0);

        case 'L':
          if (isdigit (arg[2])) return ('h');
          return (0);

        default:
          return (0);
      }
      return (0);
  }
  return (0);
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
