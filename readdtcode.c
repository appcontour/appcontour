#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "contour.h"
#include "parser.h"

extern int debug;
extern int verbose;

/*
 * the code vector is twice as the number of nodes, we allow
 * for each odd/even label
 */

#define ALG_DT 1
#define MAXDTCODELEN 200

static char *rolfsen_to_dt[] = {
#include "rolfsen_to_dt.h"
};

static int numnodes;
static int numlabels;
static int numregions;
static int *abscode;
static int *dtsign;
static int *regionsign;
static int *tagged;
static int *marknodes;

struct sketch *realize_dtcode (int numnodes, int *vecofint, int *gregionsign);
int reconstruct_sign (int which, int *gregionsign);
int nextlabel (int label);
int prevlabel (int label);
int inherit (int label);
int isinpath (int i, int start, int end);
int propagate (int sign, int label);
void display_arcs_from_arcs (struct sketch *s);
void display_arcs_from_nodes (struct sketch *s);
void display_regions (struct sketch *s);
void display_regions_from_arcs (struct sketch *s);
void display_regions_from_nodes (struct sketch *s);
int isconsistent (void);
int tour_of_region (int label, int velocity);
void walk_left (int *labelpt, int *velocitypt);

#ifdef ALG_DT
/*
 * a few preliminary definitions.  See at end for the main
 * function, based on Jim Hoste code
 */

extern int experimental;
//typedef unsigned char   Boolean;

static void dt_realize (int *anInvolution, int *aRealization, int aNumCrossings);

void uFatalError (char *function, char *file)
{ printf ("FATAL error in function %s, file %s\n", function, file); exit (3); }
#endif

struct sketch *
readdtcode (FILE *file)
{
  struct sketch *sketch;
  int i;
  int startwithlbracket = 1;
  int tok, *vecofint, *gregionsign;

  /*
   * read dt code
   */

  tok = gettoken (file);
  if (tok != TOK_LBRACE)
  {
    fprintf (stderr, "Error: left brace expected\n");
    return (0);
  }
  tok = gettoken (file);
  if (tok != TOK_LBRACKET)
  {
    startwithlbracket = 0;
    ungettoken (tok);
  }
  vecofint = (int *) malloc (MAXDTCODELEN * sizeof (int));
  gregionsign = (int *) malloc (MAXDTCODELEN * sizeof (int));
  i = 0;
  while ((tok = gettoken (file)) == ISNUMBER || tok == TOK_MINUS)
  {
    if (tok == TOK_MINUS)
    {
      tok = gettoken (file);
      assert (tok == ISNUMBER);
      vecofint[i] = - gettokennumber ();
    } else vecofint[i] = gettokennumber ();
    tok = gettoken (file);
    gregionsign[i] = 0;
    if (tok == KEY_GT || tok == KEY_LT)
    {
      /* handedness is given by the user */
      if (tok == KEY_GT) gregionsign[i] = 1;
      else gregionsign[i] = -1;
    } else ungettoken (tok);
    i++;
    if (i >= MAXDTCODELEN)
    {
      printf ("Error: dtcode exceeds maximum allowed length: %d\n", MAXDTCODELEN);
      free (vecofint);
      return (0);
    }
    if ((tok = gettoken (file)) != TOK_COMMA)
    {
      ungettoken (tok);  /* comma separation is optional */
    }
  }
  if (startwithlbracket && tok != TOK_RBRACKET)
  {
    printf ("Error: missing terminating ]\n");
    free (vecofint);
    return (0);
  }
  if (startwithlbracket) tok = gettoken (file);
  if (tok != TOK_RBRACE)
  {
    printf ("Expected terminating }\n");
    free (vecofint);
    return (0);
  }

  sketch = realize_dtcode (i, vecofint, gregionsign);
  free (vecofint);
  return (sketch);
}

/*
 * main function that reconstructs the correct crossings handedness
 */

struct sketch *
realize_dtcode (int lnumnodes, int *vecofint, int *gregionsign)
{
  struct sketch *sketch;
  int i, numconsistent, curoddnode, curevennode, rcode;
#ifdef ALG_DT
  int *dt_realization;
  int *dt_involution;
  int dt_incomplete;
  int agree, agreeneg;
  char ch;
#endif

  numnodes = lnumnodes;
  numlabels = 2*numnodes;
  assert (numnodes >= 2);

  abscode = (int *) malloc ((numlabels+1)*(sizeof (int)));
  dtsign = (int *) malloc ((numlabels+1)*(sizeof (int)));
  regionsign = (int *) malloc ((numlabels+1)*(sizeof (int)));
  tagged = (int *) malloc ( (2*numnodes+1) * sizeof(int) );
  marknodes = (int *) malloc ( (2*numnodes+1) * sizeof(int) );

  curoddnode = 1;
  for (i = 0; i < numnodes; i++)
  {
    rcode = vecofint[i];
    if ((rcode/2)*2 != rcode)
    {
      printf ("Only even values allowed!\n");
      exit (4);
    }
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
  for (i = 1; i <= numlabels; i++)
  {
    if (abscode[i] == 0)
    {
      printf ("Must use all even number from 2 to %d\n", numlabels);
      exit (5);
    }
    if (abs(i - abscode[i]) <= 1)
    {
      printf ("No tight loop allowed: %d %d\n", i, abscode[i]);
      exit (6);
    }
  }

  numregions = 2 + numlabels - numnodes;

#ifdef ALG_DT
  /* this function wants a zero-based vector */
  dt_realization = (int *) malloc (numlabels * sizeof (int));
  dt_involution = (int *) malloc (numlabels * sizeof (int));
  for (i = 0; i < numlabels; i++) {dt_involution[i] = abscode[i+1] - 1;}
  dt_realize (dt_involution, dt_realization, numnodes);
  dt_incomplete = 0;
  for (i = 0; i < numlabels; i++) if (dt_realization[i] == 0) dt_incomplete = 1;
  if (gregionsign)
  {
    /*
     * this (partial?) reconstruction must be consistent with possible informations given
     * by the user
     */
    agree = 1;
    agreeneg = 1;
    for (i = 0; i < numnodes; i++)
    {
      if (gregionsign[i]*dt_realization[2*i] > 0) agreeneg = 0;
      if (gregionsign[i]*dt_realization[2*i] < 0) agree = 0;
    }
    if ((agree == 0) && (agreeneg == 0))
    {
      printf ("Fatal: the reconstruction is inconsistent with user-supplied values.\n");
      if (verbose)
      {
        for (i = 0; i < numnodes; i++)
        {
          printf ("%d%c ", dtsign[2*i+1]*(dt_involution[2*i]+1), (dt_realization[2*i] > 0)?'>':'<');
        }
        printf ("\n");
      }
      exit (4);
    }
    if (agree == 0)
    {
      for (i = 0; i < numlabels; i++) dt_realization[i] = -dt_realization[i];
    }
    if (verbose)
    {
      for (i = 0; i < numnodes; i++)
      {
        ch = '?';
        if (dt_realization[2*i] > 0) ch = '>';
        if (dt_realization[2*i] < 0) ch = '<';
        printf ("%d%c ", dtsign[2*i+1]*(dt_involution[2*i]+1), ch);
      }
      printf ("\n");
    }

  }
  if (dt_incomplete)
  {
    /* fallback to my old code */
    numconsistent = reconstruct_sign (0, gregionsign);
    if (numconsistent <= 0)
    {
      printf ("No consistent completions found. Perhaps this is a composite knot\n");
      exit (3);
    }
    if (numconsistent > 2)
    {
      printf ("More than one (%d) consistent completion found!\n", numconsistent/2);
      printf ("maybe this is a compound knot\n");
      exit (2);
    }
    reconstruct_sign (1, gregionsign);
  } else {
    for (i = 0; i < numlabels; i++) regionsign[i+1] = dt_realization[i];
  }
  free (dt_involution);
  free (dt_realization);

#else

  numconsistent = reconstruct_sign (0, gregionsign);
  if (numconsistent <= 0)
  {
    printf ("No consistent completions found. Perhaps this is a composite knot\n");
    exit (3);
  }
  if (numconsistent > 2)
  {
    printf ("More than one (%d) consistent completion found!\n", numconsistent/2);
    printf ("maybe this is a compound knot\n");
    exit (2);
  }
  reconstruct_sign (1, gregionsign);
#endif

  sketch = newsketch ();
  // printf ("numregions: %d\n", numregions);
  //printf ("#\n# tubular knot with Dowker-Thistletwait notation\n");
  //printf ("# [");
  //for (i = 1; i <= numlabels; i += 2)
  //{
  //  if (i > 1) printf (", ");
  //  printf ("%d", dtsign[i]*abscode[i]);
  //}
  //printf ("]\n#\nsketch {\n");
  display_arcs_from_arcs (sketch);
  display_arcs_from_nodes (sketch);
  display_regions (sketch);
  display_regions_from_arcs (sketch);
  display_regions_from_nodes (sketch);
  // printf ("}\n");
  free (tagged);
  free (regionsign);
  free (dtsign);
  free (abscode);
  sketch->huffman_labelling = 1;
  if (sketch->regions->next == 0) {
    fprintf (stderr, "Warning: empty sketch!\n");
    sketch->isempty = sketch->huffman_labelling = 1;
  }
  if (debug) printsketch (sketch);
  postprocesssketch (sketch);
  return (sketch);
}

/*
 * Read into the knotscape pak files.
 * Please note: the various *.pak files have been created by Morwen Thistlethwaite and Jim Hoste
 * by use with the well-known software knotscape.
 * They are part of the 'knotscape' program.
 *
 * I would like to thank them for their work!
 */

#define MAXFILELENGTH 2000

static char tokenword[80];
static char pathname[MAXFILELENGTH];

struct sketch *
readknotscape (FILE *file)
{
  extern int quiet, verbose;
  struct sketch *sketch;
  char *knotscape_homes[]={".", "/home", "/usr/local", "/usr/local/share", "/usr/local/share/appcontour", 0};
  char *pakpaths[]={"knotTable", "knotscape/knotTable", "knotscape/knotscape_1.01/knotTable", 0};
  char basename[20];
  int nodeid, i, tok, ip1, ip2, crossings, codelen, alternate, knotnum;
  int sign, ch, tailstart, high, low;
  int sum1 = 0;
  int sum2 = 0;
  int found = 0;
  struct stat s = {0};
  char *namept;
  int *dtcode, *dtpositive;
  FILE *pakfile;

  for (ip1 = 0; knotscape_homes[ip1]; ip1++)
  {
    for (ip2 = 0; pakpaths[ip2]; ip2++)
    {
      strncpy (pathname, knotscape_homes[ip1], MAXFILELENGTH);
      strncat (pathname, "/", MAXFILELENGTH);
      strncat (pathname, pakpaths[ip2], MAXFILELENGTH);
      if (!stat(pathname, &s))
      {
        if (s.st_mode & S_IFDIR) {found = 1; break;}
      }
    }
    if (found) break;
  }
  if (!found)
  {
    printf ("Cannot find knotscape installation.  You should install the knotscape package by\n");
    printf ("Morwen Thistlethwaite and Jim Hoste in /home/knotscape or /usr/local/knotscake\n");
    printf ("The pak files should reside in a folder with a name like /home/knotscape/knotscape_1.01/knotTable/\n");
    exit (2);
  }

  tok = gettoken (file);
  assert (tok == TOK_LBRACE);

  if (getword (file, tokenword, 80) == TOK_EOF) return (0);
  for (i = 0; rolfsen_to_dt[i]; i += 2)
  {
    if (strcmp (tokenword, rolfsen_to_dt[i]) == 0)
    {
      strcpy (tokenword, rolfsen_to_dt[i+1]);
      if (*tokenword == '#')
      {
        printf ("# Warning: there is a difference in the numbering of the last four knots in Rolfsen table\n");
        printf ("# with respect to other references, e.g. the knot atlas.  To avoid this warning you should\n");
        printf ("# postfix the knot name with 'R' (Rolfsen book), or 'KA' (Knot Atlas): e.g. 10_165KA refers\n");
        printf ("# to the last knot in the Knot Atlas and is the same as 10_166R\n");
        strcpy (tokenword, rolfsen_to_dt[i+1]+1);
      }
      if (!quiet)
      {
        printf ("# Converting Rolfsen notation into dt numbering: %s -> %s\n", rolfsen_to_dt[i], tokenword);
      }
    }
  }

  if (!quiet) printf ("# ['pak' files courtesy of Morwen Thistlethwaite and Jim Hoste, in knotscape package]\n");
  //printf ("Knot name: %s\n", tokenword);

  crossings = strtol (tokenword, &namept, 10);
  assert (crossings >= 3 && crossings <= 16);
  switch (*namept++)
  {
    case 'a':
    alternate = 1;
    break;

    case 'n':
    alternate = 0;
    break;

    default:
    printf ("Invalid character after number of crossings: %c (valid values: 'a' or 'n')\n", *namept);
    exit (2);
    break;
  }
  assert (*namept++ == '_');
  assert (isdigit(*namept));
  knotnum = strtol (namept, &namept, 10);

  tok = gettoken (file);
  assert (tok == TOK_RBRACE);

  sprintf (basename, "%d%c", crossings, (alternate)?'a':'n');
  strncat (pathname, "/", MAXFILELENGTH);
  strncat (pathname, basename, MAXFILELENGTH);
  strncat (pathname, ".pak", MAXFILELENGTH);
  //printf ("Split: %d %c %d\n", crossings, alternate?'a':'n', knotnum);
  if (access (pathname, R_OK))
  {
    printf ("Cannot open file %s\n", pathname);
    exit (3);
  }

  dtcode = (int *) malloc (crossings *sizeof (int));
  dtpositive = (int *) malloc (crossings *sizeof (int));

  pakfile = fopen (pathname, "r");
  codelen = crossings - 1;
  for (i = 0; i < crossings; i++) {dtcode[i] = 0; dtpositive[i] = 1;}
  nodeid = 0;
  while ((ch = fgetc (pakfile)) != EOF)
  {
    nodeid++;
    high = ch/16;
    low = ch - 16*high;
    assert (high == low);
    tailstart = codelen - high;
    for (i = tailstart; i < codelen;)
    {
      ch = fgetc (pakfile);
      high = ch/16;
      low = ch - 16*high;
      dtcode[i++] = high;
      if (i >= codelen) break;
      dtcode[i++] = low;
    }
    if (alternate == 0)
    {
      ch = fgetc (pakfile);
      for (i = 7; i >= 0; i--)
      {
        if (i < crossings) dtpositive[i] = ch & 1;
        ch = ch >> 1;
      }
      ch = fgetc (pakfile);
      for (i = 15; i >= 8; i--)
      {
        if (i < crossings) dtpositive[i] = ch & 1;
        ch = ch >> 1;
      }
    }
    if (nodeid == knotnum) {found = 1; break;}
  }
  if (!found)
  {
    printf ("Cannot find knot %d in pak file %s\n", knotnum, pathname);
    exit (2);
  }

  for (i = 0; i < codelen; i++)
  {
    sign = 2*dtpositive[i]-1;
    sum1 += i;
    sum2 += dtcode[i];
    dtcode[i] = sign*(2*dtcode[i] + 2);
  }
  assert (i == codelen);
  sum1 += codelen;
  sign = 2*dtpositive[codelen]-1;
  dtcode[i] = sign*(2*(sum1 - sum2) + 2);
  codelen++;
  free (dtpositive);
  if (verbose)
  {
    printf ("dtcode {[");
    for (i = 0; i < codelen; i++)
    {
      if (i > 0) printf (" ");
      printf ("%d", dtcode[i]);
    }
    printf ("]}\n");
  }
  sketch = realize_dtcode (codelen, dtcode, 0);
  free (dtcode);

  //printf ("found path: %s\n", pathname);
  return (sketch);
}

/*
 * display_arcs_from_arcs
 * node arcs are numbered from 1 to numlabels, each produces two
 * apparent contour arcs, one on the right, same orientation, i -> 2*i - 1
 * one on the left, with opposite orientation, i -> 2*i
 */

void
display_arcs_from_arcs (struct sketch *s)
{
  struct arc *arc1, *arc2;
  int i;

  for (i = 1; i <= numlabels; i++)
  {
    arc1 = newarc (s);
    arc1->tag = 2*i-1;
    if (arc1->next && arc1->tag > arc1->next->tag)
    {
      s->arcs = arc1->next;
      arc1->next = 0;
      insert_arc_in_list (arc1, s->arcs);
    }
    arc2 = newarc (s);
    arc2->tag = 2*i;
    if (arc2->next && arc2->tag > arc2->next->tag)
    {
      s->arcs = arc2->next;
      arc2->next = 0;
      insert_arc_in_list (arc2, s->arcs);
    }
    arc1->depths = (int *) malloc (sizeof (int));
    arc2->depths = (int *) malloc (sizeof (int));
    arc1->depths[0] = arc2->depths[0] = 0;
    arc1->depthsdim = arc2->depthsdim = 1;
    arc1->cusps = arc2->cusps = 0;
    arc1->endpoints = arc2->endpoints = 2;
    //printf ("Arc %d: [0];\n", 2*i - 1);
    //printf ("Arc %d: [0];\n", 2*i);
  }
}

/*
 * display_arcs_from_nodes
 */

void
display_arcs_from_nodes (struct sketch *s)
{
  struct arc *arc1, *arc2, *arc3, *arc4;
  int i, overpass;
  int offset = 2*numlabels;

  //printf ("# start of node arcs\n");
  for (i = 1; i <= numnodes; i++)
  {
    arc1 = newarc (s);
    arc1->tag = offset + 4*i - 3;
    if (arc1->next && arc1->tag > arc1->next->tag)
    {
      s->arcs = arc1->next;
      arc1->next = 0;
      insert_arc_in_list (arc1, s->arcs);
    }
    arc2 = newarc (s);
    arc2->tag = offset + 4*i - 2;
    if (arc2->next && arc2->tag > arc2->next->tag)
    {
      s->arcs = arc2->next;
      arc2->next = 0;
      insert_arc_in_list (arc2, s->arcs);
    }
    arc3 = newarc (s);
    arc3->tag = offset + 4*i - 1;
    if (arc3->next && arc3->tag > arc3->next->tag)
    {
      s->arcs = arc3->next;
      arc3->next = 0;
      insert_arc_in_list (arc3, s->arcs);
    }
    arc4 = newarc (s);
    arc4->tag = offset + 4*i - 0;
    if (arc4->next && arc4->tag > arc4->next->tag)
    {
      s->arcs = arc4->next;
      arc4->next = 0;
      insert_arc_in_list (arc4, s->arcs);
    }
    arc1->depths = (int *) malloc (sizeof (int));
    arc2->depths = (int *) malloc (sizeof (int));
    arc3->depths = (int *) malloc (sizeof (int));
    arc4->depths = (int *) malloc (sizeof (int));
    overpass = 0;
    if (dtsign[2*i-1] < 0) overpass = 2;
    arc1->depths[0] = arc3->depths[0] = 2 - overpass;
    arc2->depths[0] = arc4->depths[0] = overpass;
    arc1->depthsdim = arc2->depthsdim = arc3->depthsdim = arc4->depthsdim = 1;
    arc1->cusps = arc2->cusps = arc3->cusps = arc4->cusps = 0;
    arc1->endpoints = arc2->endpoints = arc3->endpoints = arc4->endpoints = 2;
    //printf ("Arc %d: [%d];\n", offset + 4*i - 3, 2 - overpass);
    //printf ("Arc %d: [%d];\n", offset + 4*i - 2, overpass);
    //printf ("Arc %d: [%d];\n", offset + 4*i - 1, 2 - overpass);
    //printf ("Arc %d: [%d];\n", offset + 4*i - 0, overpass);
  }
}

/*
 * display_regions
 */

void
display_regions (struct sketch *s)
{
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *a;
  int *redtagged, *bluetagged;
  int redregionnum = 0;
  int blueregionnum = 0;
  int i, arc, atag;

  redtagged = (int *) malloc ( (numlabels+1) * sizeof(int) );
  bluetagged = (int *) malloc ( (numlabels+1) * sizeof(int) );
  for (i = 1; i <= numlabels; i++) redtagged[i] = bluetagged[i] = 0;

  //printf ("Displaying red regions...\n");

  for (i = 1; i <= numlabels; i++)
  {
    if (redtagged[i]) continue;
    redregionnum++;  /* found new region */
    blast = 0;
    //printf ("New red region: %d:\n", redregionnum);
    region = newregion (s);
    region->tag = redregionnum - 1;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }
    bl = newborderlist_tail (region);
    if (redregionnum == 1)
    {
      bl->isexternal = 1;
      bl = newborderlist_tail (region);
    }
    //printf ("Region %d: ", redregionnum - 1);
    //if (redregionnum == 1) printf ("() ");
    //printf ("(");
    for (arc = i;;)
    {
      redtagged[arc] = redregionnum;
      b = newborder (bl);
      if (blast == 0)
      {
        bl->sponda = b;
      } else {
        b->next = bl->sponda;
        blast->next = b;
      }
      blast = b;
      b->orientation = -1;

      if ( (arc % 2) == 1 )
      { /* odd arc */
        if (regionsign[arc] > 0)
        {
          arc = nextlabel(abscode[arc]);
          atag = 2*arc;
          //printf ("-a%d ", 2*arc);
        } else {
          arc = abscode[arc];
          atag = 2*arc - 1;
          //printf ("-a%d ", 2*arc - 1);
        }
      } else {  /* even arc */
        if (regionsign[prevlabel(arc)] > 0)
        {
          arc = abscode[prevlabel(arc)];
          atag = 2*arc - 1;
          //printf ("-a%d ", 2*arc - 1);
        } else {
          arc = nextlabel(abscode[prevlabel(arc)]);
          atag = 2*arc;
          //printf ("-a%d ", 2*arc);
        }
      }
      for (a = s->arcs; a; a = a->next)
      {
        if (a->tag == atag)
        {
          b->info = a;
          break;
        }
      }
      if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

      if (arc == i) break;
      //{
      //  printf (");\n");
      //  break;
      //}
    }
  }

  //printf ("Displaying blue regions (counterclockwise)...\n");

  for (i = 1; i <= numlabels; i++)
  {
    if (bluetagged[i]) continue;
    blueregionnum++;  /* found new region */
    blast = 0;
    //printf ("New blue region: %d:\n", blueregionnum);
    region = newregion (s);
    region->tag = redregionnum + blueregionnum - 1;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }
    bl = newborderlist_tail (region);
    //printf ("Region %d: ", redregionnum + blueregionnum - 1);
    //printf ("(");
    for (arc = i;;)
    {
      bluetagged[arc] = blueregionnum;
      b = newborder (bl);
      if (blast == 0)
      {
        bl->sponda = b;
      } else {
        b->next = bl->sponda;
        blast->next = b;
      }
      blast = b;
      b->orientation = -1;
      if ( (arc % 2) == 1 )
      { /* odd arc */
        if (regionsign[prevlabel(arc)] > 0)
        {
          arc = abscode[prevlabel(arc)];
          atag = 2*arc - 1;
          //printf ("-a%d ", 2*arc - 1);
        } else {
          arc = nextlabel(abscode[prevlabel(arc)]);
          atag = 2*arc;
          //printf ("-a%d ", 2*arc);
        }
      } else {  /* even arc */
        if (regionsign[arc] > 0)
        {
          arc = nextlabel(abscode[arc]);
          atag = 2*arc;
          //printf ("-a%d ", 2*arc);
        } else {
          arc = abscode[arc];
          atag = 2*arc - 1;
          //printf ("-a%d ", 2*arc - 1);
        }
      }
      for (a = s->arcs; a; a = a->next)
      {
        if (a->tag == atag)
        {
          b->info = a;
          break;
        }
      }
      if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

      if (arc == i) break;
      //{
      //  printf (");\n");
      //  break;
      //}
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
display_regions_from_arcs (struct sketch *s)
{
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *a;
  int i, atag;
  int offset = numregions - 1;
  int arcsoffset = 2*numlabels + 1;
  int arcahead, arcbehind;

  //printf ("# elongated regions corresponding to arcs\n");
  for (i = 1; i <= numlabels; i++)
  {
    blast = 0;
    region = newregion (s);
    region->tag = offset + i;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }

    if ( (i % 2) == 1)
    {
      arcahead = 2*(i-1);
      arcbehind = 2*(abscode[prevlabel(i)] - 1) + 3;
    } else {
      arcahead = 2*(abscode[i]-1) + 1;
      arcbehind = 2*(prevlabel(i) - 1) + 2;
    }
    //printf ("Region %d: ", offset + i);
    bl = newborderlist_tail (region);
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    atag = 2*i - 1;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = -1;
    atag = arcsoffset + arcahead;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    atag = 2*i;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = -1;
    atag = arcsoffset + arcbehind;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);
    //printf ("(+a%d -a%d +a%d -a%d);\n", 2*i - 1, arcsoffset + arcahead, 2*i, arcsoffset + arcbehind);
  }
}

/*
 * display_regions_from_nodes
 */

void
display_regions_from_nodes (struct sketch *s)
{
  struct region *region;
  struct borderlist *bl;
  struct border *b, *blast;
  struct arc *a;
  int i, atag, oddarc, ori;
  int offset = numregions + numlabels;
  int arcsoffset = 2*numlabels - 1;

  //printf ("# small nodal regions\n");
  for (i = 0; i < numnodes; i++)
  {
    blast = 0;
    region = newregion (s);
    region->tag = offset + i;
    if (s->regions->next == region) /* region was inserted in second place */
    {
      /* swap first two positions */
      s->regions->next = region->next;
      region->next = s->regions;
      s->regions = region;
    }
    assert (s->regions == region);
    if (region->next && region->tag > region->next->tag)
    {
      s->regions = region->next;
      region->next = 0;
      insert_region_in_list (region, s->regions);
    }
    bl = newborderlist_tail (region);
    oddarc = 2*i + 1;
    ori = regionsign[oddarc];
    oddarc = 2*oddarc + arcsoffset;

    atag = oddarc;
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    atag = oddarc + 2 - ori;
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    atag = oddarc + 2;
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    atag = oddarc + 2 + ori;
    b = newborder (bl);
    if (blast == 0)
    {
      bl->sponda = b;
    } else {
      b->next = bl->sponda;
      blast->next = b;
    }
    blast = b;
    b->orientation = 1;
    for (a = s->arcs; a; a = a->next)
    {
      if (a->tag == atag)
      {
        b->info = a;
        break;
      }
    }
    if (b->info == 0) printf ("cannot associate arc to a%d\n", atag);

    //printf ("Region %d: (+a%d +a%d +a%d +a%d);\n", offset + i, oddarc,
    //                                                           oddarc + 2 - ori,
    //                                                           oddarc + 2,
    //                                                           oddarc + 2 + ori);
  }
}

/*
 * reconstruct regionsign of nodes
 * if nonzero, gregionsign points to a vector of given handedness:
 * 0: value unassigned
 * +1: if arriving at the odd-numbered label, the crossed strand
 *     is oriented from right to left
 * -1: if arriving at the odd-numbered label, the crossed strand
 *     is oriented from left to right
 *
 * This doesn't work well (see Haken satellite example), perhaps we should try
 * to incorporate the routine "realize" by Jim Hoste in decode_new_DT.c (knotscape)
 */

static int *choices;
static int choicept;

int maximal_expansion (void);
void make_choice (void);

int
reconstruct_sign (int which, int *gregionsign)
{
  int i, numtagged, countreconstructions;

  choices = (int *) malloc ( numnodes*sizeof(int) );
  choices[0] = 0;
  choicept = 0;
  countreconstructions = 0;

  while (1) /* cycle on all possible reconstructions */
  {
    /* make all regionsigns unknown */
    numtagged = 0;
    for (i = 1; i <= numlabels; i++) regionsign[i] = 0;
    if (gregionsign)
    {
      for (i = 0; i < numnodes; i++)
      {
        if (gregionsign[i])
        {
          regionsign[2*i+1] = gregionsign[i];
          regionsign[abscode[2*i+1]] = -gregionsign[i];
          numtagged++;
          numtagged += maximal_expansion ();
        }
      }
    }
    while (numtagged < numnodes) /* subsequent waves of reconstruction */
    {
      make_choice ();
      numtagged++;
      numtagged += maximal_expansion ();
      if (debug) printf ("reconstruction: %d, covered %d of %d\n", countreconstructions, numtagged, numnodes);
    }
    if (isconsistent ())
    {
      countreconstructions++;
      if (countreconstructions == which)
      {
        /* found desired reconstruction */
        free (choices);
        return (countreconstructions);
      }
    }
    /* next possible choice... */
    assert (choices[choicept] == 0);
    if (choicept == 0)
    { /* exhausted all possible choices */
      free (choices);
      return (countreconstructions);
    }
    assert (choicept > 0);
    choicept--;
    while (choices[choicept] == -1)
    {
      choices[choicept] = 0;
      if (choicept == 0)
      { /* exhausted all possible choices */
        free (choices);
        return (countreconstructions);
      }
      choicept--;
    }
    assert (choices[choicept] == 1);
    choices[choicept] = -1;
    choices[choicept+1] = 0;
    choicept = 0;
  }
  free (choices);
  assert (0);
  return (countreconstructions);
}

void
make_choice (void)
{
  int i;

  for (i = 1; i <= numlabels; i++)
  {
    if (regionsign[i] == 0)
    {
      /* found an unoriented node */
      if (choices[choicept] == 0)
      {
        choices[choicept] = 1;
        choices[choicept + 1] = 0;
      }
      assert (abs(choices[choicept]) == 1);
      regionsign[i] = choices[choicept];
      regionsign[abscode[i]] = -choices[choicept];
      choicept++;
      return;
    }
  }
  assert (0);
  return;
}

int
maximal_expansion ()
{
  int insist = 1;
  int totexpansions = 0;
  int i, j, node, expansions;

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
  return (totexpansions);
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
  int i, node_to_tag;

  for (i = 1; i <= numlabels; i++) marknodes[i] = 0;  /* will mark traversed labels */

  while (1)
  {
    node_to_tag = label;
    if (velocity < 0) node_to_tag = prevlabel(label);
    if (marknodes[node_to_tag] != 0) return (0);
    marknodes[node_to_tag] = 1;
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

#ifdef ALG_DT
/*
 * this is taken verbatim from decode_new_DT.c of knotscape
 * based on Jim Hoste previous code
 */

static void dt_realize(
    int     *anInvolution,
    int     *aRealization,
    int     aNumCrossings)
{
    /*
     *  I'm hoping to read through -- and understand -- Dowker
     *  and Thistlethwaite's paper.  But for the moment I'll
     *  try to splice in some of Jim's code and hope I'm interpreting
     *  it correctly.
     */

    int N;
    int i,j;
    int *modTWO_N;
    int *modN;
    int *seq;
    int *emb;
    int *A,*D;
    int Aempty,Dempty;
    int OkSoFar;
    int *phi;
    int x;

    N = aNumCrossings;

    /*
     *  Allocate local arrays.
     */
    modTWO_N    = (int *) malloc(4 * N * sizeof(int));
    modN        = (int *) malloc(2 * N * sizeof(int));
    seq         = (int *) malloc(4 * N * sizeof(int));
    emb         = (int *) malloc(2 * N * sizeof(int));
    A           = (int *) malloc(2 * N * sizeof(int));
    D           = (int *) malloc(2 * N * sizeof(int));
    phi         = (int *) malloc(2 * N * sizeof(int));

    /*create the modTWO_N array*/
    for(i=0;i<2*N;i++){
        modTWO_N[i]=i;
        modTWO_N[i+2*N]=i;
    }
    /*create the modN array*/
    for(i=0;i<N;i++){
        modN[i]=i;
        modN[i+N]=i;
    }
    /* get seq and height from DT code*/
    /* seq is two copies of full DT involution on crossings numbered 0 to 2N-1 */

    /*
     *  Concatenate two copies of theInvolution[] to obtain seq[].
     */
    for(i = 0; i < 2*N; i++)
    {
        seq[i      ] = anInvolution[i];
        seq[i + 2*N] = anInvolution[i];
    }

    /* begin realizability routine to recover embedding of projection */
    /*zero emb, A, and D. A and D will only contain zeroes and ones*/
    for(i=0;i<2*N;i++){
        emb[i]=A[i]=D[i]=0;
    }
    /*set initial conditions*/
    OkSoFar=A[0]=A[seq[0]]=1;
    emb[0]=1;
    emb[seq[0]]=-1;

    /* see if A is empty, ie is all zeroes*/
    for(j=0;j<2*N-1 && !A[j];j++)
        /*nothing*/;
    Aempty=!A[j];

    while(!Aempty && OkSoFar){
        /* let i be least member of A*/
        for(i=0; !A[i]; i++)
            /*nothing*/;
        /*determine phi for this value of i*/
        phi[i]=1;
        for(j=i+1;j<i+2*N;j++){
            phi[modTWO_N[j]]=(seq[j]>=i && seq[j]<=seq[i]) ? -phi[modTWO_N[j-1]]:phi[modTWO_N[j-1]];
        }
        /*establish D*/
        for(j=0; j<i; j++)
            D[j]=1;
        for(j=seq[i]+1; j<2*N; j++)
            D[j]=1;
        /* see if D is empty, ie is all zeroes*/
        for(j=0;j<2*N-1 && !D[j];j++)
            ;
        Dempty=!D[j];
        while(!Dempty && OkSoFar){
            /*let x be least member of D*/
            for(x=0; !D[x]; x++)
                /*nothing*/;
            if(x<i){
                if(seq[x]<i || seq[x]>seq[i]){
                    if(phi[x]*phi[seq[x]]==1){
                        D[x]=D[seq[x]]=0;
                    }
                    else{
                        OkSoFar=0;
                    }
                }
                else{
                    if(emb[x] != 0){/* emb[x] is already defined*/
                        if(phi[x]*phi[seq[x]]*emb[i]==emb[x])
                            D[x]=0;
                        else
                            OkSoFar=0;
                    }
                    else{/* emb[x] is not yet defined. Hence x<>0 and x<>seq[0]*/
                        emb[x]=phi[x]*phi[seq[x]]*emb[i];
                        emb[seq[x]]=-emb[x];
                        D[x]=0;
                        if( modTWO_N[abs(seq[x]-seq[x-1])]==1){
                            /*nothing*/
                        }
                        else{
                            A[x]=A[seq[x]]=1;
                        }
                    }
                }
            }
            else{/*x>seq[i]*/
                if(seq[x]<i || seq[x]>seq[i]){
                    D[x]=D[seq[x]]=0;
                }
                else{
                    if(emb[x]!=0){
                        D[x]=0;
                    }
                    else{
                        emb[x]=phi[x]*phi[seq[x]]*emb[i];
                        emb[seq[x]]=-emb[x];
                        D[x]=0;
                        if( modTWO_N[abs(seq[x]-seq[x-1])]==1){
                            /*nothing*/
                        }
                        else{
                            A[x]=A[seq[x]]=1;
                        }
                    }
                }
            }
            /* see if D is empty, ie is all zeroes*/
            for(j=0;j<2*N-1 && !D[j];j++)
                ;
            Dempty=!D[j];
        }/*end of while*/
        A[i]=0;
        A[seq[i]]=0;
        /* see if A is empty, ie is all zeroes*/
        for(j=0;j<2*N-1 && !A[j];j++)
            /*nothing*/;
        Aempty=!A[j];
    }/*end of while*/

    /*end of realizability routine*/

    if(!OkSoFar){/* sequence is not realizable*/
        uFatalError("realize", "decode_DT");
    }

    /*
     *  Convert emb[] to aRealization[].
     */
    for(i = 0; i < 2*N; i++)
        //aRealization[i] = (emb[i] == +1);
        aRealization[i] = emb[i];

    /*
     *  Free local arrays.
     */
    free(modTWO_N);
    free(modN);
    free(seq);
    free(emb);
    free(A);
    free(D);
    free(phi);
}

#endif
