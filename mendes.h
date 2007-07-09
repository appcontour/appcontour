/*
 * definitions for the mendes graph. 
 */

extern int mendesge;

#define HGE_TEXT 0
#define HGE_PYKIG 1
#define HGE_KIG 2

struct mendesgraph {
    int nummendesnodes;
    int nummendesarcs;
    int **nodesdata;
    int *arcsdata;
    int *arcincplus;
    int *arcincminus;
    int *nodessign;
    int *nodesgenus;
    struct sketch *sketch;
    int *count_strata;
    double *x;
    double *y;
  };

/*
 * prototypes
 */

struct mendesgraph *compute_mendes (struct sketch *s);

int tag_mendes_strata (int **data, struct sketch *s);
int single_tag_mendes_strata (int tag, int **data, struct sketch *s);
int mendes_try_expand_node (int tag, int **data, struct sketch *s);
int local_mendes_try_expand_node (int **data, struct border *bp, int k);
int **init_mendes_strata (struct sketch *s);

int tag_mendes_arcs (int *arcdata, struct sketch *s);
int single_tag_mendes_arc (int tag, int *arcdata, struct sketch *s);
int mendes_try_expand_arc (int tag, int *arcdata, struct sketch *s);

void mendes_compute_cusps_data (struct arc *a, int *dmin, int *dmax, int *necusps);

void describe_mendes_nodes (int num, int **data, struct sketch *s);
void describe_mendes_arcs (int num, int *arcdata, struct sketch *s);

void print_mendes (struct mendesgraph *h);

void mendes_node_canonify (struct mendesgraph *h);
void l_mendes_reorder_n (int *start, int num, struct mendesgraph *h);
int h_node_compare (int tag1, int tag2, struct mendesgraph *h);
int h_count_strata (int tag, struct mendesgraph *h);

void mendes_arc_canonify (struct mendesgraph *h);
void l_mendes_reorder_a (int *start, int num, struct mendesgraph *h);
int h_arc_compare (int tag1, int tag2, struct mendesgraph *h);

void mendes_xy_alloc (struct mendesgraph *h);
void mendes_xy_compute (struct mendesgraph *h);
void mendes_xy_randomize (struct mendesgraph *h);
double mendes_energy (struct mendesgraph *h);
