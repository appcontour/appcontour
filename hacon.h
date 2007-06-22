/*
 * definitions for the hacon graph.  The "next" pointer is
 * used whenever the surface is not connected.
 */

struct hacongraph {
    int numhaconnodes;
    int numhaconarcs;
    int **nodesdata;
    int *arcsdata;
    int *arcincplus;
    int *arcincminus;
    struct sketch *sketch;
  };

/*
 * prototypes
 */

struct hacongraph *compute_hacon (struct sketch *s);

int tag_hacon_strata (int **data, struct sketch *s);
int single_tag_hacon_strata (int tag, int **data, struct sketch *s);
int hacon_try_expand_node (int tag, int **data, struct sketch *s);
int local_hacon_try_expand_node (int **data, struct border *bp, int k);
int **init_hacon_strata (struct sketch *s);

int tag_hacon_arcs (int *arcdata, struct sketch *s);
int single_tag_hacon_arc (int tag, int *arcdata, struct sketch *s);
int hacon_try_expand_arc (int tag, int *arcdata, struct sketch *s);

void describe_hacon_nodes (int num, int **data, struct sketch *s);
void describe_hacon_arcs (int num, int *arcdata, struct sketch *s);

void print_hacon (struct hacongraph *h);
