/*
 * definitions for the hacon graph.  The "next" pointer is
 * used whenever the surface is not connected.
 */

//struct hacon_strata {
//    int hacontag;
//  };

struct hacongraph {
    struct haconnode *node;
    struct haconarc *arc;
    struct haconnode **nodealloc;
    struct haconarc **arcalloc;
    struct hacongraph *next;
  };

struct haconnode {
    int flag;
    int strato;
    struct region *region;
  };

struct haconarc {
    int flag;
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

void describe_hacon_nodes (int num, int **data, struct sketch *s);

void print_hacon (struct hacongraph *h);
