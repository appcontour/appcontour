/*
 * definitions for the hacon graph.  The "next" pointer is
 * used whenever the surface is not connected.
 */

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
void print_hacon (struct hacongraph *h);
