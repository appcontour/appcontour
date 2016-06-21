/* prototypes */

void foxjacobian (struct presentation *p);
void foxderivative (int len, int *rvar, int *weightin, int *weightout, int gen);
struct laurentpoly *map_to_abelian (int *rvec, int *lincomb, int len, int offset, int indets);
int map_to_trivial (int *rvec, int *lincomb, int len);
void count_and_map (int *v, int len, int offset, int rank, int *exponvec);

struct alexanderideal *three_components_link (struct presentation *p);

void print_relator (int *rvar, int len);
void print_groupring_el (int *rvar, int *lincomb, int len);
