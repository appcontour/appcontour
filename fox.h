/* prototypes */

void foxjacobian (struct presentation *p);
struct laurentpoly *foxderivative (struct presentationrule *r, int gen, int offset, int indets);
void count_and_map (int *v, int len, int offset, int rank, int *exponvec);
