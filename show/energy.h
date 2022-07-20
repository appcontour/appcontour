extern double k1_coeff;
extern double k2_coeff;
extern double k3_coeff;
extern double k4_coeff;
extern double k5_coeff;
extern int allowrepulsion;

void energyinit (void);
void tryrepulsiveenergy (struct polyline *contour);
void compute_gradient (struct polyline *contour);
double compute_energy (struct polyline *contour);
void compute_repulsive_gradient (struct polyline *contour);
double compute_repulsive_energy (struct polyline *contour);
double repulsive_field (double dist);
double repulsive_force (double dist);
double repulsive_field2 (double dist);
double repulsive_force2 (double dist);
double get_alpha (struct vertex *a, struct vertex *p, struct vertex *b,
       double *lapt, double *lbpt);
void grad_curv (struct vertex *b, struct vertex *p, double alpha, double lb, double la,
           double *gcxpt, double *gcypt);
