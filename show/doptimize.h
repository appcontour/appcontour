void doptimize (struct polyline *contour);
int check_u_turn (struct polyline *contour, struct line *line);
int check_plateau (struct polyline *contour, struct line *line);
int check_plateau_back (struct polyline *contour, struct line *line);
int check_cross_turn (struct polyline *contour, struct line *line);
struct vertex *site_occupied (struct polyline *contour, int px, int py);
