void doptimize (struct polyline *contour);
int check_u_turn (struct polyline *contour, struct line *line);
int check_flip_left (struct polyline *contour, struct line *line);
int check_plateau (struct polyline *contour, struct line *line, int back);
int check_cross_turn (struct polyline *contour, struct line *line);
int check_cross_slide (struct polyline *contour, struct line *line);
struct vertex *site_occupied (struct polyline *contour, int px, int py);
