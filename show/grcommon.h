/*
 * there will be more flags in the future
 * these coincide with xfig colors
 */

#define BLACK 0
#define BLUE 1
#define GREEN 2
#define CYAN 3
#define RED 4
#define MAGENTA 5
#define YELLOW 6
#define WHITE 7

struct grflags {
  int onlyvisible;
  int xfigspecial;
  int visiblewidth;
  int invisiblewidth;
  int pointsize;
  int pointcolor;
  int linecolor;
  double dashlength;
  double dotspacing;
  };

extern struct grflags grflags;

/* prototypes */

void grparser (int *argcpt, char *argv[]);
