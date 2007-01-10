/*
 * there will be more flags in the future
 */

struct grflags {
  int onlyvisible;
  int xfigspecial;
  };

extern struct grflags grflags;

/* prototypes */

void grparser (int *argcpt, char *argv[]);
