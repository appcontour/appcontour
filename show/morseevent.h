#define ME_TRAN 1
#define ME_TOP 2
#define ME_BOT 3
#define ME_CROSS 4
#define ME_NEWROW 5
#define ME_LASTROW 6
#define ME_RESET 7

struct morseevent {
  int type;
  int ori;
  int ori2;
//  int cusps;
//  int cusps2;
  struct earc *arc;
  struct earc *arc2;
//  struct morseevent *next;
};

void getmorseevent (FILE *filein, struct morseevent *mev, int maxone);
void getmorseeventl (FILE *filein, struct morseevent *mev);
void getarcinfo (FILE *filein, struct morseevent *morseevent);
void getoricusps (FILE *filein, int *oript, struct earc **arcpt);

