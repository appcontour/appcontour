#define ME_TRAN 1
#define ME_TOP 2
#define ME_BOT 3
#define ME_CROSS 4
#define ME_CROSS_NWSE 5
#define ME_CROSS_NESW 6
#define ME_NEWROW 7
#define ME_LASTROW 8
#define ME_RESET 9
#define ME_RBRACE 10

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

