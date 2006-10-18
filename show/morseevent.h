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
  struct arc *arc;
  struct arc *arc2;
//  struct morseevent *next;
};

void getmorseevent (struct morseevent *mev, int maxone);
void getmorseeventl (struct morseevent *mev);
void getarcinfo (struct morseevent *morseevent);
void getoricusps (int *oript, struct arc **arcpt);

