/*
 * data types
 */

struct list_of_borders {
  struct border *border;
  struct list_of_borders *next;
};

/*
 * prototypes
 */

void giovecanonify (struct sketch *s);
struct borderlist *giovecanonifyblist (struct borderlist *bl);
struct border *giovecanonifyhole (struct border *entrypoint);
struct borderlist *gioveinsertholeinblist (struct borderlist *blhole, struct borderlist *bl);
void free_list_of_borders (struct list_of_borders *optimalborders);
struct list_of_borders *extract_optimal_borders (struct border *entrypoint);
int giove_compare_holes (struct border *b1, struct border *b2);
int giove_compare_dfs (struct border *b1, struct border *b2);
/* note: it might be the same hole with different entry points */
void giove_normalize (struct border *entrypoint);
void giove_renormalize (struct border *entrypoint);
void giove_normalize_common (struct border *entrypoint, int canonify_holes);
void giove_normalize_dfs (struct border *entrypoint, int canonify_holes);
struct borderlist *gioveinsertholeinblist (struct borderlist *blhole, struct borderlist *bl);
struct border *crossriver (struct border *b);
void reset_marks (struct border *b, short int *marks);
void reset_marks_dfs (struct border *b, short int *marks);
void giove_sort_regions (struct sketch *s);
void giove_relink_regions (struct region *r);
