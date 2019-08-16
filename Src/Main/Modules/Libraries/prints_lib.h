#define SECTION "#######################"
#define PR_LINE "--------------------------------------------------------------------------"

struct GIRD_T;
struct PATCH_T;

/* print field */
typedef struct PR_FIELD_T
{
  const struct GIRD_T *grid;
  const struct PATCH_T *patch;
  const char *par;
  const char *folder;
  int cycle;
  double time;
  unsigned ng;/* number of group */
  void *group;/* points to a group for printing */
  void *opt_patch;/* points to options for patch */
  void *opt_field;/* points to options for field */
  void *vptr;/* general pointer for different purposes */
  void *a;/* a in double or float */
  void *b;/* b in double or float */
  void *c;/* c in double or float */
  void *v;/* v in double or float */
  void *file;/* file */
  void *file2;/* file */
}Pr_Field_T;

int test_print(const Print_T f);
void pr_parameters(void);
void pr_coords(const Grid_T *const grid);
void pr_line(void);
void pr_line_custom(const char c);
void pr_half_line_custom(const char c);
void pr_comment(const char *const comment);
void pr_clock(void);
void pr_interfaces(const Grid_T *const grid);
void pr_fields(Pr_Field_T *const pr);
double get_time_sec(void);
void pr_spent_time(const double start,const char *const event);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);
double pr_derivatives_DiffByNode(const double *const numc, const double *const anac,const Patch_T *const patch,const char *const prefix);
void pr_matrix(const Matrix_T *const M,const char *const name);
void pr_field_difference(const Grid_T *const grid,const char *const fld1,const char *const fld2);
