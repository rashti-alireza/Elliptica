#define SECTION "#######################"
#define PR_LINE "--------------------------------------------------------------------------"

int test_print(const Print_T f);
void pr_parameters(void);
void pr_coords(const Grid_T *const grid);
void pr_line(void);
void pr_comment(const char *const comment);
void pr_clock(void);
void pr_interfaces(const Grid_T *const grid);
void pr_fields(Pr_Field_T *const pr);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);

Flag_T pr_derivatives_DiffByNode(const double *const numc, const double *const anac,const Patch_T *const patch,const char *const prefix);