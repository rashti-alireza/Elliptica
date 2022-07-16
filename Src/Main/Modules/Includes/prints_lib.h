#ifndef prints_LIB_H
#define prints_LIB_H
#include "elliptica_system_lib.h"



#define SECTION "#######################"
#define PR_LINE "--------------------------------------------------------------------------"

/* forward declaration structures */
struct GRID_T;
struct PATCH_T;
struct MATRIX_T;

/* print flags */
typedef enum PRINT_T
{
  UNDEFINED_PRINT = 0,
  PRINT_PARAMETERS,
  PRINT_COORDS,
  PRINT_INTERFACES
}Print_T;

/* print field ,this is a general stuct for print purposes. */
typedef struct PR_FIELD_T
{
  const struct GRID_T *grid;
  const struct PATCH_T *patch;
  const char *par;
  const char *folder;
  int cycle;
  double time;
  Uint ng;/* number of group */
  void *group;/* points to a group for printing */
  void *opt_patch;/* options for patch */
  void *opt_field;/* options for field */
  void *vptr;/* general pointer for different purposes */
  void *a;/* a in double or float */
  void *b;/* b in double or float */
  void *c;/* c in double or float */
  void *v;/* v in double or float */
  void *file;/* file */
  void *file2;/* file */
  
  /* some options and flags */
  Uint abc_f: 1;/* if 1 prints also the patches and fields in 
                    // their (a,b,c) coords ((X,Y,Z) coords ); 
                    // but it seems boring, 
                    // because it would be just a bunch of boxes. 
                    // so the default is 0, which means print only
                    // in Cartesian coords.
                    // if 0 does nothing. default is 0 */
  Uint multimesh_f: 1;/* if 1 it makes a master file of 
                          // all patches as a whole grid, 
                          // if 0 does nothing. default is 0 */
  Uint multivar_f : 1;/* if 1 it makes a master file of 
                          // all fields, if 0 does nothing. default is 0 */
  
}Pr_Field_T;

int test_print(const Print_T f);
void pr_parameters(void);
void pr_coords(const struct GRID_T *const grid);
void pr_line(void);
void pr_line_custom(const char c);
void pr_half_line_custom(const char c);
void pr_comment(const char *const comment);
void pr_clock(void);
void pr_interfaces(const struct GRID_T *const grid);
void pr_fields(Pr_Field_T *const pr);
double get_time_sec(void);
void pr_spent_time(const double start,const char *const event);
Pr_Field_T *init_PrField(const struct GRID_T *const grid);
void free_PrField(Pr_Field_T *pr);
double pr_derivatives_DiffByNode(const double *const numc, const double *const anac,const struct PATCH_T *const patch,const char *const prefix);
void pr_matrix(const struct MATRIX_T *const M,const char *const name);
void pr_field_difference(const struct GRID_T *const grid,const char *const fld1,const char *const fld2);
void pr_hdf5_silo(Pr_Field_T *const pr);
int print_fields_3D(const struct GRID_T *const grid,const int iteration,const char *const dir);
void print_fields_0D(const struct GRID_T *const grid,const int iteration,const char *const folder);
void print_fields_1D(const struct GRID_T *const grid,const int iteration,const char *const folder);
void print_fields_2D(const struct GRID_T *const grid,const int iteration,const char *const folder);
double convert_clock(int *const days,int *const hours,
                     int *const minutes,int *const seconds);


#endif


