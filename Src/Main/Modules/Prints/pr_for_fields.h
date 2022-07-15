#include "core_lib.h"/* it includes prints_lib.h too */
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "pr_for_fields_structs.h"
#include  <unistd.h>

#define STR_LEN (999)

void pr_fields(Pr_Field_T *const pr);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);
void pr_hdf5_silo(Pr_Field_T *const pr);
static void parse_parameter_3d(const char *const par,Pr_Field_T *const pr);
int print_fields_3D(const Grid_T *const grid,const int iteration,const char *const dir);
double print_fields_0D(const Grid_T *const grid,const int iteration,const char *const folder);
void print_fields_1D(const Grid_T *const grid,const int iteration,const char *const folder);
static Field_T **find_field_by_name_or_regex(const Patch_T *const patch,char **const fld_names,Uint *const Nfld);
static double map_to_ref_interval(const double X,const Patch_T *const patch,const double *const min,const double *const max,const int dir, const int map_type);




