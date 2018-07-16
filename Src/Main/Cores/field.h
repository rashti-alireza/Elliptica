#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_approximation_lib.h"

Field_T *init_field_3d(const char *const name,Grid_T *const grid);
Field_T *init_field_2d(const char *const name,Grid_T *const grid);
void add_field(Field_T *const f,Grid_T *const grid);
Field_T *get_field_S(const char *const name,Grid_T *const grid);
double *make_coeffs_1d(Field_T *const f,const unsigned dir,const Patch_T *const patch);
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2,const Patch_T *const patch);
double *make_coeffs_3d(Field_T *const f,const Patch_T *const patch);
double *make_coeffs_grid_1d(Field_T *const f,const unsigned dir);
double *make_coeffs_grid_2d(Field_T *const f,const unsigned dir1,const unsigned dir2);
double *make_coeffs_grid_3d(Field_T *const f);
static unsigned IsAvailable_1d(Field_T *const f,const Patch_T *const patch,const unsigned dir);
static unsigned IsAvailable_2d(Field_T *const f,const Patch_T *const patch,const unsigned dir1,const unsigned dir2,int *dir);
static unsigned IsAvailable_3d(Field_T *const f,const Patch_T *const patch);
static unsigned Is3d_fft(const Collocation_T *collocation,const Basis_T *basis);
static double *find_1d_coeffs_in_patch(Field_T *const f,const Patch_T *const patch,const unsigned dir);
static void coeffs_patch_Tn_Extrema_1d(const Field_T *const f,const Patch_T *const patch,const unsigned dir,double *const coeffs);
static void add_Tinfo(Field_T *const f,const unsigned pn,const unsigned dir,const Collocation_T collocation,const Basis_T basis);
