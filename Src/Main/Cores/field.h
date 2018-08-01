#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_approximation_lib.h"

Field_T *add_field(const char *const name,const char *attribute,Patch_T *const patch,const Flag_T alloc_flg);
void remove_field(Field_T *f);
void add_attribute(Field_T *const fld,const char *const attribute);
int LookUpField(const char *const name,const Patch_T *const patch);
double *make_coeffs_1d(Field_T *const f,const unsigned dir);
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2);
double *make_coeffs_3d(Field_T *const f);
static unsigned IsAvailable_1d(Field_T *const f,const unsigned dir);
static unsigned IsAvailable_2d(Field_T *const f,const unsigned dir1,const unsigned dir2,int *dir);
static unsigned IsAvailable_3d(Field_T *const f);
static unsigned Is3d_fft(const Collocation_T *collocation,const Basis_T *basis);
static double *find_1d_coeffs_in_patch(Field_T *const f,const unsigned dir);
static void coeffs_patch_Tn_Extrema_1d(Field_T *const f,const unsigned dir);
static void add_Tinfo(Field_T *const f,const unsigned dir,const Collocation_T collocation,const Basis_T basis);
