#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "maths_general_lib.h"
#include "maths_approximation_lib.h"
#include "coordinates_lib.h"

typedef void Fourier_Transformation_1d_F(double *const values,double *const coeffs,const unsigned n);

Field_T *add_field(const char *const name,const char *attribute,Patch_T *const patch,const Flag_T alloc_flg);
void remove_field(Field_T *f);
void add_attribute(Field_T *const fld,const char *const attribute);
int LookUpField(const char *const name,const Patch_T *const patch);
int LookUpField_E(const char *const name,const Patch_T *const patch);
double *make_coeffs_1d(Field_T *const f,const unsigned dir);
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2);
double *make_coeffs_3d(Field_T *const f);
static unsigned IsAvailable_1d(Field_T *const f,const unsigned dir);
static unsigned IsAvailable_2d(Field_T *const f,const unsigned dir1,const unsigned dir2,int *dir);
static unsigned IsAvailable_3d(Field_T *const f);
static double *find_1d_coeffs_in_patch(Field_T *const f,const unsigned dir);
static void coeffs_patch_Tn_Extrema_1d(Field_T *const f,const unsigned dir);
static void coeffs_patch_Tn_Nodes_1d(Field_T *const f,const unsigned dir);
static void add_Tinfo(Field_T *const f,const unsigned dir,const Collocation_T collocation,const Basis_T basis);
void enable_fields(Grid_T *const grid);
void free_v2(Field_T *f);
void free_attr(Field_T *f);
void free_v(Field_T *f);
void free_info(Field_T *f);
void free_field(Field_T *fld);
void empty_field(Field_T *fld);
void free_coeffs(Field_T *fld);

