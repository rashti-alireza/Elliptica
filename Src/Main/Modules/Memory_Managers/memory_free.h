#include "core_lib.h"
#include "error_handling_lib.h"

void free_2d(void *mem0);
void free_2d_mem(void *mem0, const unsigned long c);
void free_needle(Needle_T *needle);
void free_points(Grid_T *const grid);
void free_func_PtoV(sFunc_PtoV_T **func);
void free_field(Field_T *f);
void empty_field(Field_T *f);
void free_func_Patch2Pdouble(sFunc_Patch2Pdouble_T **func);
void free_v2(Field_T *f);
void free_v(Field_T *f);
void free_info(Field_T *f);
void free_coeffs(Field_T *fld);
void free_db_eqs(sEquation_T **db);
void free_interpolation(Interpolation_T *interp_s);
void free_matrix(Matrix_T *m);
void free_attr(Field_T *f);
