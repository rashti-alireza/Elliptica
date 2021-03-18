#include "core_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"

int matrix_by_vector(const Matrix_T *const m, const double *const v,double *const b,const Flag_T flag);
Matrix_T *matrix_by_matrix(const Matrix_T *const a, 
                           const Matrix_T *const b,
                           Matrix_T *const result,
                           const char *const dir);
Matrix_T *CCSOpCCS(Matrix_T *const ccs2,Matrix_T *const ccs1,const char Op);
double read_matrix_entry_ccs(Matrix_T *const m, const long r,const long c);
int matrix_comparison(const Matrix_T *const M1,const Matrix_T *const M2);

