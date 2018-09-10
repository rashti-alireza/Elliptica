#include "core_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"
#include <suitesparse/umfpack.h>

Matrix_T *invert_matrix(Matrix_T *const M);
Matrix_T *cast_matrix_ccs(Matrix_T *const m);
void umfpack_error_di(const double *const Control,const int status,const char *const file,const int line);
void umfpack_error_dl(const double *const Control,const long status,const char *const file,const int line);
static void test_invert_matrix(const Matrix_T *const invertM,const Matrix_T *const M);

