#include "core_lib.h"
#include "error_handling_lib.h"
#include "maths_equation_solvings_lib.h"
#include "utilities_lib.h"
#include "maths_linear_algebra_lib.h"

#define DO 1
#define DO_NOT 0

void matrix_tests(void);
Matrix_T *cast_matrix_ccs(Matrix_T *const m);
Matrix_T *cast_matrix_ccs_long(Matrix_T *const m);
Matrix_T *cast_matrix_reg(Matrix_T *const m);
int matrix_by_vector(const Matrix_T *const m, const double *const v,double *const b,const Flag_T flag);
static int cast_matrix_ccs_test(void);
static Matrix_T *make_generic_matrix(const long Nr,const long Nc);
static int read_ccs_test(void);
static int matrices_arithmetic_test(void);
static int cast_matrix_ccs_long_test(void);

