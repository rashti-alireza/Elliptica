#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "maths_solvings_lib.h"

#define DO 1
#define NOT_DO 0

void matrix_tests(void);
Matrix_T *cast_matrix_ccs(Matrix_T *const m);
Matrix_T *cast_matrix_reg(Matrix_T *const m);
static int cast_matrix_ccs_test(void);
static Matrix_T *make_generic_matrix(const long Nr,const long Nc);
static int read_ccs_test(void);

