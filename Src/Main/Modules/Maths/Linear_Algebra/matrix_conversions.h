#include "core_lib.h"
#include "macros_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"

Matrix_T *cast_matrix_ccs(Matrix_T *const m);
static void convert_reg2ccs(const Matrix_T *const reg,Matrix_T *const ccs,const double DropLimit);
