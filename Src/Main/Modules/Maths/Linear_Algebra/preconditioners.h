#include "core_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "memory_managing_lib.h"
#include "maths_linear_algebra_lib.h"

Matrix_T *invert_matrix(Matrix_T *const M);
void precondition(Matrix_T *const a,double *const b);
//static void precondition_jacobi(Matrix_T *const a,double *const b);
static void precondition_GS(Matrix_T *const A,double *const b);
