#include "core_lib.h"
#include "memory_managing_lib.h"
#include "error_handling_lib.h"
#include "maths_general_lib.h"
#include "maths_linear_algebra_lib.h"

Matrix_T *cast_matrix_ccs(Matrix_T *const m);
Matrix_T *cast_matrix_ccs_long(Matrix_T *const m);
Matrix_T *cast_matrix_reg(Matrix_T *const m);
Matrix_T *compress_stack2ccs(Matrix_T **const S,const unsigned nm,const unsigned *const nr,const unsigned Nrow,const unsigned Ncol,const Flag_T flg);
void copy_reg2reg(const Matrix_T *const reg1,Matrix_T *const reg2);
void copy_ccs2ccs(const Matrix_T *const ccs1,Matrix_T *const ccs2);
void copy_ccs_long2ccs_long(const Matrix_T *const ccs_l1,Matrix_T *const ccs_l2);
static void convert_reg2ccs(const Matrix_T *const reg,Matrix_T *const ccs,const double DropLimit);
static void convert_ccs2reg(const Matrix_T *const ccs,Matrix_T *const reg);
static void convert_ccs_long2reg(const Matrix_T *const ccs,Matrix_T *const reg);
static void convert_reg2ccs_long(const Matrix_T *const reg,Matrix_T *const ccs_l,const double DropLimit);
