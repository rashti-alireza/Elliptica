Matrix_T *cast_matrix_ccs(Matrix_T *const m);
Matrix_T *cast_matrix_reg(Matrix_T *const m);
void copy_reg2reg(const Matrix_T *const reg1,Matrix_T *const reg2);
void copy_ccs2ccs(const Matrix_T *const ccs1,Matrix_T *const ccs2);
void matrix_tests(void);
void precondition(Matrix_T *const a,double *const b);
Matrix_T *invert_matrix(Matrix_T *const M);
int matrix_by_vector(const Matrix_T *const m, const double *const v,double *const b,const Flag_T flag);
Matrix_T *matrix_by_matrix(const Matrix_T *const a, const Matrix_T *const b,const char *const dir);
Matrix_T *cast_matrix_ccs_long(Matrix_T *const m);
void copy_ccs_long2ccs_long(const Matrix_T *const ccs_l1,Matrix_T *const ccs_l2);
Matrix_T *compress_stack2ccs(Matrix_T **const S,const unsigned nm,const unsigned *const nr,const unsigned Nrow,const unsigned Ncol,const Flag_T flg);
Matrix_T *CCSOpCCS(Matrix_T *const ccs2,Matrix_T *const ccs1,const char Op);
int matrix_comparison(const Matrix_T *const M1,const Matrix_T *const M2);
Matrix_T *alloc_matrix(const Matrix_SF_T type_e,const long row,const long col);
void free_matrix(Matrix_T *m);




