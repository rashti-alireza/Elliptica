#ifndef maths_linear_algebra_lib
#define maths_linear_algebra_lib


/* matrix storage format */
typedef enum MATRIX_SF_T
{
  UNDEF_SF = 0/* undefined */,
  REG_SF/* regular storage format */,
  CCS_SF/* compressed column storage format */,
  CRS_SF/* compressed row storage format */,
  TRI_SF/* triplet storage format */,
  CCS_L_SF/* long compressed column storage format */,
  CRS_L_SF/* long compressed row storage format */,
  TRI_L_SF/* long triplet storage format */
}Matrix_SF_T;

/* matrix */
typedef struct MATRIX_T
{
  unsigned reg_f: 1;/* regular format */
  unsigned tri_f: 1;/* 1 if tripet storage format, 0 otherwise*/
  unsigned ccs_f: 1;/* 1 if compressed column storage format, 0 otherwise */
  unsigned crs_f: 1;/* 1 if compressed row storage format,0 otherwise */
  unsigned tri_l_f: 1;/* 1 if long tripet storage format, 0 otherwise*/
  unsigned ccs_l_f: 1;/* 1 if long compressed column storage format, 0 otherwise */
  unsigned crs_l_f: 1;/* 1 if long compressed row storage format,0 otherwise */
  long row;
  long col;
  struct/* triplet storage format */
  {
    int *row;
    int *col;
    double *a;
  }tri[1];
  struct/* compressed column storage format */
  {
    int *Ap;
    int *Ai;
    double *Ax;
  }ccs[1];
  struct/* compressed row storage format */
  {
    int *Ap;
    int *Aj;
    double *Ax;
  }crs[1];
  struct/* regular storage format */
  {
    double **A;
  }reg[1];
  struct/* triplet storage format */
  {
    long *row;
    long *col;
    double *a;
  }tri_long[1];
  struct/* compressed column storage format */
  {
    long *Ap;
    long *Ai;
    double *Ax;
  }ccs_long[1];
  struct/* compressed row storage format */
  {
    long *Ap;
    long *Aj;
    double *Ax;
  }crs_long[1];
}Matrix_T;


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



#endif





