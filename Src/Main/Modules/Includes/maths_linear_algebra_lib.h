#ifndef maths_linear_algebra_LIB_H
#define maths_linear_algebra_LIB_H
#include "elliptica_system_lib.h"


/* MLA stands for "maths_linear_algebra" */


/* concatenate ex: MLA_CONC(g,U,1,2) => g_U1U2 */
#define MLA_CONC(m,t,x,y)  m##_##t##x##t##y

/* compute inverse of a 3x3 symmetric matrix (back-end) for field
// more info at Matrix_Inverse_3x3_Symmetric_Field. */
#define MLA_COMPUTE_INVERSE_3x3_SYMMETRIC_FIELD(mI,tI,ijk,a00,a01,a02,a10,a11,a12,a20,a21,a22) \
  { \
  MLA_CONC(mI,tI,0,0)[ijk] = (a11*a22 - a12*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,1)[ijk] = (-a01*a22 + a02*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,2)[ijk] = (a01*a12 - a02*a11)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,1,1)[ijk] = a00*(a00*a22 - a02*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,1,2)[ijk] =-a00*(a00*a12 - a02*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,2,2)[ijk] = a00*(a00*a11 - a01*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  }

/* compute inverse of a 3x3 general matrix (back-end) for field
// more info at Matrix_Inverse_3x3_General_Field. */
#define MLA_COMPUTE_INVERSE_3x3_GENERAL_FIELD(mI,tI,ijk,a00,a01,a02,a10,a11,a12,a20,a21,a22) \
  { \
  MLA_CONC(mI,tI,0,0)[ijk] = (a11*a22 - a12*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,1)[ijk] = (-a01*a22 + a02*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,2)[ijk] = (a01*a12 - a02*a11)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,1,0)[ijk] = (-a10*a22 + a12*a20)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);\
  MLA_CONC(mI,tI,1,1)[ijk] = a00*(a00*a22 - a02*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,1,2)[ijk] =-a00*(a00*a12 - a02*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,2,0)[ijk] = (a10*a21 - a11*a20)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);\
  MLA_CONC(mI,tI,2,1)[ijk] =-a00*(a00*a21 - a01*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20));\
  MLA_CONC(mI,tI,2,2)[ijk] = a00*(a00*a11 - a01*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  }

/* compute inverse of a 3x3 symmetric matrix (back-end) for variable (no ijk given)
// more info at Matrix_Inverse_3x3_Symmetric_Var. */
#define MLA_COMPUTE_INVERSE_3x3_SYMMETRIC_VAR(mI,tI,a00,a01,a02,a10,a11,a12,a20,a21,a22) \
  { \
  MLA_CONC(mI,tI,0,0) = (a11*a22 - a12*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,1) = (-a01*a22 + a02*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,2) = (a01*a12 - a02*a11)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,1,1) = a00*(a00*a22 - a02*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,1,2) =-a00*(a00*a12 - a02*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,2,2) = a00*(a00*a11 - a01*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  }

/* compute inverse of a 3x3 general matrix (back-end) for variable (no ijk given)
// more info at Matrix_Inverse_3x3_General_Var. */
#define MLA_COMPUTE_INVERSE_3x3_GENERAL_VAR(mI,tI,a00,a01,a02,a10,a11,a12,a20,a21,a22) \
  { \
  MLA_CONC(mI,tI,0,0) = (a11*a22 - a12*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,1) = (-a01*a22 + a02*a21)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,0,2) = (a01*a12 - a02*a11)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
  MLA_CONC(mI,tI,1,0) = (-a10*a22 + a12*a20)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);\
  MLA_CONC(mI,tI,1,1) = a00*(a00*a22 - a02*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,1,2) =-a00*(a00*a12 - a02*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  MLA_CONC(mI,tI,2,0) = (a10*a21 - a11*a20)/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20);\
  MLA_CONC(mI,tI,2,1) =-a00*(a00*a21 - a01*a20)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20));\
  MLA_CONC(mI,tI,2,2) = a00*(a00*a11 - a01*a10)/((a00*a11 - a01*a10)*(a00*a22 - a02*a20) - (a00*a12 - a02*a10)*(a00*a21 - a01*a20)); \
  }


/* inverse of a 3x3 symmetric matrix, 
// it populates mI for field on each ijk point.
// ijk: the point on the patch
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2.
// mI: matrix inverse name stem, tI: matrix inverse index.
// ex: mI = gI and tI = U => MLA_CONC(mI,tI,0,0) = gI_U0U0. */
#define Matrix_Inverse_3x3_Symmetric_Field(m,t,mI,tI,ijk) \
   MLA_COMPUTE_INVERSE_3x3_SYMMETRIC_FIELD(mI,tI,ijk,\
    MLA_CONC(m,t,0,0)[ijk],MLA_CONC(m,t,0,1)[ijk],MLA_CONC(m,t,0,2)[ijk],\
    MLA_CONC(m,t,0,1)[ijk],MLA_CONC(m,t,1,1)[ijk],MLA_CONC(m,t,1,2)[ijk],\
    MLA_CONC(m,t,0,2)[ijk],MLA_CONC(m,t,1,2)[ijk],MLA_CONC(m,t,2,2)[ijk])

/* inverse of a 3x3 general matrix, 
// it populates mI for field on each ijk point.
// ijk: the point on the patch
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2.
// mI: matrix inverse name stem, tI: matrix inverse index.
// ex: mI = gI and tI = U => MLA_CONC(mI,tI,0,0) = gI_U0U0. */
#define Matrix_Inverse_3x3_General_Field(m,t,mI,tI,ijk) \
   MLA_COMPUTE_INVERSE_3x3_GENERAL_FIELD(mI,tI,ijk,\
    MLA_CONC(m,t,0,0)[ijk],MLA_CONC(m,t,0,1)[ijk],MLA_CONC(m,t,0,2)[ijk],\
    MLA_CONC(m,t,1,0)[ijk],MLA_CONC(m,t,1,1)[ijk],MLA_CONC(m,t,1,2)[ijk],\
    MLA_CONC(m,t,2,0)[ijk],MLA_CONC(m,t,2,1)[ijk],MLA_CONC(m,t,2,2)[ijk])

/* inverse of a 3x3 symmetric matrix, it populates mI variable
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2.
// mI: matrix inverse name stem, tI: matrix inverse index.
// ex: mI = gI and tI = U => MLA_CONC(mI,tI,0,0) = gI_U0U0. */
#define Matrix_Inverse_3x3_Symmetric_Var(m,t,mI,tI) \
   MLA_COMPUTE_INVERSE_3x3_SYMMETRIC_VAR(mI,tI,\
    MLA_CONC(m,t,0,0),MLA_CONC(m,t,0,1),MLA_CONC(m,t,0,2),\
    MLA_CONC(m,t,0,1),MLA_CONC(m,t,1,1),MLA_CONC(m,t,1,2),\
    MLA_CONC(m,t,0,2),MLA_CONC(m,t,1,2),MLA_CONC(m,t,2,2))

/* inverse of a 3x3 general matrix, it populates mI variable
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2.
// mI: matrix inverse name stem, tI: matrix inverse index.
// ex: mI = gI and tI = U => MLA_CONC(mI,tI,0,0) = gI_U0U0. */
#define Matrix_Inverse_3x3_General_Var(m,t,mI,tI) \
   MLA_COMPUTE_INVERSE_3x3_GENERAL_VAR(mI,tI,\
    MLA_CONC(m,t,0,0),MLA_CONC(m,t,0,1),MLA_CONC(m,t,0,2),\
    MLA_CONC(m,t,1,0),MLA_CONC(m,t,1,1),MLA_CONC(m,t,1,2),\
    MLA_CONC(m,t,2,0),MLA_CONC(m,t,2,1),MLA_CONC(m,t,2,2))


/* determinant of a 3x3 symmetric matrix for field.
// ijk: the point on the patch
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2. */
#define Matrix_Determinant_3x3_Symmetric_Field(m,t,ijk) \
 (Matrix_Determinant_3x3(\
   MLA_CONC(m,t,0,0)[ijk],MLA_CONC(m,t,0,1)[ijk],MLA_CONC(m,t,0,2)[ijk],\
   MLA_CONC(m,t,0,1)[ijk],MLA_CONC(m,t,1,1)[ijk],MLA_CONC(m,t,1,2)[ijk],\
   MLA_CONC(m,t,0,2)[ijk],MLA_CONC(m,t,1,2)[ijk],MLA_CONC(m,t,2,2)[ijk]))

/* determinant of a 3x3 symmetric matrix for variable.
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2. */
#define Matrix_Determinant_3x3_Symmetric_Var(m,t) \
 (Matrix_Determinant_3x3(\
   MLA_CONC(m,t,0,0),MLA_CONC(m,t,0,1),MLA_CONC(m,t,0,2),\
   MLA_CONC(m,t,0,1),MLA_CONC(m,t,1,1),MLA_CONC(m,t,1,2),\
   MLA_CONC(m,t,0,2),MLA_CONC(m,t,1,2),MLA_CONC(m,t,2,2)))

/* determinant of a 3x3 general matrix for field.
// ijk: the point on the patch
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2. */
#define Matrix_Determinant_3x3_General_Field(m,t,ijk) \
 (Matrix_Determinant_3x3(\
   MLA_CONC(m,t,0,0)[ijk],MLA_CONC(m,t,0,1)[ijk],MLA_CONC(m,t,0,2)[ijk],\
   MLA_CONC(m,t,1,0)[ijk],MLA_CONC(m,t,1,1)[ijk],MLA_CONC(m,t,1,2)[ijk],\
   MLA_CONC(m,t,2,0)[ijk],MLA_CONC(m,t,2,1)[ijk],MLA_CONC(m,t,2,2)[ijk]))

/* determinant of a 3x3 general matrix for variable.
// m: matrix name stem, t: matrix index.
// ex: m = g and t = D, => MLA_CONC(m,t,1,2) = g_D1D2. */
#define Matrix_Determinant_3x3_General_Var(m,t) \
 (Matrix_Determinant_3x3(\
   MLA_CONC(m,t,0,0),MLA_CONC(m,t,0,1),MLA_CONC(m,t,0,2),\
   MLA_CONC(m,t,1,0),MLA_CONC(m,t,1,1),MLA_CONC(m,t,1,2),\
   MLA_CONC(m,t,2,0),MLA_CONC(m,t,2,1),MLA_CONC(m,t,2,2)))

/* determinant of a general 3x3 matrix with component a?? */
#define Matrix_Determinant_3x3(a00,a01,a02,a10,a11,a12,a20,a21,a22) \
    (a00 *a11 *a22  - a00 *a12 *a21  -\
     a01 *a10 *a22  + a01 *a12 *a20  +\
     a02 *a10 *a21  - a02 *a11 *a20)


/* matrix storage format */
typedef enum MATRIX_SF_T
{
  UNDEF_SF = 0/* undefined */,
  REG_SF/* regular storage format */,
  RMO_SF/* row major order format */,
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
  Uint reg_f: 1;/* regular format */
  Uint rmo_f: 1;/* row major order format */
  Uint tri_f: 1;/* 1 if tripet storage format, 0 otherwise*/
  Uint ccs_f: 1;/* 1 if compressed column storage format, 0 otherwise */
  Uint crs_f: 1;/* 1 if compressed row storage format,0 otherwise */
  Uint tri_l_f: 1;/* 1 if long tripet storage format, 0 otherwise*/
  Uint ccs_l_f: 1;/* 1 if long compressed column storage format, 0 otherwise */
  Uint crs_l_f: 1;/* 1 if long compressed row storage format,0 otherwise */
  
  long row;/* number of rows */
  long col;/* number of columns */
  
  struct/* triplet storage format */
  {
    int *row;
    int *col;
    double *a;
  }tri[1];
  
  struct/* row major order format */
  {
    double *A;
  }rmo[1];
  
  struct/* compressed column storage format */
  {
    int *Ap;
    int *Ai;
    double *Ax;
    /* coarse-graining for matrix read optimization */
    int *Ap_cg;
    int *i_cg;
    /* slice each Ap[c+1]-Ap[c] into equal pieces with known 
    // value of index for Ai intervals thus for given row r we know in 
    // which interval r can be found as apposed to just look for r blindly.
    // it's kind of coarse-graining over the value Ai's.
    // the logic is:
    // i = Ap_cg[c] => 0    <= Ai[i_cg[i]]   < row0
    // i + 1        => row0 <= Ai[i_cg[i+1]] < row1
    // .
    // .
    // .
    // till (i < Ap_cg[c+1]) and naturally i < Ap[c+1].
    // thus if given row r is between row0 and row1, 
    // it looks in second piece and not other pieces.
    // NOTE: the idea is pretty much like the CCS format itself. */
    int Nslice;/* number of slices(pieces) */
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
Matrix_T *cast_matrix_rmo(Matrix_T *const m);
void copy_reg2reg(const Matrix_T *const reg1,Matrix_T *const reg2);
void copy_ccs2ccs(const Matrix_T *const ccs1,Matrix_T *const ccs2);
void matrix_tests(void);
void precondition(Matrix_T *const a,double *const b);
Matrix_T *invert_matrix(Matrix_T *const M);
int matrix_by_vector(const Matrix_T *const m, const double *const v,double *const b,const Flag_T flag);
Matrix_T *matrix_by_matrix(const Matrix_T *const a, const Matrix_T *const b,Matrix_T *const result,const char *const dir);
Matrix_T *cast_matrix_ccs_long(Matrix_T *const m);
void copy_ccs_long2ccs_long(const Matrix_T *const ccs_l1,Matrix_T *const ccs_l2);
Matrix_T *compress_stack2ccs(Matrix_T **const S,const Uint nm,const Uint *const nr,const Uint Nrow,const Uint Ncol,const Flag_T flg);
Matrix_T *CCSOpCCS(Matrix_T *const ccs2,Matrix_T *const ccs1,const char Op);
int matrix_comparison(const Matrix_T *const M1,const Matrix_T *const M2);
Matrix_T *alloc_matrix(const Matrix_SF_T type_e,const long row,const long col);
void free_matrix(Matrix_T *m);



#endif





