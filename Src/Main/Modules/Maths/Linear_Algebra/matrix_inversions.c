/*
// Alireza Rashti
// September 2018
*/

#include "matrix_inversions.h"

/* finding the invert of the given matrix and returning it.
// algorithm:
// ==========
//
// setting up the following equation:
// M.x = ej j = 0 ... row
// M is an convertible row by col matrix.
// ej a row vector that only its jth row is 1 and the rest is 0.
// x is the jth column of matrix Invert(M).
//
// ->return value: invert(M) such that invert(M).M = I
*/
Matrix_T *invert_matrix(Matrix_T *const M)
{
  Matrix_T *invertM = alloc_matrix(REG_SF,M->row,M->col);
  Matrix_T *ccsM = 0;
  const long row  = M->row;
  const long col  = M->col;
  const long unsigned max = col>=row?(long unsigned)col:(long unsigned)row;
  int *Ap; 
  int *Ai;
  double *Ax;
  double **const inv = invertM->reg->A;
  double *ej = alloc_double((unsigned)M->row);
  double *x = alloc_double((unsigned)M->row);
  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  void *Symbolic,*Numeric;
  int status;
  long r,c;
  
  /* NOTE: this function only works for 
  // dimensions in "int" range. one reason is this function uses
  // umfpack_di_*. so make sure this condition is held.
  */
  assert(max<(long unsigned)INT_MAX);
  
  if (M->ccs_f)
    ccsM = M;
  else
    ccsM = cast_matrix_ccs(M);
  
  Ap = ccsM->ccs->Ap;
  Ai = ccsM->ccs->Ai;
  Ax = ccsM->ccs->Ax;
  umfpack_di_defaults(Control);
  
  status = umfpack_di_symbolic((int)row,(int)col,Ap,Ai,Ax,&Symbolic,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);
     
  status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);
     
  umfpack_di_free_symbolic(&Symbolic);
  
  for (c = 0; c < col; ++c)
  {
    ej[c] = 1;
    
    status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x,ej,Numeric,Control,Info);
    if(status != UMFPACK_OK)
      umfpack_error_di(Control,status,__FILE__,__LINE__);
      
    for (r = 0; r < row; ++r)
      inv[r][c] = x[r];
    
    ej[c] = 0;
  }
  
  /*freeing*/
  umfpack_di_free_numeric(&Numeric);
  if (!M->ccs_f)  free_matrix(ccsM);
  free(ej);
  free(x);
  
  /* test if this algorithm works perfectly */
  if (0)
    test_invert_matrix(M,invertM);
    
  return invertM;
}

/* testing if invertM.M = I */
static void test_invert_matrix(const Matrix_T *const invertM,const Matrix_T *const M)
{
  const long row = M->row;
  const long col = M->col;
  double **const inv = invertM->reg->A;
  double **const m    = M->reg->A;
  long i,j,k;
  
  printf("Testing Inverting matrix algorithm by checking "
          "invertM.M = I:\n");

  for (i = 0; i < row; ++i)
  {
    double Iij;
    
    for (j = 0; j < col; ++j)
    {
      Iij = 0;
      for (k = 0; k < col; ++k)
        Iij += inv[i][k]*m[k][j];
      
      if (i == j)
      {
        if(!EQL(Iij,1))
          abortEr("Inverting matrix algorithm doesn't work!\n");
      }
      else
      {
        if(!EQL(Iij,0))
          abortEr("Inverting matrix algorithm doesn't work!\n");
      }
      
    }/* end of for (j = 0; j < col; ++j) */
  }/* end of for (i = 0; i < row; ++i) */

  printf("Inverting matrix algorithm works perfectly.\n");
}
