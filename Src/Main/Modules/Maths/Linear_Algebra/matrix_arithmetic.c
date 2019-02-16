/*
// Alireza Rashti
// Feb 2019
*/

#include "matrix_arithmetic.h"

/* get matirx m and multiply it by vector v and put the result to b.
// note: entries of b must be zero.
// flag options:
// o. INITIALIZE # if one wants to set the entries of b vector to zero
// o. NOT_INITIALIZE # if entries of b is already zero thus there is no
// 		     # need to initialize again.
//
 */
int matrix_by_vector(const Matrix_T *const m, const double *const v,double *const b,const Flag_T flag)
{
  if (!m)
    abortEr("Null Matrix \"m\".");
  else if (!v)
    abortEr("Null Vector \"v\".");
  else if(!b)
    abortEr("No memory for Resultant vector \"b\".");
    
  const long row = m->row;
  const long col = m->col;
  long r,c;
  
  /* if empty */
  if (!row || !col)
    abortEr("The matrix \"m\" is empty.");

  /* if the entries of the resultant vector 
    has not been initialized to zero thus initialization is needed */
  if (flag == INITIALIZE)
  {
    for (r = 0; r < row; ++r)
      b[r] = 0;
  }
  else if (flag != NOT_INITIALIZE)
  {
    abortEr("No such flag has been defined for this function.");
  }
  
  if (m->reg_f)
  {
    double **const M = m->reg->A;
    
    for (r = 0; r < row; ++r)
      for (c = 0; c < col; ++c)
        b[r] += M[r][c]*v[c];
  }
  else if (m->tri_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->ccs_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->crs_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->tri_l_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->ccs_l_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->crs_l_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else
    abortEr("No matrix format is defined for this given matrix.");

  return 0;
}

/* get matrix a and multiply it by matrix b according to 
// the directivce dir and allocate memory put the result to d.
//
// directives:
// 	1. a*b            # multiply a matrix by b matrix
// 	2. a*Transpose(b) # multiply a matrix by Transpose(b)
//
// ->return value: a.b = d. */
Matrix_T *matrix_by_matrix(const Matrix_T *const a, const Matrix_T *const b,const char *const dir)
{
  /* some checks */
  if (!a)
    abortEr("Null Matrix \"a\".");
  if (!b)
    abortEr("Null Matrix \"b\".");
  
  Matrix_T *d = 0;
  const long a_row = a->row;
  const long a_col = a->col;
  const long b_row = b->row;
  const long b_col = b->col;
  long r,i,c;
  
  /* more checks */
  if (!a_row || !a_col)
    abortEr("The matrix \"a\" is empty.");
  if (!b_row || !b_col)
    abortEr("The matrix \"b\" is empty.");
  
  if (strcmp_i("a*b",dir))
  {
    if (a_col != b_row)
      abortEr("The dimensions of matrix a and b are not matched.");
      
    d = alloc_matrix(REG_SF,a_row,b_col);
    if (a->reg_f && b->reg_f)
    {
      double **const A = a->reg->A;
      double **const B = b->reg->A;
      double **const D = d->reg->A;
      
      for (r = 0; r < a_row; ++r)
        for (c = 0; c < b_col; ++c)
          for (i = 0; i < b_row; ++i)
              D[r][c] += A[r][i]*B[i][c];
    }
    else if (a->ccs_f && b->ccs_f)
    {
      const int *const Ap_a    = a->ccs->Ap;
      const int *const Ai_a    = a->ccs->Ai;
      const double *const Ax_a = a->ccs->Ax;
      const int *const Ap_b    = b->ccs->Ap;
      const int *const Ai_b    = b->ccs->Ai;
      double *const Ax_b       = b->ccs->Ax;
      double **const D         = d->reg->A;
      
      for (c = 0; c < b_col; ++c)
        for (r = Ap_b[c]; r < Ap_b[c+1]; ++r)
          for (i = Ap_a[Ai_b[r]]; i < Ap_a[Ai_b[r]+1]; ++i)
            D[Ai_a[i]][c] += Ax_a[i]*Ax_b[r];
    }
    else
      abortEr(INCOMPLETE_FUNC);
  }/* end of if (strcmp_i("AxB",dir)) */
  else if (strcmp_i("a*Transpose(b)",dir))
  {
    if (a_col != b_col)
      abortEr("The dimensions of matrix a and Transpose(b) are not matched.\n");
      
    d = alloc_matrix(REG_SF,a_row,b_row);
    if (a->reg_f && b->reg_f)
    {
      double **const A = a->reg->A;
      double **const B = b->reg->A;
      double **const D = d->reg->A;
      
      for (r = 0; r < a_row; ++r)
        for (c = 0; c < b_row; ++c)
          for (i = 0; i < b_col; ++i)
              D[r][c] += A[r][i]*B[c][i];
    }
    else
      abortEr(INCOMPLETE_FUNC);
  }
  else
    abortEr(INCOMPLETE_FUNC);
  
  return d;
}
