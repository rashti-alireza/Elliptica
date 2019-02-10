/*
// Alireza Rashti
// Feb 2019
*/

#include "matrix_arithmetic.h"

/* get matirx m and multiply it by vector v and put the result to b */
int matrix_by_vector(const Matrix_T *const m, const double *const v,double *const b)
{
  if (!m)
    abortEr("Null Matrix \"m\".\n");
  else if (!v)
    abortEr("Null Vector \"v\".\n");
  else if(!b)
    abortEr("No memory for Resultant vector \"b\".\n");
    
  const long row = m->row;
  const long col = m->col;
  long r,c;
  
  /* if empty */
  if (!row || !col)
    abortEr("The matrix \"m\" is empty.\n");
  
  if (m->reg_f)
  {
    double **M = m->reg->A;
    
    for (r = 0; r < row; ++r)
      for (c = 0; c < col; ++c)
        b[r] = M[r][c]*v[c];
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
    abortEr("No matrix format is defined for this given matrix.\n");

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
    abortEr("Null Matrix \"a\".\n");
  if (!b)
    abortEr("Null Matrix \"b\".\n");
  
  Matrix_T *d = 0;
  const long a_row = a->row;
  const long a_col = a->col;
  const long b_row = b->row;
  const long b_col = b->col;
  long r,i,c;
  
  /* more checks */
  if (!a_row || !a_col)
    abortEr("The matrix \"a\" is empty.\n");
  if (!b_row || !b_col)
    abortEr("The matrix \"b\" is empty.\n");
  
  if (strcmp_i("a*b",dir))
  {
    if (a_col != b_row)
      abortEr("The dimensions of matrix a and b are not matched.\n");
      
    d = alloc_matrix(REG_SF,a_row,b_col);
    if (a->reg_f && b->reg_f)
    {
      double **A = a->reg->A;
      double **B = b->reg->A;
      double **D = d->reg->A;
      
      for (r = 0; r < a_row; ++r)
        for (c = 0; c < b_col; ++c)
          for (i = 0; i < b_row; ++i)
              D[r][c] = A[r][i]*B[i][c];
    }
    else if (a->ccs_f && b->ccs_f)
    {
      int *Ap_a    = a->ccs->Ap;
      int *Ai_a    = a->ccs->Ai;
      double *Ax_a = a->ccs->Ax;
      int *Ap_b    = b->ccs->Ap;
      int *Ai_b    = b->ccs->Ai;
      double *Ax_b = b->ccs->Ax;
      double **D = d->reg->A;
      
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
      double **A = a->reg->A;
      double **B = b->reg->A;
      double **D = d->reg->A;
      
      for (r = 0; r < a_row; ++r)
        for (c = 0; c < b_row; ++c)
          for (i = 0; i < b_col; ++i)
              D[r][c] = A[r][i]*B[c][i];
    }
    else
      abortEr(INCOMPLETE_FUNC);
  }
  else
    abortEr(INCOMPLETE_FUNC);
  
  return d;
}
