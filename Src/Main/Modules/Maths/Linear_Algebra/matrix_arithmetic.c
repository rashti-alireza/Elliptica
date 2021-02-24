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
    Error0("Null Matrix \"m\".");
  else if (!v)
    Error0("Null Vector \"v\".");
  else if(!b)
    Error0("No memory for Resultant vector \"b\".");
    
  const long row = m->row;
  const long col = m->col;
  long r,c;
  
  /* if empty */
  if (!row || !col)
    Error0("The matrix \"m\" is empty.");

  /* if the entries of the resultant vector 
    has not been initialized to zero thus initialization is needed */
  if (flag == INITIALIZE)
  {
    for (r = 0; r < row; ++r)
      b[r] = 0;
  }
  else if (flag != NOT_INITIALIZE)
  {
    Error0("No such flag has been defined for this function.");
  }
  
  if (m->reg_f)
  {
    double **const M = m->reg->A;
    
    for (r = 0; r < row; ++r)
      for (c = 0; c < col; ++c)
        b[r] += M[r][c]*v[c];
  }
  else if (m->rmo_f)
  {
    double *const M = m->rmo->A;
    
    for (r = 0; r < row; ++r)
      for (c = 0; c < col; ++c)
        b[r] += M[i_j_to_ij(col,r,c)]*v[c];
  }
  else if (m->tri_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->ccs_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->crs_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->tri_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->ccs_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->crs_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else
    Error0("No matrix format is defined for this given matrix.");

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
    Error0("Null Matrix \"a\".");
  if (!b)
    Error0("Null Matrix \"b\".");
  
  Matrix_T *d = 0;
  const long a_row = a->row;
  const long a_col = a->col;
  const long b_row = b->row;
  const long b_col = b->col;
  long r,i,c;
  
  /* more checks */
  if (!a_row || !a_col)
    Error0("The matrix \"a\" is empty.");
  if (!b_row || !b_col)
    Error0("The matrix \"b\" is empty.");
  
  if (strcmp_i("a*b",dir))
  {
    if (a_col != b_row)
      Error0("The dimensions of matrix a and b are not matched.");
      
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
      Error0(INCOMPLETE_FUNC);
  }/* end of if (strcmp_i("AxB",dir)) */
  else if (strcmp_i("a*Transpose(b)",dir))
  {
    if (a_col != b_col)
      Error0("The dimensions of matrix a and Transpose(b) are not matched.\n");
      
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
    else if (a->rmo_f && b->rmo_f)
    {
      double *const A = a->rmo->A;
      double *const B = b->rmo->A;
      double **const D = d->reg->A;
      
      for (r = 0; r < a_row; ++r)
        for (c = 0; c < b_row; ++c)
          for (i = 0; i < b_col; ++i)
              D[r][c] += A[i_j_to_ij(b_col,r,i)]*B[i_j_to_ij(b_col,c,i)];
    }
    else
      Error0(INCOMPLETE_FUNC);
  }
  else
    Error0(INCOMPLETE_FUNC);
  
  return d;
}

/* add or subtract to CCS matrices as ccs2-ccs1 or ccs2+ccs1 
// and cast the resultant to CCS too and return it.
// note the char (+,-) determines the operator.
// note if the rows or columns are inconsistent it gives error.
// ->return value: resultant in CCS format */
Matrix_T *CCSOpCCS(Matrix_T *const ccs2,Matrix_T *const ccs1,const char Op)
{
  const long Nr = ccs2->row;
  const long Nc = ccs2->col;
  Matrix_T *res = 0;
  int *Ap   = 0;
  int *Ai   = 0;
  double *Ax = 0;
  double ax;
  long tNN0 = 0;/* total number of none zero entries */
  long NN0;/* number of none zero entries in each column */
  const double DropLimit = 0; 
  long r,c;/* row and column */
  
  if (Nr != ccs1->row)
    Error0("Rows of the matrices are not matched.\n");
  if (Nc != ccs1->col)
    Error0("columns of the matrices are not matched.\n");
  
  Ap = calloc((long Uint)Nc+1,sizeof(*Ap));
  IsNull(Ap);
  res = alloc_matrix(CCS_SF,Nr,Nc);
  
  if (Op == '+')
  {
    for (c = 0; c < Nc; ++c)
    {
      NN0 = 0;
      for (r = 0; r < Nr; ++r)
      {
        ax = read_matrix_entry_ccs(ccs2,r,c)+read_matrix_entry_ccs(ccs1,r,c);
        if (GRT(ABSd(ax),DropLimit))
        {
          Ai = realloc(Ai,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ai));
          IsNull(Ai);
          Ax = realloc(Ax,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ax));
          IsNull(Ax);
          Ai[Ap[c]+NN0] = (int)r;
          Ax[Ap[c]+NN0] = ax;
          NN0++;
          tNN0++;
        }
      }
      Ap[c+1] = (int)tNN0;
    }
  }
  else if (Op == '-')
  {
    for (c = 0; c < Nc; ++c)
    {
      NN0 = 0;
      for (r = 0; r < Nr; ++r)
      {
        
        ax = read_matrix_entry_ccs(ccs2,r,c)-read_matrix_entry_ccs(ccs1,r,c);
        if (GRT(ABSd(ax),DropLimit))
        {
          Ai = realloc(Ai,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ai));
          IsNull(Ai);
          Ax = realloc(Ax,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ax));
          IsNull(Ax);
          Ai[Ap[c]+NN0] = (int)r;
          Ax[Ap[c]+NN0] = ax;
          NN0++;
          tNN0++;
        }
      }
      Ap[c+1] = (int)tNN0;
    }
  }
  else
    Error0(NO_OPTION);

  res->ccs->Ap = Ap;
  res->ccs->Ai = Ai;
  res->ccs->Ax = Ax;  
  
  return res;
}

/* compare if M1 == M2.
// ->return value: 1 if they are equal, 0 if not. */
int matrix_comparison(const Matrix_T *const M1,const Matrix_T *const M2)
{
  const long Nr = M1->row;
  const long Nc = M1->col;
  long r,c;
  
  if (M1->row != M2->row)
    return 0;
  if (M1->col != M2->col)
    return 0;
  
  if (M1->reg_f && M2->reg_f)
  {
    double **const m1 = M1->reg->A;
    double **const m2 = M2->reg->A;
    
    for (r = 0; r < Nr; ++r)
    {
      for (c = 0; c < Nc; ++c)
      {
        if (!EQL(m1[r][c],m2[r][c]))
          return 0;
      }
    }
  }/* end of if (M1->reg_f && M2->reg_f) */
  else
    Error0(NO_OPTION);
  
  return 1;
}



