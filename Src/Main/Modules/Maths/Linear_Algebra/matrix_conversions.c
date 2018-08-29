/*
// Alireza Rashti
// August 2018
*/

#include "matrix_conversions.h"

/* casting given matrix m to Compressed Column Storage Format.
// it keeps the given matrix and makes a new matrix with the specified
// format.
// note 1: if the given matrix is ccs itself, it returns
// a copy of the given matrix and returns it.
// note 2: if the given matrix is 0 by 0, it returns null.
// ->return value: a matirx with ccs format.
*/
Matrix_T *cast_matrix_ccs(Matrix_T *const m)
{
  Matrix_T *ccs = 0;
  
  /* if empty, return null */
  if (!m->row || !m->col)
    return ccs;
  
  ccs = alloc_matrix(CCS_SF,m->row,m->col);
  
  if (m->reg_f)
  {
    double DropLimit = 0;
    if (get_parameter("Ignore_Number_Less_Than_in_CCS_format"))	
      DropLimit = GetParameterD("Ignore_Number_Less_Than_in_CCS_format");
      
    convert_reg2ccs(m,ccs,DropLimit);
  }
  else if (m->tri_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->ccs_f)
  {
    copy_ccs2ccs(m,ccs);
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
  
  return ccs;
}

/* casting given matrix m to Regular Format.
// it keeps the given matrix and makes a new matrix with the specified
// format.
// note 1: if the given matrix is regular itself, it returns
// a copy of the given matrix and returns it.
// note 2: if the given matrix is 0 by 0, it returns null.
// ->return value: a matirx with reg format.
*/
Matrix_T *cast_matrix_reg(Matrix_T *const m)
{
  Matrix_T *reg = 0;
  
  /* if empty, return null */
  if (!m->row || !m->col)
    return reg;
  
  reg = alloc_matrix(REG_SF,m->row,m->col);
  
  if (m->reg_f)
  {
    copy_reg2reg(m,reg);
  }
  else if (m->tri_f)
  {
    abortEr(INCOMPLETE_FUNC);
  }
  else if (m->ccs_f)
  {
    convert_ccs2reg(m,reg);
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
  
  return reg;
}

/* converting a ccs format to reg format. */
static void convert_ccs2reg(const Matrix_T *const ccs,Matrix_T *const reg)
{
  const long Nc = ccs->col;
  double **const m = reg->reg->A;
  const int *const Ap   = ccs->ccs->Ap;
  const int *const Ai   = ccs->ccs->Ai;
  const double *const Ax = ccs->ccs->Ax;
  long c,i;
  
  for (c = 0; c < Nc; ++c)
    for (i = Ap[c]; i < Ap[c+1]; ++i)
      m[Ai[i]][c] = Ax[i];
  
}

/* copying a regular matrix like reg2 = reg1.
// note 1: reg2 must already have been allocated memory.
// note 2: if either of reg1->col or reg1->row is zero it won't copy anything.
*/
void copy_reg2reg(const Matrix_T *const reg1,Matrix_T *const reg2)
{
  if (reg1->row && reg1->col)
    return;
    
  const long Nr = reg1->row;
  const long Nc = reg1->col;
  double **const m1 = reg1->reg->A;
  double **const m2 = reg2->reg->A;
  long r,c;
  
  for (r = 0 ; r < Nr; ++r)
    for (c = 0; c < Nc; ++c)
      m2[r][c] = m1[r][c];
}


/* copying a ccs matrix like ccs2 = ccs1.
// note 1: ccs2 must already have been allocated memory.
// note 2: if cc1->row or cc1->col is zero, it won't copy anything.
*/
void copy_ccs2ccs(const Matrix_T *const ccs1,Matrix_T *const ccs2)
{
  if (ccs1->col && ccs1->row)
    return;
    
  const long Nc = ccs1->col;
  const int *const Ap1   = ccs1->ccs->Ap;/* note: 
                                        // Ap1[Nc] = total number of 
                                        // none zero entries.
                                        */
  const int *const Ai1   = ccs1->ccs->Ai;
  const double *const Ax1 = ccs1->ccs->Ax;
  int *const Ap2   = calloc((long unsigned)Nc+1,sizeof(*Ap2));
  int *const Ai2   = calloc((long unsigned)Ap1[Nc],sizeof(*Ai2));
  double *const Ax2 = calloc((long unsigned)Ap1[Nc],sizeof(*Ax2));
  long c,i;
  
  for (c = 0; c < Nc; c++)
  {
    Ap2[c] = Ap1[c];
    for (i = Ap1[c]; i < Ap1[c+1]; ++i)
    {
      Ai2[i] = Ai1[i];
      Ax2[i] = Ax1[i];
    }
  }
  Ap2[Nc] = Ap1[Nc];
  
  ccs2->ccs->Ap = Ap2;
  ccs2->ccs->Ai = Ai2;
  ccs2->ccs->Ax = Ax2;
}

/* converting a regular storage format to ccs format.
// note: the matrix will be made in CCS such that be in order in row.
*/
static void convert_reg2ccs(const Matrix_T *const reg,Matrix_T *const ccs,const double DropLimit)
{
  const long Nr = reg->row;
  const long Nc = reg->col;
  double **const m = reg->reg->A;
  int *Ap   = calloc((long unsigned)Nc+1,sizeof(*Ap));
  int *Ai   = 0;
  double *Ax = 0;
  long tNN0 = 0;/* total number of none zero entries */
  long NN0;/* number of none zero entries in each column */
  long r,c;/* row and column */
  
  for (c = 0; c < Nc; ++c)
  {
    NN0 = 0;
    for (r = 0; r < Nr; ++r)
    {
      if (GRT(ABS(m[r][c]),DropLimit))
      {
        Ai = realloc(Ai,(long unsigned)(Ap[c]+NN0+1)*sizeof(*Ai));
        pointerEr(Ai);
        Ax = realloc(Ax,(long unsigned)(Ap[c]+NN0+1)*sizeof(*Ax));
        pointerEr(Ax);
        Ai[Ap[c]+NN0] = (int)r;
        Ax[Ap[c]+NN0] = m[r][c];
        NN0++;
        tNN0++;
      }
    }
    Ap[c+1] = (int)tNN0;
  }
  
  ccs->ccs->Ap = Ap;
  ccs->ccs->Ai = Ai;
  ccs->ccs->Ax = Ax;
}
