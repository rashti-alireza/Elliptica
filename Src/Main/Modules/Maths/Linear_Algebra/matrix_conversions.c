/*
// Alireza Rashti
// August 2018
*/

#include "matrix_conversions.h"


/* casting a stack of matrices S to Compressed Column Storage Format.
// the stack is collection of matrices such that each of them has
// the same column but maybe different number of rows.
// note: all of the members of the stack must have same number of columns.
// note: some of the matrices in the stack might have all zero entries
// in which case one might determine this matrix with null pointer.
// note: if the flag is YES it will free all of the matrices of the stack
// it won't free the stack itself.
// ->return value: the compressed stack matrix in ccs format. */
Matrix_T *compress_stack2ccs
( Matrix_T **const S/* stack of matrices */,
  const Uint nm/* number of matrices in the stack */,
  const Uint *const nr/* number of rows in each matrix */,
  const Uint Nrow/* total number of rows */,
  const Uint Ncol/* common number of columns */,
  const Flag_T flg/* if YES means free all of the matrices from the stack, do nothing otherwise */
)
{
  Matrix_T *ccs = alloc_matrix(CCS_SF,Nrow,Ncol);
  int *Ap   = calloc((long Uint)Ncol+1,sizeof(*Ap));
  int *Ai   = 0;
  double *Ax = 0;
  
  long tNN0 = 0;/* total number of none zero entries */
  long NN0;/* number of none zero entries in each column */
  long r,c;/* row and column */
  long R = 0;
  const double DropLimit = 0;
  double **m;
  Uint i;
  
  assert(Nrow < INT_MAX);
  
  for (c = 0; c < Ncol; ++c)
  {
    NN0 = 0;
    R = 0;
    /* go thru all of the rows */
    for (i = 0; i < nm; ++i)
    {
      if (S[i])
      {
        m = S[i]->reg->A;
        for (r = 0; r < nr[i]; ++r)
        {
          if (GRT(ABSd(m[r][c]),DropLimit))
          {
            Ai = realloc(Ai,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ai));
            IsNull(Ai);
            Ax = realloc(Ax,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ax));
            IsNull(Ax);
            Ai[Ap[c]+NN0] = (int)(r+R);
            Ax[Ap[c]+NN0] = m[r][c];
            NN0++;
            tNN0++;
          }
        }
      }/* end of if(C) */
      R += nr[i];
    }/* end of for (p = 0; p < np; ++p) */
    Ap[c+1] = (int)tNN0;
  }/* end of for (c = 0; c < Nc; ++c) */
  
  /* free S matrices */
  if (flg == YES)
  {
    for (i = 0; i < nm; ++i)
      if (S[i])
        free_matrix(S[i]);
  }
      
  ccs->ccs->Ap = Ap;
  ccs->ccs->Ai = Ai;
  ccs->ccs->Ax = Ax;
  
  return ccs;
}

/* casting given matrix m to Compressed Column Storage Format.
// it keeps the given matrix and makes a new matrix with the specified
// format.
// note 1: if the given matrix is ccs itself, it returns
// a copy of the given matrix and returns it.
// note 2: if the given matrix is 0 by 0, it returns null.
// ->return value: the matirx in ccs format.
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
    if (get_parameter("matrix_ccs_drop_below"))	
      DropLimit = PgetdEZ("matrix_ccs_drop_below");
      
    convert_reg2ccs(m,ccs,DropLimit);
  }
  else if (m->tri_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->ccs_f)
  {
    copy_ccs2ccs(m,ccs);
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
    Error0("No matrix format is defined for this given matrix.\n");
  
  return ccs;
}

/* casting given matrix m to row major order Storage Format.
// it keeps the given matrix and makes a new matrix with the specified
// format.
// note 1: if the given matrix is 0 by 0, it returns null.
// ->return value: the matirx in rmo format. */
Matrix_T *cast_matrix_rmo(Matrix_T *const m)
{
  Matrix_T *cast = 0;
  
  /* if empty, return null */
  if (!m->row || !m->col)
    return 0;
  
  if (m->reg_f)
  {
    cast                 = alloc_matrix(RMO_SF,m->row,m->col);
    double *const rmo_a  = cast->rmo->A;
    double **const reg_a = m->reg->A;
    const long Nr        = m->row;
    const long Nc        = m->col;
    long r,c;
  
    for (r = 0; r < Nr; ++r)
      for (c = 0; c < Nc; ++c)
        rmo_a[i_j_to_ij(Nc,r,c)] = reg_a[r][c];
  }
  else
    Error0(NO_OPTION);
  
  return cast;
}

/* casting given matrix m to Long Compressed Column Storage Format.
// it keeps the given matrix and makes a new matrix with the specified
// format.
// note 1: if the given matrix is ccs itself, it returns
// a copy of the given matrix and returns it.
// note 2: if the given matrix is 0 by 0, it returns null.
// ->return value: a matirx with ccs long format.
*/
Matrix_T *cast_matrix_ccs_long(Matrix_T *const m)
{
  Matrix_T *ccs_l = 0;
  
  /* if empty, return null */
  if (!m->row || !m->col)
    return ccs_l;
  
  ccs_l = alloc_matrix(CCS_L_SF,m->row,m->col);
  
  if (m->reg_f)
  {
    double DropLimit = 0;
    if (get_parameter("matrix_ccs_drop_below"))	
      DropLimit = PgetdEZ("matrix_ccs_drop_below");
      
    convert_reg2ccs_long(m,ccs_l,DropLimit);
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
    copy_ccs_long2ccs_long(m,ccs_l);
  }
  else if (m->crs_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else
    Error0("No matrix format is defined for this given matrix.\n");
  
  return ccs_l;
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
    Error0(INCOMPLETE_FUNC);
  }
  else if (m->ccs_f)
  {
    convert_ccs2reg(m,reg);
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
    convert_ccs_long2reg(m,reg);
  }
  else if (m->crs_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else
    Error0("No matrix format is defined for this given matrix.\n");
  
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
  
  reg->row = ccs->row;
  reg->col = ccs->col;
  
  for (c = 0; c < Nc; ++c)
    for (i = Ap[c]; i < Ap[c+1]; ++i)
      m[Ai[i]][c] = Ax[i];
  
}

/* converting a ccs_long format to reg format. */
static void convert_ccs_long2reg(const Matrix_T *const ccs,Matrix_T *const reg)
{
  const long Nc = ccs->col;
  double **const m = reg->reg->A;
  const long *const Ap   = ccs->ccs_long->Ap;
  const long *const Ai   = ccs->ccs_long->Ai;
  const double *const Ax = ccs->ccs_long->Ax;
  long c,i;
  
  reg->row = ccs->row;
  reg->col = ccs->col;
  
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
  if (!ccs1->col || !ccs1->row)
    return;
    
  const long Nc = ccs1->col;
  const int *const Ap1   = ccs1->ccs->Ap;/* note: 
                                        // Ap1[Nc] = total number of 
                                        // none zero entries.
                                        */
  const int *const Ai1   = ccs1->ccs->Ai;
  const double *const Ax1 = ccs1->ccs->Ax;
  int *const Ap2   = calloc((long Uint)Nc+1,sizeof(*Ap2));
  int *const Ai2   = calloc((long Uint)Ap1[Nc],sizeof(*Ai2));
  double *const Ax2 = calloc((long Uint)Ap1[Nc],sizeof(*Ax2));
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
  
  ccs2->row = ccs1->row;
  ccs2->col = ccs1->col;
  ccs2->ccs->Ap = Ap2;
  ccs2->ccs->Ai = Ai2;
  ccs2->ccs->Ax = Ax2;
}

/* copying a ccs long matrix like ccs_l2 = ccs_l1.
// note 1: ccs_l2 must already have been allocated memory.
// note 2: if ccs_l1->row or ccs_l1->col is zero, it won't copy anything.
*/
void copy_ccs_long2ccs_long(const Matrix_T *const ccs_l1,Matrix_T *const ccs_l2)
{
  if (!ccs_l1->col || !ccs_l1->row)
    return;
    
  const long Nc = ccs_l1->col;
  const long *const Ap1   = ccs_l1->ccs_long->Ap;/* note: 
                                        // Ap1[Nc] = total number of 
                                        // none zero entries.
                                        */
  const long *const Ai1   = ccs_l1->ccs_long->Ai;
  const double *const Ax1 = ccs_l1->ccs_long->Ax;
  long *const Ap2   = calloc((long Uint)Nc+1,sizeof(*Ap2));
  long *const Ai2   = calloc((long Uint)Ap1[Nc],sizeof(*Ai2));
  double *const Ax2 = calloc((long Uint)Ap1[Nc],sizeof(*Ax2));
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
  
  ccs_l2->row = ccs_l1->row;
  ccs_l2->col = ccs_l1->col;
  ccs_l2->ccs_long->Ap = Ap2;
  ccs_l2->ccs_long->Ai = Ai2;
  ccs_l2->ccs_long->Ax = Ax2;
}


/* converting a regular storage format to ccs format.
// note: the matrix will be made in CCS such that be in order in row.
*/
static void convert_reg2ccs(const Matrix_T *const reg,Matrix_T *const ccs,const double DropLimit)
{
  const long Nr = reg->row;
  const long Nc = reg->col;
  double **const m = reg->reg->A;
  int *Ap   = calloc((long Uint)Nc+1,sizeof(*Ap));
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
      if (GRT(ABSd(m[r][c]),DropLimit))
      {
        Ai = realloc(Ai,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ai));
        IsNull(Ai);
        Ax = realloc(Ax,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ax));
        IsNull(Ax);
        Ai[Ap[c]+NN0] = (int)r;
        Ax[Ap[c]+NN0] = m[r][c];
        NN0++;
        tNN0++;
      }
    }
    Ap[c+1] = (int)tNN0;
  }
  
  ccs->row = reg->row;
  ccs->col = reg->col;
  ccs->ccs->Ap = Ap;
  ccs->ccs->Ai = Ai;
  ccs->ccs->Ax = Ax;
}

/* converting a regular storage format to ccs long format.
// note: the matrix will be made in CCS long such that be in order in row. */
static void convert_reg2ccs_long(const Matrix_T *const reg,Matrix_T *const ccs_l,const double DropLimit)
{
  const long Nr = reg->row;
  const long Nc = reg->col;
  double **const m = reg->reg->A;
  long *Ap   = calloc((long Uint)Nc+1,sizeof(*Ap));
  long *Ai   = 0;
  double *Ax = 0;
  long tNN0 = 0;/* total number of none zero entries */
  long NN0;/* number of none zero entries in each column */
  long r,c;/* row and column */
  
  for (c = 0; c < Nc; ++c)
  {
    NN0 = 0;
    for (r = 0; r < Nr; ++r)
    {
      if (GRT(ABSd(m[r][c]),DropLimit))
      {
        Ai = realloc(Ai,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ai));
        IsNull(Ai);
        Ax = realloc(Ax,(long Uint)(Ap[c]+NN0+1)*sizeof(*Ax));
        IsNull(Ax);
        Ai[Ap[c]+NN0] = r;
        Ax[Ap[c]+NN0] = m[r][c];
        NN0++;
        tNN0++;
      }
    }
    Ap[c+1] = tNN0;
  }
  
  ccs_l->row = reg->row;
  ccs_l->col = reg->col;
  ccs_l->ccs_long->Ap = Ap;
  ccs_l->ccs_long->Ai = Ai;
  ccs_l->ccs_long->Ax = Ax;
}

