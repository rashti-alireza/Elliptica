/*
// Alireza Rashti
// August 2018
*/

#include "matrix_conversions.h"

/* casting given matrix m to Compressed Column Storage Format.
// it keeps the given matrix and makes a new matrix with the specified
// format.
// ->return value: a matirx with ccs format.
*/
Matrix_T *cast_matrix_ccs(Matrix_T *const m)
{

  Matrix_T *ccs = alloc_matrix(CCS_SF,0,0);
  
  if (m->reg_f)
  {
    double DropLimit = 0;
    if (get_parameter("Ignore_Number_Less_Than"))	
      DropLimit = GetParameterD("Ignore_Number_Less_Than");
      
    convert_reg2ccs(m,ccs,DropLimit);
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
  else
    abortEr("No matrix format is defined for this given matrix.\n");
  
  return ccs;
}

/* converting a regular storag format to ccs format.
// note: the matrix will be made in CCS such that be in order in row.
*/
static void convert_reg2ccs(const Matrix_T *const reg,Matrix_T *const ccs,const double DropLimit)
{
  const long unsigned Nr = reg->row;
  const long unsigned Nc = reg->col;
  double **const m = reg->reg->A;
  long *Ap   = calloc(Nc+1,sizeof(*Ap));
  long *Ai   = 0;
  double *Ax = 0;
  long tNN0 = 0;/* total number of none zero entries */
  long NN0;/* number of none zero entries in each columnt */
  long unsigned r,c;/* row and column */
  
  for (c = 0; c < Nc; ++c)
  {
    NN0 = 0;
    for (r = 0; r < Nr; ++r)
    {
      if (GRT(m[r][c],DropLimit))
      {
        Ai = realloc(Ai,(long unsigned)(Ap[c]+NN0+1)*sizeof(*Ai));
        pointerEr(Ai);
        Ax = realloc(Ax,(long unsigned)(Ap[c]+NN0+1)*sizeof(*Ax));
        pointerEr(Ax);
        Ai[Ap[c]+NN0] = (long)r;
        Ax[Ap[c]+NN0] = m[r][c];
        NN0++;
        tNN0++;
      }
    }
    Ap[c+1] = tNN0;
  }
  
  ccs->ccs->Ap = Ap;
  ccs->ccs->Ai = Ai;
  ccs->ccs->Ax = Ax;
}
