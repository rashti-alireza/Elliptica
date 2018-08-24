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
    convert_reg2ccs(m,ccs);
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

/* converting a regular storag format to ccs format */
static void convert_reg2ccs(const Matrix_T *const reg,Matrix_T *const ccs)
{
}
