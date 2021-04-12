/*
// Alireza Rashti
// Feb 2018
*/

#include "matrix_memory_managers.h"

/* given row and column of a asked matrix it allocates matrix 
// according to the given type using calloc and return the result.
// note : if either of col or row is zero, it returns null.
// ->return value: new empty matrix. null if something not defined correctly. */
Matrix_T *alloc_matrix(const Matrix_SF_T type_e,const long row,const long col)
{
  Matrix_T *m = 0;
  
  /* returns null if either of row or col is zero */
  if (row == 0 || col == 0)
    return m;
  
  if (type_e == UNDEF_SF)
    return m;
    
  m = calloc(1,sizeof(*m));
  IsNull(m);
  m->row = row;
  m->col = col;
  
  switch(type_e)
  {
    case REG_SF:
      m->reg_f = 1;
      m->reg->A = alloc_2D_double((long Uint)row,(long Uint)col);
      break;
    case RMO_SF:
      m->rmo_f  = 1;
      m->rmo->A = alloc_double((Uint)(row*col));
      break;
    case TRI_SF:
      m->tri_f = 1;
      break;
    case CCS_SF:
      m->ccs_f = 1;
      break;
    case CRS_SF:
      m->crs_f = 1;
      break;
    case TRI_L_SF:
      m->tri_l_f = 1;
      break;
    case CCS_L_SF:
      m->ccs_l_f = 1;
      break;
    case CRS_L_SF:
      m->crs_l_f = 1;
      break;
    default:
      Error0("The specified type is undefined.\n");
  }
  
  return m;
}

/* freeing matrix memroy */
void free_matrix(Matrix_T *m)
{
  if (!m)
    return;
    
  if (m->reg_f)
  {
    if (m->reg->A)
      free_2d_mem(m->reg->A,(long Uint)m->row);
  }
  else if (m->rmo_f)
  {
    if (m->rmo->A)
      free(m->rmo->A);
  }
  else if (m->tri_f)
  {
    if (m->tri->row)
      free(m->tri->row);
    if (m->tri->col)
      free(m->tri->col);
    if (m->tri->a)
      free(m->tri->a); 
  }
  else if (m->ccs_f)
  {
    if (m->ccs->Ap)
      free(m->ccs->Ap);
    if (m->ccs->Ai)
      free(m->ccs->Ai);
    if (m->ccs->Ax)
      free(m->ccs->Ax);
    if (m->ccs->Ap_cg)
      free(m->ccs->Ap_cg);
    if (m->ccs->i_cg)
      free(m->ccs->i_cg);
  }
  else if (m->crs_f)
  {
    if (m->crs->Ap)
      free(m->crs->Ap);
    if (m->crs->Aj)
      free(m->crs->Aj);
    if (m->crs->Ax)
      free(m->crs->Ax);
  }
  else if (m->tri_l_f)
  {
    if (m->tri_long->row)
      free(m->tri_long->row);
    if (m->tri_long->col)
      free(m->tri_long->col);
    if (m->tri_long->a)
      free(m->tri_long->a); 
  }
  else if (m->ccs_l_f)
  {
    if (m->ccs_long->Ap)
      free(m->ccs_long->Ap);
    if (m->ccs_long->Ai)
      free(m->ccs_long->Ai);
    if (m->ccs_long->Ax)
      free(m->ccs_long->Ax);
  }
  else if (m->crs_l_f)
  {
    if (m->crs_long->Ap)
      free(m->crs_long->Ap);
    if (m->crs_long->Aj)
      free(m->crs_long->Aj);
    if (m->crs_long->Ax)
      free(m->crs_long->Ax);
  }
  else
    Error0("No matrix format is defined for this given matrix.\n");
    
  free(m);
}
