/*
// Alireza Rashti
// September 2018
*/

#include "preconditioners.h"

/* preconditioning a.x = b equations.
// NOTE: the implementaion is very crude and its only for testes!
*/
void precondition(Matrix_T *const a,double *const b)
{
  precondition_GS(a,b);
}

/* Gauss-Seidel preconditioned */
static void precondition_GS(Matrix_T *const A,double *const b)
{
  Matrix_T *L = alloc_matrix(REG_SF,A->row,A->col);
  Matrix_T *invL = 0;
  double **a = A->reg->A;
  double **l = L->reg->A;
  double **invl;
  long i,j,k;
  
  for(i = 0; i < A->row; i++)
  {
    for (j = 0; j <= i; ++j)
      l[i][j] = a[i][j];
  }
  
  invL = invert_matrix(L);
  invl = invL->reg->A;
  
  /* M^-1 * a.x = M^-1 b */
  for(i = 0; i < A->row; i++)
  {
    double aij;
    double b_new;
    
    for (j = 0; j < A->col; ++j)
    {
      aij = 0;
      for (k = 0; k < A->row; ++k)
        aij += invl[i][k]*a[k][j];
      
      a[i][j] = aij;
    }
    
    b_new = 0;
    for (k = 0; k < A->col; ++k)
      b_new += invl[i][k]*b[k];
      
    b[i] = b_new;
  }
  
  free_matrix(invL);
  free_matrix(L);
}
