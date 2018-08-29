/*
// Alireza Rashti
// August 2018
*/

#include "solver_tests.h"

/* testing various solvers */
void solver_tests(void)
{
  int status;
  
  if (DO)
  {
    printf("Testing umfpack_di solver:");
    status = test_solver_umfpack_di();
    check_test_result(status);
  }
  
  
}

/* solving a generic system of equation a.x = b 
// to test umfpack_di+* solver.
*/
static int test_solver_umfpack_di(void)
{
  UmfPack_T eg[1];
  Matrix_T *A = alloc_matrix(REG_SF,4,4);
  Matrix_T *ccs;
  double **const a = A->reg->A;
  double *b = alloc_double(4);
  double *x = alloc_double(4);
  double *ans = alloc_double(4);
  Flag_T flg = NONE;
  int i;
  
  /* filling system with a known answer one */
  a[0][0] = 1;
  a[0][1] = -1;
  a[0][2] = 2;
  a[0][3] = -1;
  a[1][0] = 2;
  a[1][1] = -2;
  a[1][2] = 3;
  a[1][3] = -3;
  a[2][0] = 1;
  a[2][1] = 1;
  a[2][2] = 1;
  a[2][3] = 0;
  a[3][0] = 1;
  a[3][1] = -1;
  a[3][2] = 4;
  a[3][3] = 3;
  
  b[0] = -8;
  b[1] = -20;
  b[2] = -2;
  b[3] = 4;
  
  ans[0] = -7;
  ans[1] = 3;
  ans[2] = 2;
  ans[3] = 2;
  
  ccs = cast_matrix_ccs(A);
  /* filling umfpack input */
  eg->a = ccs;
  eg->b = b;
  eg->x = x;
  /* solving */
  direct_solver_umfpack_di(eg);
  
  /* comparing x vs ans*/
  flg = NONE;
  for (i = 0 ; i < 4; ++i)
    if (!EQL(x[i],ans[i]))
    {
      flg = FOUND;
      break;
    }
  
  free_matrix(A);
  free_matrix(ccs);
  free(x);
  free(ans);
  free(b);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}
