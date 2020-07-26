/*
// Alireza Rashti
// August 2018
*/

#include "solver_tests.h"

/* testing various solvers */
void solver_tests(void)
{
  int status;
  
  if (DO_NOT)
  {
    status = test_solver_umfpack_di();
    printf("Testing umfpack_di solver:");
    check_test_result(status);
  }
  if (DO)
  {
    status = test_solver_series_umfpack_di();
    printf("Testing umfpack_series_di solver:");
    check_test_result(status);
  }
  
  
}

/* solving a generic system of equation a.x = b 
// to test direct_solver_umfpack_di solver. */
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

/* solving a generic system of equation a.x = b
// for a series of b. it tests 3 b's.
// to test direct_solver_series_umfpack_di */
static int test_solver_series_umfpack_di(void)
{
  UmfPack_T eg[1];
  Matrix_T *A = alloc_matrix(REG_SF,4,4);
  Matrix_T *ccs;
  const int N = 3;
  const unsigned int dim = 4;
  double **const a = A->reg->A;
  double *b[N];
  double *x[N];
  double *ans[N];
  Flag_T flg = NONE;
  int i,n;
  
  /* filling system with a known answer one */
  a[0][0] = 1.1;
  a[0][1] = -1;
  a[0][2] = 2;
  a[0][3] = -1;
  a[1][0] = 2.1;
  a[1][1] = -2;
  a[1][2] = 3.7;
  a[1][3] = -3;
  a[2][0] = 1;
  a[2][1] = 1;
  a[2][2] = 1.9;
  a[2][3] = 0;
  a[3][0] = 1;
  a[3][1] = -1;
  a[3][2] = 4.3;
  a[3][3] = 3.6;
  
  for (i = 0; i < N; i++)
  {
    b[i]   = alloc_double(dim);
    x[i]   = alloc_double(dim);
    ans[i] = alloc_double(dim);
  }
  
  /* N = 0 */
  b[0][0] = -8;
  b[0][1] = -20;
  b[0][2] = -2;
  b[0][3] = 4;
  
  ans[0][0] = 4.881170018281531;
  ans[0][1] = 0.6910420475319937;
  ans[0][2] = -3.985374771480802;
  ans[0][3] = 4.7074954296160865;
  
  /* N = 1 */
  b[1][0] = -4.4;
  b[1][1] = -30.9;
  b[1][2] = 2.6;
  b[1][3] = 300.67;
  
  ans[1][0] = -190.9611517367457;
  ans[1][1] = -5.673994515539325;
  ans[1][2] = 104.86060329067634;
  ans[1][3] = 9.737934186471705;
  
  /* N = 2 */
  b[2][0] = 4.4;
  b[2][1] = 30.9;
  b[2][2] = -2.6;
  b[2][3] = -10.67;
  
  ans[2][0] = -87.3752285191955;
  ans[2][1] = 6.999405850091391;
  ans[2][2] = 40.93464351005479;
  ans[2][3] = -25.642870201096866;
  
  ccs = cast_matrix_ccs(A);
  /* filling umfpack input */
  eg->a = ccs;
  eg->bs = b;
  eg->xs = x;
  eg->ns = (unsigned)N;
  /* solving */
  direct_solver_series_umfpack_di(eg);
  
  /* comparing x vs ans*/
  flg = NONE;
  for (n = 0; n < N; ++n)
  {
    for (i = 0 ; i < 4; ++i)
      if (!EQL(x[n][i],ans[n][i]))
      {
        flg = FOUND;
        //test
        printf("answer[%d](numeric,analytic) = (%0.15f,%0.15f)\n"
          ,n,x[n][i],ans[n][i]);
        //end
        break;
      }
  }
  
  free_matrix(A);
  free_matrix(ccs);
  for (n = 0; n < N; ++n)
  {
    free(x[n]);
    free(ans[n]);
    free(b[n]);
  }
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}
