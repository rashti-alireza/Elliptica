/*
// Alireza Rashti
// August 2018
*/

#include "linear_algebra_tests.h"
#define ROW 400
#define COL 300

/* testing various matrix routines */
void matrix_tests(void)
{
  int status;
  
  if (NOT_DO)
  {
    printf("Casting of regular matrix to CCS one:");
    status = cast_matrix_ccs_test();
    check_test_result(status);
  }
  if (DO)
  {
    printf("Reading a CCS matrix directly one:");
    status = read_ccs_test();
    check_test_result(status);
  } 
}

/* test: reading a CCS matrix directly rather casting 
// it to regular format.
// ->return value: TEST_SUCCESSFUL or TEST_UNSUCCESSFUL
*/
static int read_ccs_test(void)
{
  const long Nr = ROW;
  const long Nc = COL;
  Matrix_T *reg = make_generic_matrix(Nr,Nc); 
  Matrix_T *ccs  = cast_matrix_ccs(reg);
  double **const m = reg->reg->A;
  long r,c;
  Flag_T flg = NONE;
  
  for (r = 0; r < Nr; ++r)
  {
    for (c = 0; c < Nc; ++c)
      if (!EQL(m[r][c],read_matrix_entry_ccs(ccs,r,c)))
      {
        flg = FOUND;
        break;
      }
      
    if (flg == FOUND)
      break;
  }
  
  
  free_matrix(reg);
  free_matrix(ccs);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* testing function "cast_matrix_ccs" 
// ->return value: TEST_SUCCESSFUL or TEST_UNSUCCESSFUL.
*/
static int cast_matrix_ccs_test(void)
{
  const long Nr = ROW;
  const long Nc = COL;
  Matrix_T *reg = make_generic_matrix(Nr,Nc);
  Matrix_T *ccs = 0;
  Matrix_T *reg2 = 0;
  double **m1,**m2;
  long r,c;
  Flag_T flg;
  
  ccs  = cast_matrix_ccs(reg);
  reg2 = cast_matrix_reg(ccs);
  
  m1 = reg->reg->A;
  m2 = reg2->reg->A;
  
  flg = NONE;
  for (r = 0; r < Nr; ++r)
  {
    for (c = 0; c < Nc; ++c)
      if (!EQL(m1[r][c],m2[r][c]))
      {
        flg = FOUND;
        break;
      }
      
    if (flg == FOUND)
      break;
  }
  
  free_matrix(ccs);
  free_matrix(reg);
  free_matrix(reg2);
  
  if (flg == FOUND)
    return TEST_UNSUCCESSFUL;
  
  return TEST_SUCCESSFUL;
}

/* making a generic matrix at regular format for test purposes.
// ->return value: made matrix.
*/
static Matrix_T *make_generic_matrix(const long Nr,const long Nc)
{
  Matrix_T *reg = alloc_matrix(REG_SF,Nr,Nc);
  double **const A = reg->reg->A;
  long r,c;
  
  /* making a generic matrix */
  for (r = 0; r < Nr; ++r)
  {
    for (c = 0; c < Nc; ++c)
    {
      A[r][c] = (double)Nr*sin((double)c)+(double)Nr*cos((double)r);
      if ((r+c)%4 == 1)
        A[r][c] = 0;
    }
  }
  
  return reg;
}