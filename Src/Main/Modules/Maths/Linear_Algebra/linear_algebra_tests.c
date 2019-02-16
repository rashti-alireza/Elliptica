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
  
  if(NOT_DO)
  {
    status = cast_matrix_ccs_test();
    printf("Casting of regular matrix to CCS format:");
    check_test_result(status);
  }
  if (NOT_DO)
  {
    status = read_ccs_test();
    printf("Testing the reader of CCS matrix:");
    check_test_result(status);
  }
  if (DO)
  {
    status = matrices_arithmetic_test();
    printf("Testing matrices arithmetic:");
    check_test_result(status);
  }
}

/* test: matrices arithmetic like multiplication and etc.
// ->return value: TEST_SUCCESSFUL or TEST_UNSUCCESSFUL */
static int matrices_arithmetic_test(void)
{
  const long Nr = 4;
  const long Nc = 4;
  Matrix_T *M1 = alloc_matrix(REG_SF,Nr,Nc);
  Matrix_T *M2 = alloc_matrix(REG_SF,Nr,Nc);
  Matrix_T *C  = alloc_matrix(REG_SF,Nr,Nc);
  Matrix_T *numeric_C;
  double **const m1 = M1->reg->A;
  double **const m2 = M2->reg->A;
  double **const c  = C->reg->A;
  double **numeric_c;
  double v[Nr],mv[Nr],numeric_mv[Nr];
  Flag_T flg_mv,flg_regxreg_ab,flg_ccsxccs_ab,flg_ccsxccs_aTransB;/* flags for various tests */
  long i,j;
  
  /* fill up some arbitrary matrices */
  m1[0][0] = 1.1;
  m1[0][1] = -1;
  m1[0][2] = 2;
  m1[0][3] = -1;
  m1[1][0] = 2.1;
  m1[1][1] = -2;
  m1[1][2] = 3.7;
  m1[1][3] = -3;
  m1[2][0] = 1;
  m1[2][1] = 1;
  m1[2][2] = 1.9;
  m1[2][3] = 0;
  m1[3][0] = 1;
  m1[3][1] = -1;
  m1[3][2] = 4.3;
  m1[3][3] = 3.6;
  
  m2[0][0] = 5; 
  m2[0][1] = 8;
  m2[0][2] = 80;
  m2[0][3] = 800;
  m2[1][0] = 1.4;
  m2[1][1] = 66.98;
  m2[1][2] = 24;
  m2[1][3] = 241.9;
  m2[2][0] = 55.888;
  m2[2][1] = 231;
  m2[2][2] = 23;
  m2[2][3] = 10.5;
  m2[3][0] = 0;
  m2[3][1] = 454;
  m2[3][2] = 6.7;
  m2[3][3] = 0.001;
  
  /* c = m1.m2 */
  c[0][0] = 115.87599999999999;
  c[0][1] = -50.18000000000001;
  c[0][2] = 103.3;
  c[0][3] = 659.0990000000002;
  c[1][0] = 214.4856;
  c[1][1] = -624.4599999999999;
  c[1][2] = 185.00000000000003;
  c[1][3] = 1235.047;
  c[2][0] = 112.5872;
  c[2][1] = 513.88;
  c[2][2] = 147.7;
  c[2][3] = 1061.8500000000001;
  c[3][0] = 243.91839999999996;
  c[3][1] = 2568.7200000000003;
  c[3][2] = 179.01999999999998;
  c[3][3] = 603.2536;
  
  /* an arbitrary vector */
  v[0] = 1.3;
  v[1] = -0.4555;
  v[2] = -200.222;
  v[3] = 65.454;
  
  /* mv = m1.v */
  mv[0] = -464.01250000000005;
  mv[1] = -933.5424;
  mv[2] = -379.57730000000004;
  mv[3] = -623.5647000000001;
  
  /* test matrix_by_vector funciont */
  matrix_by_vector(M1,v,numeric_mv,INITIALIZE);
  flg_mv = NONE;
  for (i = 0; i < Nr; ++i)
    if (!EQL(mv[i],numeric_mv[i]))
    {
      flg_mv = FOUND;
      break;
    }

  /* test matrix_by_matrix function (regular*regular) and "a*b" */
  numeric_C = matrix_by_matrix(M1,M2,"a*b");
  numeric_c = numeric_C->reg->A;
  flg_regxreg_ab = NONE;
  for (i = 0; i < Nr; ++i)
  {
    for (j = 0; j < Nc; ++j)
      if (!EQL(c[i][j],numeric_c[i][j]))
      {
        flg_regxreg_ab = FOUND;
        break;
      }
    if (flg_regxreg_ab == FOUND)
      break;
  }
  free_matrix(numeric_C);
  
  /* test matrix_by_matrix function (CCS*CCS) and "a*b" */
  flg_ccsxccs_ab = NONE;
  Matrix_T *M1_ccs = cast_matrix_ccs(M1);
  Matrix_T *M2_ccs = cast_matrix_ccs(M2);
  numeric_C = matrix_by_matrix(M1_ccs,M2_ccs,"a*b");
  numeric_c = numeric_C->reg->A;
  for (i = 0; i < Nr; ++i)
  {
    for (j = 0; j < Nc; ++j)
      if (!EQL(c[i][j],numeric_c[i][j]))
      {
        flg_ccsxccs_ab = FOUND;
        printf("matrix_by_matrix function for (ccs*ccs) "
            "and directive \"a*b\" failed.\n");
        break;
      }
    if (flg_ccsxccs_ab == FOUND)
      break;
  }
  free_matrix(numeric_C);
  free_matrix(M1_ccs);
  free_matrix(M2_ccs);
  
  /* test matrix_by_matrix function (reg*reg) and "a*Transpose(b)" */
  /* c = m1 * Trans(m2) */
  c[0][0] = -642.5;
  c[0][1] = -259.34000000000003;
  c[0][2] = -134.0232;
  c[0][3] = -440.601;
  c[1][0] = -2109.5;
  c[1][1] = -767.9200000000001;
  c[1][2] = -291.0352;
  c[1][3] = -883.2130000000001;
  c[2][0] = 165.0;
  c[2][1] = 113.98;
  c[2][2] = 330.58799999999997;
  c[2][3] = 466.73;
  c[3][0] = 3221.0;
  c[3][1] = 908.46;
  c[3][2] = -38.412;
  c[3][3] = -425.1864;
  flg_ccsxccs_aTransB = NONE;
  numeric_C = matrix_by_matrix(M1,M2,"a*Transpose(b)");
  numeric_c = numeric_C->reg->A;
  for (i = 0; i < Nr; ++i)
  {
    for (j = 0; j < Nc; ++j)
      if (!EQL(c[i][j],numeric_c[i][j]))
      {
        flg_ccsxccs_aTransB = FOUND;
        break;
      }
    if (flg_ccsxccs_aTransB == FOUND)
      break;
  }
  free_matrix(numeric_C);
  
  free_matrix(M1);
  free_matrix(M2);
  free_matrix(C);
  
  /* results of tests */
  if (flg_mv == FOUND)
  {
    printf("matrix_by_vector function failed.\n");
    return TEST_UNSUCCESSFUL;
  }
  if (flg_regxreg_ab == FOUND)
  {
    printf("matrix_by_matrix function for (reg*reg) "
            "and directive \"a*b\" function failed.\n");
    return TEST_UNSUCCESSFUL;
  }
  if (flg_ccsxccs_ab == FOUND)
  {
    printf("matrix_by_matrix function for (ccs*ccs) "
            "and directive \"a*b\" failed.\n");
    return TEST_UNSUCCESSFUL;
  }
  if (flg_ccsxccs_aTransB == FOUND)
  {
    printf("matrix_by_matrix function for (ccs*ccs) "
            "and directive \"a*Transpose(b)\" failed.\n");
    return TEST_UNSUCCESSFUL;
  }

  return TEST_SUCCESSFUL;
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
      A[r][c] = random_double((double)(Nr-Nc),(double)(Nr+Nc),(unsigned)r);
      if ((r+c)%4 == 1)
        A[r][c] = 0;
    }
  }
  
  return reg;
}
