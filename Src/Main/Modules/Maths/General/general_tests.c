/*
// Alireza Rashti
// Feb 2018
*/

#include "general_tests.h"

/* testing functions used for summing various quantities */
void summation_tests(void)
{
  int status;

  if(DO)
  {
    status = various_sums_of_cos_test();
    printf("Various sums of cosine functions:");
    check_test_result(status);
  }

}

/* tests:
// o. sum_1_N_cos_ia
// o. d_dq_sum_1_N_cos_ixb_cos_ixa
// ->return value: TEST_SUCCESSFUL or TEST_UNSUCCESSFUL */
static int various_sums_of_cos_test(void)
{
  const unsigned M = 40;/* number of different numbers */
  const unsigned N = (unsigned)floor(random_double(10,60,0));
  double a[M],sum_cal,sum_real;
  Flag_T flg_sum_1_N_cos_ia;
  unsigned i;
  
  /* filling a */
  a[0] = 0;
  a[1] = 2*M_PI;
  for (i = 2; i < M; ++i)
    a[i] = random_double(0,2*M_PI,i-2);
  
  /* testing sum_1_N_cos_ia */
  flg_sum_1_N_cos_ia = NONE;
  for (i = 0; i < M; ++i)
  {
    unsigned j;
    /* real calculation */
    sum_real = 0;
    
    for (j = 1; j <= N; ++j)
      sum_real += cos(j*a[i]);
    sum_cal = sum_1_N_cos_ia(N,a[i]);
    
    if (!EQL(sum_real,sum_cal))
    {
      flg_sum_1_N_cos_ia = FOUND;
      break;
    }
    
  }
  
  if (flg_sum_1_N_cos_ia == FOUND)
  {
    printf("sum_1_N_cos_ia function failed!\n");
    return TEST_UNSUCCESSFUL;
  }
  
  return TEST_SUCCESSFUL;
}
