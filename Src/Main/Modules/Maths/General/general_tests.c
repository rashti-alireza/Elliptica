/*
// Alireza Rashti
// Feb 2018
*/

#include "general_tests.h"

/* testing functions used for summing various quantities */
void summation_tests(void)
{
  int status;
  
  if(DO_NOT)
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
  const Uint M = 60;/* number of different numbers MUST BE > 4 */
  const Uint N = (Uint)floor(random_double(10,60,0));
  double a[M],b[M],sum_cal,sum_real;
  Flag_T flg_sum_1_N_cos_ia,flg_d_dq_sum_1_N_cos_ixb_cos_ixa;
  Uint i,j,k;
  
  /* filling a */
  a[0] = 0;
  a[1] = 2*M_PI;
  a[2] = M_PI;
  a[3] = -M_PI;
  for (i = 4; i < M; ++i)
    a[i] = random_double(-M_PI,2*M_PI,i-4);
  
  /* testing sum_1_N_cos_ia */
  flg_sum_1_N_cos_ia = NONE;
  for (i = 0; i < M; ++i)
  {
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
  
  /* filling a and b */
  a[0] = 0;
  a[1] = M_PI;
  a[2] = M_PI/3;
  a[3] = M_PI/11;
  b[0] = 0;
  b[1] = M_PI;
  b[2] = M_PI/3;
  b[3] = M_PI/11;
  
  for (i = 4; i < M; ++i)
  {
    a[i] = random_double(0,M_PI,i-4);
    b[i] = random_double(0,M_PI,i);
  }
  
  /* testing d_dq_sum_1_N_cos_ixb_cos_ixa */
  flg_d_dq_sum_1_N_cos_ixb_cos_ixa = NONE;
  for (i = 0; i < M; ++i)
  {
    for (j = 0; j < M; ++j)
    {
      /* real calculation */
      sum_real = 0;
      
      for (k = 1; k <= N; ++k)
      {
        double dT_dx = k*Cheb_Un((int)k-1,cos(a[i]));
        sum_real += cos(k*b[j])*dT_dx;
      }
        
      sum_cal = d_dq_sum_1_N_cos_ixb_cos_ixa((int)N,b[j],a[i]);
      
      if (!EQL(ABSd(sum_real-sum_cal)/MaxMag_d(sum_real,sum_cal),0))
      {
        /* more details */
        printf("Discrepancy as follows:\n"
        "(N,a,b,real,cal,diff)=(%u,%f,%f,%f,%f,%0.15f)\n"
        "Note:\nWhen using summation method a huge round off error incurred.\n"
          ,N,a[i],b[j],sum_real,sum_cal,sum_real-sum_cal);
        flg_d_dq_sum_1_N_cos_ixb_cos_ixa = FOUND;
        break;
      }
      
    }
    if (flg_d_dq_sum_1_N_cos_ixb_cos_ixa == FOUND)
      break;
  }
  
  if (flg_sum_1_N_cos_ia == FOUND)
  {
    printf("sum_1_N_cos_ia function failed!\n");
    return TEST_UNSUCCESSFUL;
  }
  if (flg_d_dq_sum_1_N_cos_ixb_cos_ixa == FOUND)
  {
    printf("d_dq_sum_1_N_cos_ixb_cos_ixa function failed!\n");
    return TEST_UNSUCCESSFUL;
  }
  
  return TEST_SUCCESSFUL;
}
