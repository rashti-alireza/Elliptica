 /* direct output of this file to spherical_harmonic_Ylm.c */
#include "core_lib.h"
#include "maths_analytic_lib.h"
#include "maths_general_lib.h"
#include <complex.h>

static const int lmax = 15;
typedef double complex fYlm_T (const double theta, const double phi);

double complex Ylm(const int l, const int m, const double theta, const double phi);
double complex Ylm_l0m0(const double theta, const double phi);
double complex Ylm_l1m0(const double theta, const double phi);
double complex Ylm_l1m1(const double theta, const double phi);
double complex Ylm_l2m0(const double theta, const double phi);
double complex Ylm_l2m1(const double theta, const double phi);
double complex Ylm_l2m2(const double theta, const double phi);
double complex Ylm_l3m0(const double theta, const double phi);
double complex Ylm_l3m1(const double theta, const double phi);
double complex Ylm_l3m2(const double theta, const double phi);
double complex Ylm_l3m3(const double theta, const double phi);
double complex Ylm_l4m0(const double theta, const double phi);
double complex Ylm_l4m1(const double theta, const double phi);
double complex Ylm_l4m2(const double theta, const double phi);
double complex Ylm_l4m3(const double theta, const double phi);
double complex Ylm_l4m4(const double theta, const double phi);
double complex Ylm_l5m0(const double theta, const double phi);
double complex Ylm_l5m1(const double theta, const double phi);
double complex Ylm_l5m2(const double theta, const double phi);
double complex Ylm_l5m3(const double theta, const double phi);
double complex Ylm_l5m4(const double theta, const double phi);
double complex Ylm_l5m5(const double theta, const double phi);
double complex Ylm_l6m0(const double theta, const double phi);
double complex Ylm_l6m1(const double theta, const double phi);
double complex Ylm_l6m2(const double theta, const double phi);
double complex Ylm_l6m3(const double theta, const double phi);
double complex Ylm_l6m4(const double theta, const double phi);
double complex Ylm_l6m5(const double theta, const double phi);
double complex Ylm_l6m6(const double theta, const double phi);
double complex Ylm_l7m0(const double theta, const double phi);
double complex Ylm_l7m1(const double theta, const double phi);
double complex Ylm_l7m2(const double theta, const double phi);
double complex Ylm_l7m3(const double theta, const double phi);
double complex Ylm_l7m4(const double theta, const double phi);
double complex Ylm_l7m5(const double theta, const double phi);
double complex Ylm_l7m6(const double theta, const double phi);
double complex Ylm_l7m7(const double theta, const double phi);
double complex Ylm_l8m0(const double theta, const double phi);
double complex Ylm_l8m1(const double theta, const double phi);
double complex Ylm_l8m2(const double theta, const double phi);
double complex Ylm_l8m3(const double theta, const double phi);
double complex Ylm_l8m4(const double theta, const double phi);
double complex Ylm_l8m5(const double theta, const double phi);
double complex Ylm_l8m6(const double theta, const double phi);
double complex Ylm_l8m7(const double theta, const double phi);
double complex Ylm_l8m8(const double theta, const double phi);
double complex Ylm_l9m0(const double theta, const double phi);
double complex Ylm_l9m1(const double theta, const double phi);
double complex Ylm_l9m2(const double theta, const double phi);
double complex Ylm_l9m3(const double theta, const double phi);
double complex Ylm_l9m4(const double theta, const double phi);
double complex Ylm_l9m5(const double theta, const double phi);
double complex Ylm_l9m6(const double theta, const double phi);
double complex Ylm_l9m7(const double theta, const double phi);
double complex Ylm_l9m8(const double theta, const double phi);
double complex Ylm_l9m9(const double theta, const double phi);
double complex Ylm_l10m0(const double theta, const double phi);
double complex Ylm_l10m1(const double theta, const double phi);
double complex Ylm_l10m2(const double theta, const double phi);
double complex Ylm_l10m3(const double theta, const double phi);
double complex Ylm_l10m4(const double theta, const double phi);
double complex Ylm_l10m5(const double theta, const double phi);
double complex Ylm_l10m6(const double theta, const double phi);
double complex Ylm_l10m7(const double theta, const double phi);
double complex Ylm_l10m8(const double theta, const double phi);
double complex Ylm_l10m9(const double theta, const double phi);
double complex Ylm_l10m10(const double theta, const double phi);
double complex Ylm_l11m0(const double theta, const double phi);
double complex Ylm_l11m1(const double theta, const double phi);
double complex Ylm_l11m2(const double theta, const double phi);
double complex Ylm_l11m3(const double theta, const double phi);
double complex Ylm_l11m4(const double theta, const double phi);
double complex Ylm_l11m5(const double theta, const double phi);
double complex Ylm_l11m6(const double theta, const double phi);
double complex Ylm_l11m7(const double theta, const double phi);
double complex Ylm_l11m8(const double theta, const double phi);
double complex Ylm_l11m9(const double theta, const double phi);
double complex Ylm_l11m10(const double theta, const double phi);
double complex Ylm_l11m11(const double theta, const double phi);
double complex Ylm_l12m0(const double theta, const double phi);
double complex Ylm_l12m1(const double theta, const double phi);
double complex Ylm_l12m2(const double theta, const double phi);
double complex Ylm_l12m3(const double theta, const double phi);
double complex Ylm_l12m4(const double theta, const double phi);
double complex Ylm_l12m5(const double theta, const double phi);
double complex Ylm_l12m6(const double theta, const double phi);
double complex Ylm_l12m7(const double theta, const double phi);
double complex Ylm_l12m8(const double theta, const double phi);
double complex Ylm_l12m9(const double theta, const double phi);
double complex Ylm_l12m10(const double theta, const double phi);
double complex Ylm_l12m11(const double theta, const double phi);
double complex Ylm_l12m12(const double theta, const double phi);
double complex Ylm_l13m0(const double theta, const double phi);
double complex Ylm_l13m1(const double theta, const double phi);
double complex Ylm_l13m2(const double theta, const double phi);
double complex Ylm_l13m3(const double theta, const double phi);
double complex Ylm_l13m4(const double theta, const double phi);
double complex Ylm_l13m5(const double theta, const double phi);
double complex Ylm_l13m6(const double theta, const double phi);
double complex Ylm_l13m7(const double theta, const double phi);
double complex Ylm_l13m8(const double theta, const double phi);
double complex Ylm_l13m9(const double theta, const double phi);
double complex Ylm_l13m10(const double theta, const double phi);
double complex Ylm_l13m11(const double theta, const double phi);
double complex Ylm_l13m12(const double theta, const double phi);
double complex Ylm_l13m13(const double theta, const double phi);
double complex Ylm_l14m0(const double theta, const double phi);
double complex Ylm_l14m1(const double theta, const double phi);
double complex Ylm_l14m2(const double theta, const double phi);
double complex Ylm_l14m3(const double theta, const double phi);
double complex Ylm_l14m4(const double theta, const double phi);
double complex Ylm_l14m5(const double theta, const double phi);
double complex Ylm_l14m6(const double theta, const double phi);
double complex Ylm_l14m7(const double theta, const double phi);
double complex Ylm_l14m8(const double theta, const double phi);
double complex Ylm_l14m9(const double theta, const double phi);
double complex Ylm_l14m10(const double theta, const double phi);
double complex Ylm_l14m11(const double theta, const double phi);
double complex Ylm_l14m12(const double theta, const double phi);
double complex Ylm_l14m13(const double theta, const double phi);
double complex Ylm_l14m14(const double theta, const double phi);
double complex Ylm_l1m_1(const double theta, const double phi);
double complex Ylm_l2m_1(const double theta, const double phi);
double complex Ylm_l2m_2(const double theta, const double phi);
double complex Ylm_l3m_1(const double theta, const double phi);
double complex Ylm_l3m_2(const double theta, const double phi);
double complex Ylm_l3m_3(const double theta, const double phi);
double complex Ylm_l4m_1(const double theta, const double phi);
double complex Ylm_l4m_2(const double theta, const double phi);
double complex Ylm_l4m_3(const double theta, const double phi);
double complex Ylm_l4m_4(const double theta, const double phi);
double complex Ylm_l5m_1(const double theta, const double phi);
double complex Ylm_l5m_2(const double theta, const double phi);
double complex Ylm_l5m_3(const double theta, const double phi);
double complex Ylm_l5m_4(const double theta, const double phi);
double complex Ylm_l5m_5(const double theta, const double phi);
double complex Ylm_l6m_1(const double theta, const double phi);
double complex Ylm_l6m_2(const double theta, const double phi);
double complex Ylm_l6m_3(const double theta, const double phi);
double complex Ylm_l6m_4(const double theta, const double phi);
double complex Ylm_l6m_5(const double theta, const double phi);
double complex Ylm_l6m_6(const double theta, const double phi);
double complex Ylm_l7m_1(const double theta, const double phi);
double complex Ylm_l7m_2(const double theta, const double phi);
double complex Ylm_l7m_3(const double theta, const double phi);
double complex Ylm_l7m_4(const double theta, const double phi);
double complex Ylm_l7m_5(const double theta, const double phi);
double complex Ylm_l7m_6(const double theta, const double phi);
double complex Ylm_l7m_7(const double theta, const double phi);
double complex Ylm_l8m_1(const double theta, const double phi);
double complex Ylm_l8m_2(const double theta, const double phi);
double complex Ylm_l8m_3(const double theta, const double phi);
double complex Ylm_l8m_4(const double theta, const double phi);
double complex Ylm_l8m_5(const double theta, const double phi);
double complex Ylm_l8m_6(const double theta, const double phi);
double complex Ylm_l8m_7(const double theta, const double phi);
double complex Ylm_l8m_8(const double theta, const double phi);
double complex Ylm_l9m_1(const double theta, const double phi);
double complex Ylm_l9m_2(const double theta, const double phi);
double complex Ylm_l9m_3(const double theta, const double phi);
double complex Ylm_l9m_4(const double theta, const double phi);
double complex Ylm_l9m_5(const double theta, const double phi);
double complex Ylm_l9m_6(const double theta, const double phi);
double complex Ylm_l9m_7(const double theta, const double phi);
double complex Ylm_l9m_8(const double theta, const double phi);
double complex Ylm_l9m_9(const double theta, const double phi);
double complex Ylm_l10m_1(const double theta, const double phi);
double complex Ylm_l10m_2(const double theta, const double phi);
double complex Ylm_l10m_3(const double theta, const double phi);
double complex Ylm_l10m_4(const double theta, const double phi);
double complex Ylm_l10m_5(const double theta, const double phi);
double complex Ylm_l10m_6(const double theta, const double phi);
double complex Ylm_l10m_7(const double theta, const double phi);
double complex Ylm_l10m_8(const double theta, const double phi);
double complex Ylm_l10m_9(const double theta, const double phi);
double complex Ylm_l10m_10(const double theta, const double phi);
double complex Ylm_l11m_1(const double theta, const double phi);
double complex Ylm_l11m_2(const double theta, const double phi);
double complex Ylm_l11m_3(const double theta, const double phi);
double complex Ylm_l11m_4(const double theta, const double phi);
double complex Ylm_l11m_5(const double theta, const double phi);
double complex Ylm_l11m_6(const double theta, const double phi);
double complex Ylm_l11m_7(const double theta, const double phi);
double complex Ylm_l11m_8(const double theta, const double phi);
double complex Ylm_l11m_9(const double theta, const double phi);
double complex Ylm_l11m_10(const double theta, const double phi);
double complex Ylm_l11m_11(const double theta, const double phi);
double complex Ylm_l12m_1(const double theta, const double phi);
double complex Ylm_l12m_2(const double theta, const double phi);
double complex Ylm_l12m_3(const double theta, const double phi);
double complex Ylm_l12m_4(const double theta, const double phi);
double complex Ylm_l12m_5(const double theta, const double phi);
double complex Ylm_l12m_6(const double theta, const double phi);
double complex Ylm_l12m_7(const double theta, const double phi);
double complex Ylm_l12m_8(const double theta, const double phi);
double complex Ylm_l12m_9(const double theta, const double phi);
double complex Ylm_l12m_10(const double theta, const double phi);
double complex Ylm_l12m_11(const double theta, const double phi);
double complex Ylm_l12m_12(const double theta, const double phi);
double complex Ylm_l13m_1(const double theta, const double phi);
double complex Ylm_l13m_2(const double theta, const double phi);
double complex Ylm_l13m_3(const double theta, const double phi);
double complex Ylm_l13m_4(const double theta, const double phi);
double complex Ylm_l13m_5(const double theta, const double phi);
double complex Ylm_l13m_6(const double theta, const double phi);
double complex Ylm_l13m_7(const double theta, const double phi);
double complex Ylm_l13m_8(const double theta, const double phi);
double complex Ylm_l13m_9(const double theta, const double phi);
double complex Ylm_l13m_10(const double theta, const double phi);
double complex Ylm_l13m_11(const double theta, const double phi);
double complex Ylm_l13m_12(const double theta, const double phi);
double complex Ylm_l13m_13(const double theta, const double phi);
double complex Ylm_l14m_1(const double theta, const double phi);
double complex Ylm_l14m_2(const double theta, const double phi);
double complex Ylm_l14m_3(const double theta, const double phi);
double complex Ylm_l14m_4(const double theta, const double phi);
double complex Ylm_l14m_5(const double theta, const double phi);
double complex Ylm_l14m_6(const double theta, const double phi);
double complex Ylm_l14m_7(const double theta, const double phi);
double complex Ylm_l14m_8(const double theta, const double phi);
double complex Ylm_l14m_9(const double theta, const double phi);
double complex Ylm_l14m_10(const double theta, const double phi);
double complex Ylm_l14m_11(const double theta, const double phi);
double complex Ylm_l14m_12(const double theta, const double phi);
double complex Ylm_l14m_13(const double theta, const double phi);
double complex Ylm_l14m_14(const double theta, const double phi);



/* Y_n^m(\theta, \varphi) := \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}} exp(i m \varphi)\mathrm{P}_n^m\left(\cos(\theta)\right) 
// ->return value: Y_{l}^{m}(x) */
double complex Ylm(const int l, int m, const double theta, const double phi)
{
  fYlm_T *Y[lmax][lmax];
  fYlm_T *Y_[lmax][lmax];
  if (theta > M_PI || theta < 0)
    abortEr("theta exceeds from [0,pi] interval.\n");
  if (phi > 2*M_PI || phi < 0)
    abortEr("phi exceeds from [0,2*pi] interval.\n");
  if (l >= 15)
    abortEr("l exceeds the maximum.\n");
  if (l < 0)
    abortEr("l is negative.\n");
  if (abs(m) > l)
    return 0;
  Y[0][0] = Ylm_l0m0;
  Y[1][0] = Ylm_l1m0;
  Y[1][1] = Ylm_l1m1;
  Y[2][0] = Ylm_l2m0;
  Y[2][1] = Ylm_l2m1;
  Y[2][2] = Ylm_l2m2;
  Y[3][0] = Ylm_l3m0;
  Y[3][1] = Ylm_l3m1;
  Y[3][2] = Ylm_l3m2;
  Y[3][3] = Ylm_l3m3;
  Y[4][0] = Ylm_l4m0;
  Y[4][1] = Ylm_l4m1;
  Y[4][2] = Ylm_l4m2;
  Y[4][3] = Ylm_l4m3;
  Y[4][4] = Ylm_l4m4;
  Y[5][0] = Ylm_l5m0;
  Y[5][1] = Ylm_l5m1;
  Y[5][2] = Ylm_l5m2;
  Y[5][3] = Ylm_l5m3;
  Y[5][4] = Ylm_l5m4;
  Y[5][5] = Ylm_l5m5;
  Y[6][0] = Ylm_l6m0;
  Y[6][1] = Ylm_l6m1;
  Y[6][2] = Ylm_l6m2;
  Y[6][3] = Ylm_l6m3;
  Y[6][4] = Ylm_l6m4;
  Y[6][5] = Ylm_l6m5;
  Y[6][6] = Ylm_l6m6;
  Y[7][0] = Ylm_l7m0;
  Y[7][1] = Ylm_l7m1;
  Y[7][2] = Ylm_l7m2;
  Y[7][3] = Ylm_l7m3;
  Y[7][4] = Ylm_l7m4;
  Y[7][5] = Ylm_l7m5;
  Y[7][6] = Ylm_l7m6;
  Y[7][7] = Ylm_l7m7;
  Y[8][0] = Ylm_l8m0;
  Y[8][1] = Ylm_l8m1;
  Y[8][2] = Ylm_l8m2;
  Y[8][3] = Ylm_l8m3;
  Y[8][4] = Ylm_l8m4;
  Y[8][5] = Ylm_l8m5;
  Y[8][6] = Ylm_l8m6;
  Y[8][7] = Ylm_l8m7;
  Y[8][8] = Ylm_l8m8;
  Y[9][0] = Ylm_l9m0;
  Y[9][1] = Ylm_l9m1;
  Y[9][2] = Ylm_l9m2;
  Y[9][3] = Ylm_l9m3;
  Y[9][4] = Ylm_l9m4;
  Y[9][5] = Ylm_l9m5;
  Y[9][6] = Ylm_l9m6;
  Y[9][7] = Ylm_l9m7;
  Y[9][8] = Ylm_l9m8;
  Y[9][9] = Ylm_l9m9;
  Y[10][0] = Ylm_l10m0;
  Y[10][1] = Ylm_l10m1;
  Y[10][2] = Ylm_l10m2;
  Y[10][3] = Ylm_l10m3;
  Y[10][4] = Ylm_l10m4;
  Y[10][5] = Ylm_l10m5;
  Y[10][6] = Ylm_l10m6;
  Y[10][7] = Ylm_l10m7;
  Y[10][8] = Ylm_l10m8;
  Y[10][9] = Ylm_l10m9;
  Y[10][10] = Ylm_l10m10;
  Y[11][0] = Ylm_l11m0;
  Y[11][1] = Ylm_l11m1;
  Y[11][2] = Ylm_l11m2;
  Y[11][3] = Ylm_l11m3;
  Y[11][4] = Ylm_l11m4;
  Y[11][5] = Ylm_l11m5;
  Y[11][6] = Ylm_l11m6;
  Y[11][7] = Ylm_l11m7;
  Y[11][8] = Ylm_l11m8;
  Y[11][9] = Ylm_l11m9;
  Y[11][10] = Ylm_l11m10;
  Y[11][11] = Ylm_l11m11;
  Y[12][0] = Ylm_l12m0;
  Y[12][1] = Ylm_l12m1;
  Y[12][2] = Ylm_l12m2;
  Y[12][3] = Ylm_l12m3;
  Y[12][4] = Ylm_l12m4;
  Y[12][5] = Ylm_l12m5;
  Y[12][6] = Ylm_l12m6;
  Y[12][7] = Ylm_l12m7;
  Y[12][8] = Ylm_l12m8;
  Y[12][9] = Ylm_l12m9;
  Y[12][10] = Ylm_l12m10;
  Y[12][11] = Ylm_l12m11;
  Y[12][12] = Ylm_l12m12;
  Y[13][0] = Ylm_l13m0;
  Y[13][1] = Ylm_l13m1;
  Y[13][2] = Ylm_l13m2;
  Y[13][3] = Ylm_l13m3;
  Y[13][4] = Ylm_l13m4;
  Y[13][5] = Ylm_l13m5;
  Y[13][6] = Ylm_l13m6;
  Y[13][7] = Ylm_l13m7;
  Y[13][8] = Ylm_l13m8;
  Y[13][9] = Ylm_l13m9;
  Y[13][10] = Ylm_l13m10;
  Y[13][11] = Ylm_l13m11;
  Y[13][12] = Ylm_l13m12;
  Y[13][13] = Ylm_l13m13;
  Y[14][0] = Ylm_l14m0;
  Y[14][1] = Ylm_l14m1;
  Y[14][2] = Ylm_l14m2;
  Y[14][3] = Ylm_l14m3;
  Y[14][4] = Ylm_l14m4;
  Y[14][5] = Ylm_l14m5;
  Y[14][6] = Ylm_l14m6;
  Y[14][7] = Ylm_l14m7;
  Y[14][8] = Ylm_l14m8;
  Y[14][9] = Ylm_l14m9;
  Y[14][10] = Ylm_l14m10;
  Y[14][11] = Ylm_l14m11;
  Y[14][12] = Ylm_l14m12;
  Y[14][13] = Ylm_l14m13;
  Y[14][14] = Ylm_l14m14;
  Y_[1][1] = Ylm_l1m_1;
  Y_[2][1] = Ylm_l2m_1;
  Y_[2][2] = Ylm_l2m_2;
  Y_[3][1] = Ylm_l3m_1;
  Y_[3][2] = Ylm_l3m_2;
  Y_[3][3] = Ylm_l3m_3;
  Y_[4][1] = Ylm_l4m_1;
  Y_[4][2] = Ylm_l4m_2;
  Y_[4][3] = Ylm_l4m_3;
  Y_[4][4] = Ylm_l4m_4;
  Y_[5][1] = Ylm_l5m_1;
  Y_[5][2] = Ylm_l5m_2;
  Y_[5][3] = Ylm_l5m_3;
  Y_[5][4] = Ylm_l5m_4;
  Y_[5][5] = Ylm_l5m_5;
  Y_[6][1] = Ylm_l6m_1;
  Y_[6][2] = Ylm_l6m_2;
  Y_[6][3] = Ylm_l6m_3;
  Y_[6][4] = Ylm_l6m_4;
  Y_[6][5] = Ylm_l6m_5;
  Y_[6][6] = Ylm_l6m_6;
  Y_[7][1] = Ylm_l7m_1;
  Y_[7][2] = Ylm_l7m_2;
  Y_[7][3] = Ylm_l7m_3;
  Y_[7][4] = Ylm_l7m_4;
  Y_[7][5] = Ylm_l7m_5;
  Y_[7][6] = Ylm_l7m_6;
  Y_[7][7] = Ylm_l7m_7;
  Y_[8][1] = Ylm_l8m_1;
  Y_[8][2] = Ylm_l8m_2;
  Y_[8][3] = Ylm_l8m_3;
  Y_[8][4] = Ylm_l8m_4;
  Y_[8][5] = Ylm_l8m_5;
  Y_[8][6] = Ylm_l8m_6;
  Y_[8][7] = Ylm_l8m_7;
  Y_[8][8] = Ylm_l8m_8;
  Y_[9][1] = Ylm_l9m_1;
  Y_[9][2] = Ylm_l9m_2;
  Y_[9][3] = Ylm_l9m_3;
  Y_[9][4] = Ylm_l9m_4;
  Y_[9][5] = Ylm_l9m_5;
  Y_[9][6] = Ylm_l9m_6;
  Y_[9][7] = Ylm_l9m_7;
  Y_[9][8] = Ylm_l9m_8;
  Y_[9][9] = Ylm_l9m_9;
  Y_[10][1] = Ylm_l10m_1;
  Y_[10][2] = Ylm_l10m_2;
  Y_[10][3] = Ylm_l10m_3;
  Y_[10][4] = Ylm_l10m_4;
  Y_[10][5] = Ylm_l10m_5;
  Y_[10][6] = Ylm_l10m_6;
  Y_[10][7] = Ylm_l10m_7;
  Y_[10][8] = Ylm_l10m_8;
  Y_[10][9] = Ylm_l10m_9;
  Y_[10][10] = Ylm_l10m_10;
  Y_[11][1] = Ylm_l11m_1;
  Y_[11][2] = Ylm_l11m_2;
  Y_[11][3] = Ylm_l11m_3;
  Y_[11][4] = Ylm_l11m_4;
  Y_[11][5] = Ylm_l11m_5;
  Y_[11][6] = Ylm_l11m_6;
  Y_[11][7] = Ylm_l11m_7;
  Y_[11][8] = Ylm_l11m_8;
  Y_[11][9] = Ylm_l11m_9;
  Y_[11][10] = Ylm_l11m_10;
  Y_[11][11] = Ylm_l11m_11;
  Y_[12][1] = Ylm_l12m_1;
  Y_[12][2] = Ylm_l12m_2;
  Y_[12][3] = Ylm_l12m_3;
  Y_[12][4] = Ylm_l12m_4;
  Y_[12][5] = Ylm_l12m_5;
  Y_[12][6] = Ylm_l12m_6;
  Y_[12][7] = Ylm_l12m_7;
  Y_[12][8] = Ylm_l12m_8;
  Y_[12][9] = Ylm_l12m_9;
  Y_[12][10] = Ylm_l12m_10;
  Y_[12][11] = Ylm_l12m_11;
  Y_[12][12] = Ylm_l12m_12;
  Y_[13][1] = Ylm_l13m_1;
  Y_[13][2] = Ylm_l13m_2;
  Y_[13][3] = Ylm_l13m_3;
  Y_[13][4] = Ylm_l13m_4;
  Y_[13][5] = Ylm_l13m_5;
  Y_[13][6] = Ylm_l13m_6;
  Y_[13][7] = Ylm_l13m_7;
  Y_[13][8] = Ylm_l13m_8;
  Y_[13][9] = Ylm_l13m_9;
  Y_[13][10] = Ylm_l13m_10;
  Y_[13][11] = Ylm_l13m_11;
  Y_[13][12] = Ylm_l13m_12;
  Y_[13][13] = Ylm_l13m_13;
  Y_[14][1] = Ylm_l14m_1;
  Y_[14][2] = Ylm_l14m_2;
  Y_[14][3] = Ylm_l14m_3;
  Y_[14][4] = Ylm_l14m_4;
  Y_[14][5] = Ylm_l14m_5;
  Y_[14][6] = Ylm_l14m_6;
  Y_[14][7] = Ylm_l14m_7;
  Y_[14][8] = Ylm_l14m_8;
  Y_[14][9] = Ylm_l14m_9;
  Y_[14][10] = Ylm_l14m_10;
  Y_[14][11] = Ylm_l14m_11;
  Y_[14][12] = Ylm_l14m_12;
  Y_[14][13] = Ylm_l14m_13;
  Y_[14][14] = Ylm_l14m_14;
  if (m < 0)
    return Y_[l][-m](theta,phi);
  return Y[l][m](theta,phi);
}

/* Y_{0}^{0} */
double complex Ylm_l0m0(const double theta, const double phi)
{
  UNUSED(theta);
  UNUSED(phi);
  return (1.0/2.0)/sqrt(M_PI);
}


/* Y_{1}^{0} */
double complex Ylm_l1m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/2.0)*sqrt(3)*cos(theta)/sqrt(M_PI);
}


/* Y_{1}^{1} */
double complex Ylm_l1m1(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(6)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{2}^{0} */
double complex Ylm_l2m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/4.0)*sqrt(5)*(3*pow(cos(theta), 2) - 1)/sqrt(M_PI);
}


/* Y_{2}^{1} */
double complex Ylm_l2m1(const double theta, const double phi)
{
  return -1.0/8.0*sqrt(30)*cexp(I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* Y_{2}^{2} */
double complex Ylm_l2m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{3}^{0} */
double complex Ylm_l3m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/4.0)*sqrt(7)*(5*pow(cos(theta), 2) - 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{3}^{1} */
double complex Ylm_l3m1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*(1 - 5*pow(cos(theta), 2))*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{3}^{2} */
double complex Ylm_l3m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(210)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{3}^{3} */
double complex Ylm_l3m3(const double theta, const double phi)
{
  return -1.0/8.0*sqrt(35)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{4}^{0} */
double complex Ylm_l4m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (3.0/16.0)*(35*pow(cos(theta), 4) - 30*pow(cos(theta), 2) + 3)/sqrt(M_PI);
}


/* Y_{4}^{1} */
double complex Ylm_l4m1(const double theta, const double phi)
{
  return -sqrt(5)*((3.0/32.0)*sin(2*theta) + (21.0/64.0)*sin(4*theta))*cexp(I*phi)/sqrt(M_PI);
}


/* Y_{4}^{2} */
double complex Ylm_l4m2(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(10)*(6 - 7*pow(sin(theta), 2))*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{4}^{3} */
double complex Ylm_l4m3(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(35)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{4}^{4} */
double complex Ylm_l4m4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(70)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{5}^{0} */
double complex Ylm_l5m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/16.0)*sqrt(11)*(63*pow(cos(theta), 4) - 70*pow(cos(theta), 2) + 15)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{1} */
double complex Ylm_l5m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*(-21*pow(cos(theta), 4) + 14*pow(cos(theta), 2) - 1)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{5}^{2} */
double complex Ylm_l5m2(const double theta, const double phi)
{
  return (1.0/16.0)*sqrt(2310)*(2 - 3*pow(sin(theta), 2))*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{3} */
double complex Ylm_l5m3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(385)*(1 - 9*pow(cos(theta), 2))*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{5}^{4} */
double complex Ylm_l5m4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(770)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{5} */
double complex Ylm_l5m5(const double theta, const double phi)
{
  return -3.0/32.0*sqrt(77)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{6}^{0} */
double complex Ylm_l6m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/32.0)*sqrt(13)*(231*pow(cos(theta), 6) - 315*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)/sqrt(M_PI);
}


/* Y_{6}^{1} */
double complex Ylm_l6m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*(-33*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 5)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{2} */
double complex Ylm_l6m2(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(1365)*(33*pow(sin(theta), 4) - 48*pow(sin(theta), 2) + 16)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{6}^{3} */
double complex Ylm_l6m3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*(3 - 11*pow(cos(theta), 2))*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{4} */
double complex Ylm_l6m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(182)*(11*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{6}^{5} */
double complex Ylm_l6m5(const double theta, const double phi)
{
  return -3.0/32.0*sqrt(1001)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{6} */
double complex Ylm_l6m6(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(3003)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{7}^{0} */
double complex Ylm_l7m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/32.0)*sqrt(15)*(429*pow(cos(theta), 6) - 693*pow(cos(theta), 4) + 315*pow(cos(theta), 2) - 35)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{1} */
double complex Ylm_l7m1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*(-429*pow(cos(theta), 6) + 495*pow(cos(theta), 4) - 135*pow(cos(theta), 2) + 5)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{7}^{2} */
double complex Ylm_l7m2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(35)*(143*pow(sin(theta), 4) - 176*pow(sin(theta), 2) + 48)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{3} */
double complex Ylm_l7m3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(70)*(-143*pow(cos(theta), 4) + 66*pow(cos(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{7}^{4} */
double complex Ylm_l7m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(13*pow(cos(theta), 2) - 3)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{5} */
double complex Ylm_l7m5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(770)*(1 - 13*pow(cos(theta), 2))*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{7}^{6} */
double complex Ylm_l7m6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(5005)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{7} */
double complex Ylm_l7m7(const double theta, const double phi)
{
  return -3.0/128.0*sqrt(1430)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{8}^{0} */
double complex Ylm_l8m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/256.0)*sqrt(17)*(6435*pow(cos(theta), 8) - 12012*pow(cos(theta), 6) + 6930*pow(cos(theta), 4) - 1260*pow(cos(theta), 2) + 35)/sqrt(M_PI);
}


/* Y_{8}^{1} */
double complex Ylm_l8m1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*(-715*pow(cos(theta), 6) + 1001*pow(cos(theta), 4) - 385*pow(cos(theta), 2) + 35)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{2} */
double complex Ylm_l8m2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(595)*(-143*pow(sin(theta), 4) + 253*pow(sin(theta), 2) + 143*pow(cos(theta), 6) - 111)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{8}^{3} */
double complex Ylm_l8m3(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(39270)*(-39*pow(cos(theta), 4) + 26*pow(cos(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{4} */
double complex Ylm_l8m4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2618)*(65*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{8}^{5} */
double complex Ylm_l8m5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34034)*(1 - 5*pow(cos(theta), 2))*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{6} */
double complex Ylm_l8m6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(7293)*(15*pow(cos(theta), 2) - 1)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{8}^{7} */
double complex Ylm_l8m7(const double theta, const double phi)
{
  return -3.0/128.0*sqrt(24310)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{8} */
double complex Ylm_l8m8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(24310)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{9}^{0} */
double complex Ylm_l9m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/256.0)*sqrt(19)*(12155*pow(cos(theta), 8) - 25740*pow(cos(theta), 6) + 18018*pow(cos(theta), 4) - 4620*pow(cos(theta), 2) + 315)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{1} */
double complex Ylm_l9m1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*(-2431*pow(cos(theta), 8) + 4004*pow(cos(theta), 6) - 2002*pow(cos(theta), 4) + 308*pow(cos(theta), 2) - 7)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{9}^{2} */
double complex Ylm_l9m2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1045)*(-273*pow(sin(theta), 4) + 455*pow(sin(theta), 2) + 221*pow(cos(theta), 6) - 189)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{3} */
double complex Ylm_l9m3(const double theta, const double phi)
{
  return (1.0/256.0)*sqrt(21945)*(-221*pow(cos(theta), 6) + 195*pow(cos(theta), 4) - 39*pow(cos(theta), 2) + 1)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{9}^{4} */
double complex Ylm_l9m4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(190190)*(17*pow(cos(theta), 4) - 10*pow(cos(theta), 2) + 1)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{5} */
double complex Ylm_l9m5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2717)*(-85*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 1)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{9}^{6} */
double complex Ylm_l9m6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(40755)*(17*pow(cos(theta), 2) - 3)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{7} */
double complex Ylm_l9m7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(13585)*(1 - 17*pow(cos(theta), 2))*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{9}^{8} */
double complex Ylm_l9m8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(461890)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{9} */
double complex Ylm_l9m9(const double theta, const double phi)
{
  return -1.0/512.0*sqrt(230945)*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{10}^{0} */
double complex Ylm_l10m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/512.0)*sqrt(21)*(46189*pow(cos(theta), 10) - 109395*pow(cos(theta), 8) + 90090*pow(cos(theta), 6) - 30030*pow(cos(theta), 4) + 3465*pow(cos(theta), 2) - 63)/sqrt(M_PI);
}


/* Y_{10}^{1} */
double complex Ylm_l10m1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*(-4199*pow(cos(theta), 8) + 7956*pow(cos(theta), 6) - 4914*pow(cos(theta), 4) + 1092*pow(cos(theta), 2) - 63)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{2} */
double complex Ylm_l10m2(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(770)*(2730*pow(sin(theta), 4) - 5096*pow(sin(theta), 2) + 4199*pow(cos(theta), 8) - 6188*pow(cos(theta), 6) + 2373)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{10}^{3} */
double complex Ylm_l10m3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(5005)*(-323*pow(cos(theta), 6) + 357*pow(cos(theta), 4) - 105*pow(cos(theta), 2) + 7)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{4} */
double complex Ylm_l10m4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(10010)*(323*pow(cos(theta), 6) - 255*pow(cos(theta), 4) + 45*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{10}^{5} */
double complex Ylm_l10m5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(1001)*(-323*pow(cos(theta), 4) + 170*pow(cos(theta), 2) - 15)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{6} */
double complex Ylm_l10m6(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(5005)*(323*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 3)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{10}^{7} */
double complex Ylm_l10m7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(85085)*(3 - 19*pow(cos(theta), 2))*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{8} */
double complex Ylm_l10m8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(510510)*(19*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{10}^{9} */
double complex Ylm_l10m9(const double theta, const double phi)
{
  return -1.0/512.0*sqrt(4849845)*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{10} */
double complex Ylm_l10m10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(969969)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{11}^{0} */
double complex Ylm_l11m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/512.0)*sqrt(23)*(88179*pow(cos(theta), 10) - 230945*pow(cos(theta), 8) + 218790*pow(cos(theta), 6) - 90090*pow(cos(theta), 4) + 15015*pow(cos(theta), 2) - 693)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{1} */
double complex Ylm_l11m1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*(-29393*pow(cos(theta), 10) + 62985*pow(cos(theta), 8) - 46410*pow(cos(theta), 6) + 13650*pow(cos(theta), 4) - 1365*pow(cos(theta), 2) + 21)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{11}^{2} */
double complex Ylm_l11m2(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(98670)*(2142*pow(sin(theta), 4) - 3864*pow(sin(theta), 2) + 2261*pow(cos(theta), 8) - 3876*pow(cos(theta), 6) + 1743)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{3} */
double complex Ylm_l11m3(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(345345)*(-969*pow(cos(theta), 8) + 1292*pow(cos(theta), 6) - 510*pow(cos(theta), 4) + 60*pow(cos(theta), 2) - 1)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{11}^{4} */
double complex Ylm_l11m4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(46046)*(323*pow(cos(theta), 6) - 323*pow(cos(theta), 4) + 85*pow(cos(theta), 2) - 5)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{5} */
double complex Ylm_l11m5(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(6578)*(-2261*pow(cos(theta), 6) + 1615*pow(cos(theta), 4) - 255*pow(cos(theta), 2) + 5)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{11}^{6} */
double complex Ylm_l11m6(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(167739)*(399*pow(cos(theta), 4) - 190*pow(cos(theta), 2) + 15)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{7} */
double complex Ylm_l11m7(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(1677390)*(-133*pow(cos(theta), 4) + 38*pow(cos(theta), 2) - 1)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{11}^{8} */
double complex Ylm_l11m8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(31870410)*(7*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{9} */
double complex Ylm_l11m9(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(2124694)*(1 - 21*pow(cos(theta), 2))*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{11}^{10} */
double complex Ylm_l11m10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(22309287)*cexp(10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{11} */
double complex Ylm_l11m11(const double theta, const double phi)
{
  return -1.0/2048.0*sqrt(4056234)*cexp(11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{12}^{0} */
double complex Ylm_l12m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (5.0/2048.0)*(676039*pow(cos(theta), 12) - 1939938*pow(cos(theta), 10) + 2078505*pow(cos(theta), 8) - 1021020*pow(cos(theta), 6) + 225225*pow(cos(theta), 4) - 18018*pow(cos(theta), 2) + 231)/sqrt(M_PI);
}


/* Y_{12}^{1} */
double complex Ylm_l12m1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*(-52003*pow(cos(theta), 10) + 124355*pow(cos(theta), 8) - 106590*pow(cos(theta), 6) + 39270*pow(cos(theta), 4) - 5775*pow(cos(theta), 2) + 231)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{2} */
double complex Ylm_l12m2(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(6006)*(-2550*pow(sin(theta), 4) + 4875*pow(sin(theta), 2) + 7429*pow(cos(theta), 10) - 14535*pow(cos(theta), 8) + 9690*pow(cos(theta), 6) - 2328)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{12}^{3} */
double complex Ylm_l12m3(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(1001)*(-7429*pow(cos(theta), 8) + 11628*pow(cos(theta), 6) - 5814*pow(cos(theta), 4) + 1020*pow(cos(theta), 2) - 45)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{4} */
double complex Ylm_l12m4(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(1001)*(7429*pow(cos(theta), 8) - 9044*pow(cos(theta), 6) + 3230*pow(cos(theta), 4) - 340*pow(cos(theta), 2) + 5)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{12}^{5} */
double complex Ylm_l12m5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(34034)*(-437*pow(cos(theta), 6) + 399*pow(cos(theta), 4) - 95*pow(cos(theta), 2) + 5)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{6} */
double complex Ylm_l12m6(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(2431)*(3059*pow(cos(theta), 6) - 1995*pow(cos(theta), 4) + 285*pow(cos(theta), 2) - 5)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{12}^{7} */
double complex Ylm_l12m7(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(277134)*(-161*pow(cos(theta), 4) + 70*pow(cos(theta), 2) - 5)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{8} */
double complex Ylm_l12m8(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(277134)*(161*pow(cos(theta), 4) - 42*pow(cos(theta), 2) + 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{12}^{9} */
double complex Ylm_l12m9(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(646646)*(3 - 23*pow(cos(theta), 2))*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{10} */
double complex Ylm_l12m10(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(88179)*(23*pow(cos(theta), 2) - 1)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{12}^{11} */
double complex Ylm_l12m11(const double theta, const double phi)
{
  return -5.0/2048.0*sqrt(4056234)*cexp(11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{12} */
double complex Ylm_l12m12(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(676039)*cexp(12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{13}^{0} */
double complex Ylm_l13m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (3.0/2048.0)*sqrt(3)*(1300075*pow(cos(theta), 12) - 4056234*pow(cos(theta), 10) + 4849845*pow(cos(theta), 8) - 2771340*pow(cos(theta), 6) + 765765*pow(cos(theta), 4) - 90090*pow(cos(theta), 2) + 3003)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{1} */
double complex Ylm_l13m1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*(-185725*pow(cos(theta), 12) + 490314*pow(cos(theta), 10) - 479655*pow(cos(theta), 8) + 213180*pow(cos(theta), 6) - 42075*pow(cos(theta), 4) + 2970*pow(cos(theta), 2) - 33)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{13}^{2} */
double complex Ylm_l13m2(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2730)*(-21318*pow(sin(theta), 4) + 39831*pow(sin(theta), 2) + 37145*pow(cos(theta), 10) - 81719*pow(cos(theta), 8) + 63954*pow(cos(theta), 6) - 18612)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{3} */
double complex Ylm_l13m3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(30030)*(-37145*pow(cos(theta), 10) + 66861*pow(cos(theta), 8) - 40698*pow(cos(theta), 6) + 9690*pow(cos(theta), 4) - 765*pow(cos(theta), 2) + 9)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{13}^{4} */
double complex Ylm_l13m4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(51051)*(10925*pow(cos(theta), 8) - 15732*pow(cos(theta), 6) + 7182*pow(cos(theta), 4) - 1140*pow(cos(theta), 2) + 45)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{5} */
double complex Ylm_l13m5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(102102)*(-10925*pow(cos(theta), 8) + 12236*pow(cos(theta), 6) - 3990*pow(cos(theta), 4) + 380*pow(cos(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{13}^{6} */
double complex Ylm_l13m6(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(969969)*(575*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{7} */
double complex Ylm_l13m7(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(692835)*(-805*pow(cos(theta), 6) + 483*pow(cos(theta), 4) - 63*pow(cos(theta), 2) + 1)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{13}^{8} */
double complex Ylm_l13m8(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(9699690)*(115*pow(cos(theta), 4) - 46*pow(cos(theta), 2) + 3)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{9} */
double complex Ylm_l13m9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(88179)*(-575*pow(cos(theta), 4) + 138*pow(cos(theta), 2) - 3)*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{13}^{10} */
double complex Ylm_l13m10(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2028117)*(25*pow(cos(theta), 2) - 3)*cexp(10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{11} */
double complex Ylm_l13m11(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(4056234)*(1 - 25*pow(cos(theta), 2))*cexp(11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{13}^{12} */
double complex Ylm_l13m12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(2028117)*cexp(12*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{13} */
double complex Ylm_l13m13(const double theta, const double phi)
{
  return -15.0/8192.0*sqrt(312018)*cexp(13*I*phi)*pow(sin(theta), 13)/sqrt(M_PI);
}


/* Y_{14}^{0} */
double complex Ylm_l14m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/4096.0)*sqrt(29)*(5014575*pow(cos(theta), 14) - 16900975*pow(cos(theta), 12) + 22309287*pow(cos(theta), 10) - 14549535*pow(cos(theta), 8) + 4849845*pow(cos(theta), 6) - 765765*pow(cos(theta), 4) + 45045*pow(cos(theta), 2) - 429)/sqrt(M_PI);
}


/* Y_{14}^{1} */
double complex Ylm_l14m1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*(-334305*pow(cos(theta), 12) + 965770*pow(cos(theta), 10) - 1062347*pow(cos(theta), 8) + 554268*pow(cos(theta), 6) - 138567*pow(cos(theta), 4) + 14586*pow(cos(theta), 2) - 429)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{2} */
double complex Ylm_l14m2(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(79170)*(-334305*pow(cos(theta), 14) + 1151495*pow(cos(theta), 12) - 1552661*pow(cos(theta), 10) + 1033923*pow(cos(theta), 8) - 351747*pow(cos(theta), 6) + 56661*pow(cos(theta), 4) - 3399*pow(cos(theta), 2) + 33)*cexp(2*I*phi)/sqrt(M_PI);
}


/* Y_{14}^{3} */
double complex Ylm_l14m3(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(448630)*(-58995*pow(cos(theta), 10) + 120175*pow(cos(theta), 8) - 86526*pow(cos(theta), 6) + 26334*pow(cos(theta), 4) - 3135*pow(cos(theta), 2) + 99)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{4} */
double complex Ylm_l14m4(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(2467465)*(6555*pow(cos(theta), 10) - 10925*pow(cos(theta), 8) + 6118*pow(cos(theta), 6) - 1330*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{14}^{5} */
double complex Ylm_l14m5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(18752734)*(-1725*pow(cos(theta), 8) + 2300*pow(cos(theta), 6) - 966*pow(cos(theta), 4) + 140*pow(cos(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{6} */
double complex Ylm_l14m6(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(93763670)*(3105*pow(cos(theta), 8) - 3220*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 84*pow(cos(theta), 2) + 1)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{14}^{7} */
double complex Ylm_l14m7(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(20092215)*(-1035*pow(cos(theta), 6) + 805*pow(cos(theta), 4) - 161*pow(cos(theta), 2) + 7)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{8} */
double complex Ylm_l14m8(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(25571910)*(1035*pow(cos(theta), 6) - 575*pow(cos(theta), 4) + 69*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{14}^{9} */
double complex Ylm_l14m9(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(98025655)*(-135*pow(cos(theta), 4) + 50*pow(cos(theta), 2) - 3)*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{10} */
double complex Ylm_l14m10(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(117630786)*(225*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 1)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{14}^{11} */
double complex Ylm_l14m11(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(1 - 9*pow(cos(theta), 2))*cexp(11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{12} */
double complex Ylm_l14m12(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(1508087)*(27*pow(cos(theta), 2) - 1)*cexp(12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{14}^{13} */
double complex Ylm_l14m13(const double theta, const double phi)
{
  return -15.0/8192.0*sqrt(9048522)*cexp(13*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{14} */
double complex Ylm_l14m14(const double theta, const double phi)
{
  return (15.0/16384.0)*sqrt(1292646)*cexp(14*I*phi)*pow(sin(theta), 14)/sqrt(M_PI);
}


/* Y_{1}^{-1} */
double complex Ylm_l1m_1(const double theta, const double phi)
{
  return (1.0/4.0)*sqrt(6)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{2}^{-1} */
double complex Ylm_l2m_1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(-I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* Y_{2}^{-2} */
double complex Ylm_l2m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{3}^{-1} */
double complex Ylm_l3m_1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*(5*pow(cos(theta), 2) - 1)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{3}^{-2} */
double complex Ylm_l3m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(210)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{3}^{-3} */
double complex Ylm_l3m_3(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(35)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{4}^{-1} */
double complex Ylm_l4m_1(const double theta, const double phi)
{
  return sqrt(5)*((3.0/32.0)*sin(2*theta) + (21.0/64.0)*sin(4*theta))*cexp(-I*phi)/sqrt(M_PI);
}


/* Y_{4}^{-2} */
double complex Ylm_l4m_2(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(10)*(6 - 7*pow(sin(theta), 2))*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{4}^{-3} */
double complex Ylm_l4m_3(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(35)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{4}^{-4} */
double complex Ylm_l4m_4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(70)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{5}^{-1} */
double complex Ylm_l5m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*(21*pow(cos(theta), 4) - 14*pow(cos(theta), 2) + 1)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{5}^{-2} */
double complex Ylm_l5m_2(const double theta, const double phi)
{
  return (1.0/16.0)*sqrt(2310)*(2 - 3*pow(sin(theta), 2))*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{-3} */
double complex Ylm_l5m_3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(385)*(9*pow(cos(theta), 2) - 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{5}^{-4} */
double complex Ylm_l5m_4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(770)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{-5} */
double complex Ylm_l5m_5(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(77)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{6}^{-1} */
double complex Ylm_l6m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*(33*pow(cos(theta), 4) - 30*pow(cos(theta), 2) + 5)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{-2} */
double complex Ylm_l6m_2(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(1365)*(33*pow(sin(theta), 4) - 48*pow(sin(theta), 2) + 16)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{6}^{-3} */
double complex Ylm_l6m_3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*(11*pow(cos(theta), 2) - 3)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{-4} */
double complex Ylm_l6m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(182)*(11*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{6}^{-5} */
double complex Ylm_l6m_5(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1001)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{-6} */
double complex Ylm_l6m_6(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(3003)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{7}^{-1} */
double complex Ylm_l7m_1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*(429*pow(cos(theta), 6) - 495*pow(cos(theta), 4) + 135*pow(cos(theta), 2) - 5)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{7}^{-2} */
double complex Ylm_l7m_2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(35)*(143*pow(sin(theta), 4) - 176*pow(sin(theta), 2) + 48)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{-3} */
double complex Ylm_l7m_3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(70)*(143*pow(cos(theta), 4) - 66*pow(cos(theta), 2) + 3)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{7}^{-4} */
double complex Ylm_l7m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(13*pow(cos(theta), 2) - 3)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{-5} */
double complex Ylm_l7m_5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(770)*(13*pow(cos(theta), 2) - 1)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{7}^{-6} */
double complex Ylm_l7m_6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(5005)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{-7} */
double complex Ylm_l7m_7(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1430)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{8}^{-1} */
double complex Ylm_l8m_1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*(715*pow(cos(theta), 6) - 1001*pow(cos(theta), 4) + 385*pow(cos(theta), 2) - 35)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-2} */
double complex Ylm_l8m_2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(595)*(-143*pow(sin(theta), 4) + 253*pow(sin(theta), 2) + 143*pow(cos(theta), 6) - 111)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{8}^{-3} */
double complex Ylm_l8m_3(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(39270)*(39*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 3)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-4} */
double complex Ylm_l8m_4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2618)*(65*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{8}^{-5} */
double complex Ylm_l8m_5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34034)*(5*pow(cos(theta), 2) - 1)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-6} */
double complex Ylm_l8m_6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(7293)*(15*pow(cos(theta), 2) - 1)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{8}^{-7} */
double complex Ylm_l8m_7(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(24310)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-8} */
double complex Ylm_l8m_8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(24310)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{9}^{-1} */
double complex Ylm_l9m_1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*(2431*pow(cos(theta), 8) - 4004*pow(cos(theta), 6) + 2002*pow(cos(theta), 4) - 308*pow(cos(theta), 2) + 7)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{9}^{-2} */
double complex Ylm_l9m_2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1045)*(-273*pow(sin(theta), 4) + 455*pow(sin(theta), 2) + 221*pow(cos(theta), 6) - 189)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-3} */
double complex Ylm_l9m_3(const double theta, const double phi)
{
  return (1.0/256.0)*sqrt(21945)*(221*pow(cos(theta), 6) - 195*pow(cos(theta), 4) + 39*pow(cos(theta), 2) - 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{9}^{-4} */
double complex Ylm_l9m_4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(190190)*(17*pow(cos(theta), 4) - 10*pow(cos(theta), 2) + 1)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-5} */
double complex Ylm_l9m_5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2717)*(85*pow(cos(theta), 4) - 30*pow(cos(theta), 2) + 1)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{9}^{-6} */
double complex Ylm_l9m_6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(40755)*(17*pow(cos(theta), 2) - 3)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-7} */
double complex Ylm_l9m_7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(13585)*(17*pow(cos(theta), 2) - 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{9}^{-8} */
double complex Ylm_l9m_8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(461890)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-9} */
double complex Ylm_l9m_9(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(230945)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{10}^{-1} */
double complex Ylm_l10m_1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*(4199*pow(cos(theta), 8) - 7956*pow(cos(theta), 6) + 4914*pow(cos(theta), 4) - 1092*pow(cos(theta), 2) + 63)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-2} */
double complex Ylm_l10m_2(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(770)*(2730*pow(sin(theta), 4) - 5096*pow(sin(theta), 2) + 4199*pow(cos(theta), 8) - 6188*pow(cos(theta), 6) + 2373)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{10}^{-3} */
double complex Ylm_l10m_3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(5005)*(323*pow(cos(theta), 6) - 357*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 7)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-4} */
double complex Ylm_l10m_4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(10010)*(323*pow(cos(theta), 6) - 255*pow(cos(theta), 4) + 45*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{10}^{-5} */
double complex Ylm_l10m_5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(1001)*(323*pow(cos(theta), 4) - 170*pow(cos(theta), 2) + 15)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-6} */
double complex Ylm_l10m_6(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(5005)*(323*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 3)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{10}^{-7} */
double complex Ylm_l10m_7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(85085)*(19*pow(cos(theta), 2) - 3)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-8} */
double complex Ylm_l10m_8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(510510)*(19*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{10}^{-9} */
double complex Ylm_l10m_9(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(4849845)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-10} */
double complex Ylm_l10m_10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(969969)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{11}^{-1} */
double complex Ylm_l11m_1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*(29393*pow(cos(theta), 10) - 62985*pow(cos(theta), 8) + 46410*pow(cos(theta), 6) - 13650*pow(cos(theta), 4) + 1365*pow(cos(theta), 2) - 21)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{11}^{-2} */
double complex Ylm_l11m_2(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(98670)*(2142*pow(sin(theta), 4) - 3864*pow(sin(theta), 2) + 2261*pow(cos(theta), 8) - 3876*pow(cos(theta), 6) + 1743)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-3} */
double complex Ylm_l11m_3(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(345345)*(969*pow(cos(theta), 8) - 1292*pow(cos(theta), 6) + 510*pow(cos(theta), 4) - 60*pow(cos(theta), 2) + 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{11}^{-4} */
double complex Ylm_l11m_4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(46046)*(323*pow(cos(theta), 6) - 323*pow(cos(theta), 4) + 85*pow(cos(theta), 2) - 5)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-5} */
double complex Ylm_l11m_5(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(6578)*(2261*pow(cos(theta), 6) - 1615*pow(cos(theta), 4) + 255*pow(cos(theta), 2) - 5)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{11}^{-6} */
double complex Ylm_l11m_6(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(167739)*(399*pow(cos(theta), 4) - 190*pow(cos(theta), 2) + 15)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-7} */
double complex Ylm_l11m_7(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(1677390)*(133*pow(cos(theta), 4) - 38*pow(cos(theta), 2) + 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{11}^{-8} */
double complex Ylm_l11m_8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(31870410)*(7*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-9} */
double complex Ylm_l11m_9(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(2124694)*(21*pow(cos(theta), 2) - 1)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{11}^{-10} */
double complex Ylm_l11m_10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(22309287)*cexp(-10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-11} */
double complex Ylm_l11m_11(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(4056234)*cexp(-11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{12}^{-1} */
double complex Ylm_l12m_1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*(52003*pow(cos(theta), 10) - 124355*pow(cos(theta), 8) + 106590*pow(cos(theta), 6) - 39270*pow(cos(theta), 4) + 5775*pow(cos(theta), 2) - 231)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-2} */
double complex Ylm_l12m_2(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(6006)*(-2550*pow(sin(theta), 4) + 4875*pow(sin(theta), 2) + 7429*pow(cos(theta), 10) - 14535*pow(cos(theta), 8) + 9690*pow(cos(theta), 6) - 2328)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{12}^{-3} */
double complex Ylm_l12m_3(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(1001)*(7429*pow(cos(theta), 8) - 11628*pow(cos(theta), 6) + 5814*pow(cos(theta), 4) - 1020*pow(cos(theta), 2) + 45)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-4} */
double complex Ylm_l12m_4(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(1001)*(7429*pow(cos(theta), 8) - 9044*pow(cos(theta), 6) + 3230*pow(cos(theta), 4) - 340*pow(cos(theta), 2) + 5)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{12}^{-5} */
double complex Ylm_l12m_5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(34034)*(437*pow(cos(theta), 6) - 399*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 5)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-6} */
double complex Ylm_l12m_6(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(2431)*(3059*pow(cos(theta), 6) - 1995*pow(cos(theta), 4) + 285*pow(cos(theta), 2) - 5)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{12}^{-7} */
double complex Ylm_l12m_7(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(277134)*(161*pow(cos(theta), 4) - 70*pow(cos(theta), 2) + 5)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-8} */
double complex Ylm_l12m_8(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(277134)*(161*pow(cos(theta), 4) - 42*pow(cos(theta), 2) + 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{12}^{-9} */
double complex Ylm_l12m_9(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(646646)*(23*pow(cos(theta), 2) - 3)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-10} */
double complex Ylm_l12m_10(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(88179)*(23*pow(cos(theta), 2) - 1)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{12}^{-11} */
double complex Ylm_l12m_11(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(4056234)*cexp(-11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-12} */
double complex Ylm_l12m_12(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(676039)*cexp(-12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{13}^{-1} */
double complex Ylm_l13m_1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*(185725*pow(cos(theta), 12) - 490314*pow(cos(theta), 10) + 479655*pow(cos(theta), 8) - 213180*pow(cos(theta), 6) + 42075*pow(cos(theta), 4) - 2970*pow(cos(theta), 2) + 33)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{13}^{-2} */
double complex Ylm_l13m_2(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2730)*(-21318*pow(sin(theta), 4) + 39831*pow(sin(theta), 2) + 37145*pow(cos(theta), 10) - 81719*pow(cos(theta), 8) + 63954*pow(cos(theta), 6) - 18612)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-3} */
double complex Ylm_l13m_3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(30030)*(37145*pow(cos(theta), 10) - 66861*pow(cos(theta), 8) + 40698*pow(cos(theta), 6) - 9690*pow(cos(theta), 4) + 765*pow(cos(theta), 2) - 9)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{13}^{-4} */
double complex Ylm_l13m_4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(51051)*(10925*pow(cos(theta), 8) - 15732*pow(cos(theta), 6) + 7182*pow(cos(theta), 4) - 1140*pow(cos(theta), 2) + 45)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-5} */
double complex Ylm_l13m_5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(102102)*(10925*pow(cos(theta), 8) - 12236*pow(cos(theta), 6) + 3990*pow(cos(theta), 4) - 380*pow(cos(theta), 2) + 5)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{13}^{-6} */
double complex Ylm_l13m_6(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(969969)*(575*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-7} */
double complex Ylm_l13m_7(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(692835)*(805*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 63*pow(cos(theta), 2) - 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{13}^{-8} */
double complex Ylm_l13m_8(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(9699690)*(115*pow(cos(theta), 4) - 46*pow(cos(theta), 2) + 3)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-9} */
double complex Ylm_l13m_9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(88179)*(575*pow(cos(theta), 4) - 138*pow(cos(theta), 2) + 3)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{13}^{-10} */
double complex Ylm_l13m_10(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2028117)*(25*pow(cos(theta), 2) - 3)*cexp(-10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-11} */
double complex Ylm_l13m_11(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(4056234)*(25*pow(cos(theta), 2) - 1)*cexp(-11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{13}^{-12} */
double complex Ylm_l13m_12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(2028117)*cexp(-12*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-13} */
double complex Ylm_l13m_13(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(312018)*cexp(-13*I*phi)*pow(sin(theta), 13)/sqrt(M_PI);
}


/* Y_{14}^{-1} */
double complex Ylm_l14m_1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*(334305*pow(cos(theta), 12) - 965770*pow(cos(theta), 10) + 1062347*pow(cos(theta), 8) - 554268*pow(cos(theta), 6) + 138567*pow(cos(theta), 4) - 14586*pow(cos(theta), 2) + 429)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-2} */
double complex Ylm_l14m_2(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(79170)*(-334305*pow(cos(theta), 14) + 1151495*pow(cos(theta), 12) - 1552661*pow(cos(theta), 10) + 1033923*pow(cos(theta), 8) - 351747*pow(cos(theta), 6) + 56661*pow(cos(theta), 4) - 3399*pow(cos(theta), 2) + 33)*cexp(-2*I*phi)/sqrt(M_PI);
}


/* Y_{14}^{-3} */
double complex Ylm_l14m_3(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(448630)*(58995*pow(cos(theta), 10) - 120175*pow(cos(theta), 8) + 86526*pow(cos(theta), 6) - 26334*pow(cos(theta), 4) + 3135*pow(cos(theta), 2) - 99)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-4} */
double complex Ylm_l14m_4(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(2467465)*(6555*pow(cos(theta), 10) - 10925*pow(cos(theta), 8) + 6118*pow(cos(theta), 6) - 1330*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{14}^{-5} */
double complex Ylm_l14m_5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(18752734)*(1725*pow(cos(theta), 8) - 2300*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 140*pow(cos(theta), 2) + 5)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-6} */
double complex Ylm_l14m_6(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(93763670)*(3105*pow(cos(theta), 8) - 3220*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 84*pow(cos(theta), 2) + 1)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{14}^{-7} */
double complex Ylm_l14m_7(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(20092215)*(1035*pow(cos(theta), 6) - 805*pow(cos(theta), 4) + 161*pow(cos(theta), 2) - 7)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-8} */
double complex Ylm_l14m_8(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(25571910)*(1035*pow(cos(theta), 6) - 575*pow(cos(theta), 4) + 69*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{14}^{-9} */
double complex Ylm_l14m_9(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(98025655)*(135*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 3)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-10} */
double complex Ylm_l14m_10(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(117630786)*(225*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 1)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{14}^{-11} */
double complex Ylm_l14m_11(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(9*pow(cos(theta), 2) - 1)*cexp(-11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-12} */
double complex Ylm_l14m_12(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(1508087)*(27*pow(cos(theta), 2) - 1)*cexp(-12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{14}^{-13} */
double complex Ylm_l14m_13(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(9048522)*cexp(-13*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-14} */
double complex Ylm_l14m_14(const double theta, const double phi)
{
  return (15.0/16384.0)*sqrt(1292646)*cexp(-14*I*phi)*pow(sin(theta), 14)/sqrt(M_PI);
}

