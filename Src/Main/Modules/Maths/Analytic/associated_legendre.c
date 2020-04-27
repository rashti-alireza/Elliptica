 /* direct output of this file to associated_legendre.c */
#include "core_lib.h"
#include "maths_analytic_lib.h"

#define LMAX_ARRAY_SIZE_PLM 15
typedef double fPlm_T (const double x);
fPlm_T *P[LMAX_ARRAY_SIZE_PLM][LMAX_ARRAY_SIZE_PLM];
fPlm_T *P_[LMAX_ARRAY_SIZE_PLM][LMAX_ARRAY_SIZE_PLM];
void init_associated_legendre(void);
double associated_legendre(const int l, const int m, const double x);
static double associated_legendre_P_l0m0(const double x);
static double associated_legendre_P_l1m0(const double x);
static double associated_legendre_P_l1m1(const double x);
static double associated_legendre_P_l2m0(const double x);
static double associated_legendre_P_l2m1(const double x);
static double associated_legendre_P_l2m2(const double x);
static double associated_legendre_P_l3m0(const double x);
static double associated_legendre_P_l3m1(const double x);
static double associated_legendre_P_l3m2(const double x);
static double associated_legendre_P_l3m3(const double x);
static double associated_legendre_P_l4m0(const double x);
static double associated_legendre_P_l4m1(const double x);
static double associated_legendre_P_l4m2(const double x);
static double associated_legendre_P_l4m3(const double x);
static double associated_legendre_P_l4m4(const double x);
static double associated_legendre_P_l5m0(const double x);
static double associated_legendre_P_l5m1(const double x);
static double associated_legendre_P_l5m2(const double x);
static double associated_legendre_P_l5m3(const double x);
static double associated_legendre_P_l5m4(const double x);
static double associated_legendre_P_l5m5(const double x);
static double associated_legendre_P_l6m0(const double x);
static double associated_legendre_P_l6m1(const double x);
static double associated_legendre_P_l6m2(const double x);
static double associated_legendre_P_l6m3(const double x);
static double associated_legendre_P_l6m4(const double x);
static double associated_legendre_P_l6m5(const double x);
static double associated_legendre_P_l6m6(const double x);
static double associated_legendre_P_l7m0(const double x);
static double associated_legendre_P_l7m1(const double x);
static double associated_legendre_P_l7m2(const double x);
static double associated_legendre_P_l7m3(const double x);
static double associated_legendre_P_l7m4(const double x);
static double associated_legendre_P_l7m5(const double x);
static double associated_legendre_P_l7m6(const double x);
static double associated_legendre_P_l7m7(const double x);
static double associated_legendre_P_l8m0(const double x);
static double associated_legendre_P_l8m1(const double x);
static double associated_legendre_P_l8m2(const double x);
static double associated_legendre_P_l8m3(const double x);
static double associated_legendre_P_l8m4(const double x);
static double associated_legendre_P_l8m5(const double x);
static double associated_legendre_P_l8m6(const double x);
static double associated_legendre_P_l8m7(const double x);
static double associated_legendre_P_l8m8(const double x);
static double associated_legendre_P_l9m0(const double x);
static double associated_legendre_P_l9m1(const double x);
static double associated_legendre_P_l9m2(const double x);
static double associated_legendre_P_l9m3(const double x);
static double associated_legendre_P_l9m4(const double x);
static double associated_legendre_P_l9m5(const double x);
static double associated_legendre_P_l9m6(const double x);
static double associated_legendre_P_l9m7(const double x);
static double associated_legendre_P_l9m8(const double x);
static double associated_legendre_P_l9m9(const double x);
static double associated_legendre_P_l10m0(const double x);
static double associated_legendre_P_l10m1(const double x);
static double associated_legendre_P_l10m2(const double x);
static double associated_legendre_P_l10m3(const double x);
static double associated_legendre_P_l10m4(const double x);
static double associated_legendre_P_l10m5(const double x);
static double associated_legendre_P_l10m6(const double x);
static double associated_legendre_P_l10m7(const double x);
static double associated_legendre_P_l10m8(const double x);
static double associated_legendre_P_l10m9(const double x);
static double associated_legendre_P_l10m10(const double x);
static double associated_legendre_P_l11m0(const double x);
static double associated_legendre_P_l11m1(const double x);
static double associated_legendre_P_l11m2(const double x);
static double associated_legendre_P_l11m3(const double x);
static double associated_legendre_P_l11m4(const double x);
static double associated_legendre_P_l11m5(const double x);
static double associated_legendre_P_l11m6(const double x);
static double associated_legendre_P_l11m7(const double x);
static double associated_legendre_P_l11m8(const double x);
static double associated_legendre_P_l11m9(const double x);
static double associated_legendre_P_l11m10(const double x);
static double associated_legendre_P_l11m11(const double x);
static double associated_legendre_P_l12m0(const double x);
static double associated_legendre_P_l12m1(const double x);
static double associated_legendre_P_l12m2(const double x);
static double associated_legendre_P_l12m3(const double x);
static double associated_legendre_P_l12m4(const double x);
static double associated_legendre_P_l12m5(const double x);
static double associated_legendre_P_l12m6(const double x);
static double associated_legendre_P_l12m7(const double x);
static double associated_legendre_P_l12m8(const double x);
static double associated_legendre_P_l12m9(const double x);
static double associated_legendre_P_l12m10(const double x);
static double associated_legendre_P_l12m11(const double x);
static double associated_legendre_P_l12m12(const double x);
static double associated_legendre_P_l13m0(const double x);
static double associated_legendre_P_l13m1(const double x);
static double associated_legendre_P_l13m2(const double x);
static double associated_legendre_P_l13m3(const double x);
static double associated_legendre_P_l13m4(const double x);
static double associated_legendre_P_l13m5(const double x);
static double associated_legendre_P_l13m6(const double x);
static double associated_legendre_P_l13m7(const double x);
static double associated_legendre_P_l13m8(const double x);
static double associated_legendre_P_l13m9(const double x);
static double associated_legendre_P_l13m10(const double x);
static double associated_legendre_P_l13m11(const double x);
static double associated_legendre_P_l13m12(const double x);
static double associated_legendre_P_l13m13(const double x);
static double associated_legendre_P_l14m0(const double x);
static double associated_legendre_P_l14m1(const double x);
static double associated_legendre_P_l14m2(const double x);
static double associated_legendre_P_l14m3(const double x);
static double associated_legendre_P_l14m4(const double x);
static double associated_legendre_P_l14m5(const double x);
static double associated_legendre_P_l14m6(const double x);
static double associated_legendre_P_l14m7(const double x);
static double associated_legendre_P_l14m8(const double x);
static double associated_legendre_P_l14m9(const double x);
static double associated_legendre_P_l14m10(const double x);
static double associated_legendre_P_l14m11(const double x);
static double associated_legendre_P_l14m12(const double x);
static double associated_legendre_P_l14m13(const double x);
static double associated_legendre_P_l14m14(const double x);
static double associated_legendre_P_l1m_1(const double x);
static double associated_legendre_P_l2m_1(const double x);
static double associated_legendre_P_l2m_2(const double x);
static double associated_legendre_P_l3m_1(const double x);
static double associated_legendre_P_l3m_2(const double x);
static double associated_legendre_P_l3m_3(const double x);
static double associated_legendre_P_l4m_1(const double x);
static double associated_legendre_P_l4m_2(const double x);
static double associated_legendre_P_l4m_3(const double x);
static double associated_legendre_P_l4m_4(const double x);
static double associated_legendre_P_l5m_1(const double x);
static double associated_legendre_P_l5m_2(const double x);
static double associated_legendre_P_l5m_3(const double x);
static double associated_legendre_P_l5m_4(const double x);
static double associated_legendre_P_l5m_5(const double x);
static double associated_legendre_P_l6m_1(const double x);
static double associated_legendre_P_l6m_2(const double x);
static double associated_legendre_P_l6m_3(const double x);
static double associated_legendre_P_l6m_4(const double x);
static double associated_legendre_P_l6m_5(const double x);
static double associated_legendre_P_l6m_6(const double x);
static double associated_legendre_P_l7m_1(const double x);
static double associated_legendre_P_l7m_2(const double x);
static double associated_legendre_P_l7m_3(const double x);
static double associated_legendre_P_l7m_4(const double x);
static double associated_legendre_P_l7m_5(const double x);
static double associated_legendre_P_l7m_6(const double x);
static double associated_legendre_P_l7m_7(const double x);
static double associated_legendre_P_l8m_1(const double x);
static double associated_legendre_P_l8m_2(const double x);
static double associated_legendre_P_l8m_3(const double x);
static double associated_legendre_P_l8m_4(const double x);
static double associated_legendre_P_l8m_5(const double x);
static double associated_legendre_P_l8m_6(const double x);
static double associated_legendre_P_l8m_7(const double x);
static double associated_legendre_P_l8m_8(const double x);
static double associated_legendre_P_l9m_1(const double x);
static double associated_legendre_P_l9m_2(const double x);
static double associated_legendre_P_l9m_3(const double x);
static double associated_legendre_P_l9m_4(const double x);
static double associated_legendre_P_l9m_5(const double x);
static double associated_legendre_P_l9m_6(const double x);
static double associated_legendre_P_l9m_7(const double x);
static double associated_legendre_P_l9m_8(const double x);
static double associated_legendre_P_l9m_9(const double x);
static double associated_legendre_P_l10m_1(const double x);
static double associated_legendre_P_l10m_2(const double x);
static double associated_legendre_P_l10m_3(const double x);
static double associated_legendre_P_l10m_4(const double x);
static double associated_legendre_P_l10m_5(const double x);
static double associated_legendre_P_l10m_6(const double x);
static double associated_legendre_P_l10m_7(const double x);
static double associated_legendre_P_l10m_8(const double x);
static double associated_legendre_P_l10m_9(const double x);
static double associated_legendre_P_l10m_10(const double x);
static double associated_legendre_P_l11m_1(const double x);
static double associated_legendre_P_l11m_2(const double x);
static double associated_legendre_P_l11m_3(const double x);
static double associated_legendre_P_l11m_4(const double x);
static double associated_legendre_P_l11m_5(const double x);
static double associated_legendre_P_l11m_6(const double x);
static double associated_legendre_P_l11m_7(const double x);
static double associated_legendre_P_l11m_8(const double x);
static double associated_legendre_P_l11m_9(const double x);
static double associated_legendre_P_l11m_10(const double x);
static double associated_legendre_P_l11m_11(const double x);
static double associated_legendre_P_l12m_1(const double x);
static double associated_legendre_P_l12m_2(const double x);
static double associated_legendre_P_l12m_3(const double x);
static double associated_legendre_P_l12m_4(const double x);
static double associated_legendre_P_l12m_5(const double x);
static double associated_legendre_P_l12m_6(const double x);
static double associated_legendre_P_l12m_7(const double x);
static double associated_legendre_P_l12m_8(const double x);
static double associated_legendre_P_l12m_9(const double x);
static double associated_legendre_P_l12m_10(const double x);
static double associated_legendre_P_l12m_11(const double x);
static double associated_legendre_P_l12m_12(const double x);
static double associated_legendre_P_l13m_1(const double x);
static double associated_legendre_P_l13m_2(const double x);
static double associated_legendre_P_l13m_3(const double x);
static double associated_legendre_P_l13m_4(const double x);
static double associated_legendre_P_l13m_5(const double x);
static double associated_legendre_P_l13m_6(const double x);
static double associated_legendre_P_l13m_7(const double x);
static double associated_legendre_P_l13m_8(const double x);
static double associated_legendre_P_l13m_9(const double x);
static double associated_legendre_P_l13m_10(const double x);
static double associated_legendre_P_l13m_11(const double x);
static double associated_legendre_P_l13m_12(const double x);
static double associated_legendre_P_l13m_13(const double x);
static double associated_legendre_P_l14m_1(const double x);
static double associated_legendre_P_l14m_2(const double x);
static double associated_legendre_P_l14m_3(const double x);
static double associated_legendre_P_l14m_4(const double x);
static double associated_legendre_P_l14m_5(const double x);
static double associated_legendre_P_l14m_6(const double x);
static double associated_legendre_P_l14m_7(const double x);
static double associated_legendre_P_l14m_8(const double x);
static double associated_legendre_P_l14m_9(const double x);
static double associated_legendre_P_l14m_10(const double x);
static double associated_legendre_P_l14m_11(const double x);
static double associated_legendre_P_l14m_12(const double x);
static double associated_legendre_P_l14m_13(const double x);
static double associated_legendre_P_l14m_14(const double x);



/* initializing P_{l}^{m}(x) table */
void init_associated_legendre(void)
{
  P[0][0] = associated_legendre_P_l0m0;
  P[1][0] = associated_legendre_P_l1m0;
  P[1][1] = associated_legendre_P_l1m1;
  P[2][0] = associated_legendre_P_l2m0;
  P[2][1] = associated_legendre_P_l2m1;
  P[2][2] = associated_legendre_P_l2m2;
  P[3][0] = associated_legendre_P_l3m0;
  P[3][1] = associated_legendre_P_l3m1;
  P[3][2] = associated_legendre_P_l3m2;
  P[3][3] = associated_legendre_P_l3m3;
  P[4][0] = associated_legendre_P_l4m0;
  P[4][1] = associated_legendre_P_l4m1;
  P[4][2] = associated_legendre_P_l4m2;
  P[4][3] = associated_legendre_P_l4m3;
  P[4][4] = associated_legendre_P_l4m4;
  P[5][0] = associated_legendre_P_l5m0;
  P[5][1] = associated_legendre_P_l5m1;
  P[5][2] = associated_legendre_P_l5m2;
  P[5][3] = associated_legendre_P_l5m3;
  P[5][4] = associated_legendre_P_l5m4;
  P[5][5] = associated_legendre_P_l5m5;
  P[6][0] = associated_legendre_P_l6m0;
  P[6][1] = associated_legendre_P_l6m1;
  P[6][2] = associated_legendre_P_l6m2;
  P[6][3] = associated_legendre_P_l6m3;
  P[6][4] = associated_legendre_P_l6m4;
  P[6][5] = associated_legendre_P_l6m5;
  P[6][6] = associated_legendre_P_l6m6;
  P[7][0] = associated_legendre_P_l7m0;
  P[7][1] = associated_legendre_P_l7m1;
  P[7][2] = associated_legendre_P_l7m2;
  P[7][3] = associated_legendre_P_l7m3;
  P[7][4] = associated_legendre_P_l7m4;
  P[7][5] = associated_legendre_P_l7m5;
  P[7][6] = associated_legendre_P_l7m6;
  P[7][7] = associated_legendre_P_l7m7;
  P[8][0] = associated_legendre_P_l8m0;
  P[8][1] = associated_legendre_P_l8m1;
  P[8][2] = associated_legendre_P_l8m2;
  P[8][3] = associated_legendre_P_l8m3;
  P[8][4] = associated_legendre_P_l8m4;
  P[8][5] = associated_legendre_P_l8m5;
  P[8][6] = associated_legendre_P_l8m6;
  P[8][7] = associated_legendre_P_l8m7;
  P[8][8] = associated_legendre_P_l8m8;
  P[9][0] = associated_legendre_P_l9m0;
  P[9][1] = associated_legendre_P_l9m1;
  P[9][2] = associated_legendre_P_l9m2;
  P[9][3] = associated_legendre_P_l9m3;
  P[9][4] = associated_legendre_P_l9m4;
  P[9][5] = associated_legendre_P_l9m5;
  P[9][6] = associated_legendre_P_l9m6;
  P[9][7] = associated_legendre_P_l9m7;
  P[9][8] = associated_legendre_P_l9m8;
  P[9][9] = associated_legendre_P_l9m9;
  P[10][0] = associated_legendre_P_l10m0;
  P[10][1] = associated_legendre_P_l10m1;
  P[10][2] = associated_legendre_P_l10m2;
  P[10][3] = associated_legendre_P_l10m3;
  P[10][4] = associated_legendre_P_l10m4;
  P[10][5] = associated_legendre_P_l10m5;
  P[10][6] = associated_legendre_P_l10m6;
  P[10][7] = associated_legendre_P_l10m7;
  P[10][8] = associated_legendre_P_l10m8;
  P[10][9] = associated_legendre_P_l10m9;
  P[10][10] = associated_legendre_P_l10m10;
  P[11][0] = associated_legendre_P_l11m0;
  P[11][1] = associated_legendre_P_l11m1;
  P[11][2] = associated_legendre_P_l11m2;
  P[11][3] = associated_legendre_P_l11m3;
  P[11][4] = associated_legendre_P_l11m4;
  P[11][5] = associated_legendre_P_l11m5;
  P[11][6] = associated_legendre_P_l11m6;
  P[11][7] = associated_legendre_P_l11m7;
  P[11][8] = associated_legendre_P_l11m8;
  P[11][9] = associated_legendre_P_l11m9;
  P[11][10] = associated_legendre_P_l11m10;
  P[11][11] = associated_legendre_P_l11m11;
  P[12][0] = associated_legendre_P_l12m0;
  P[12][1] = associated_legendre_P_l12m1;
  P[12][2] = associated_legendre_P_l12m2;
  P[12][3] = associated_legendre_P_l12m3;
  P[12][4] = associated_legendre_P_l12m4;
  P[12][5] = associated_legendre_P_l12m5;
  P[12][6] = associated_legendre_P_l12m6;
  P[12][7] = associated_legendre_P_l12m7;
  P[12][8] = associated_legendre_P_l12m8;
  P[12][9] = associated_legendre_P_l12m9;
  P[12][10] = associated_legendre_P_l12m10;
  P[12][11] = associated_legendre_P_l12m11;
  P[12][12] = associated_legendre_P_l12m12;
  P[13][0] = associated_legendre_P_l13m0;
  P[13][1] = associated_legendre_P_l13m1;
  P[13][2] = associated_legendre_P_l13m2;
  P[13][3] = associated_legendre_P_l13m3;
  P[13][4] = associated_legendre_P_l13m4;
  P[13][5] = associated_legendre_P_l13m5;
  P[13][6] = associated_legendre_P_l13m6;
  P[13][7] = associated_legendre_P_l13m7;
  P[13][8] = associated_legendre_P_l13m8;
  P[13][9] = associated_legendre_P_l13m9;
  P[13][10] = associated_legendre_P_l13m10;
  P[13][11] = associated_legendre_P_l13m11;
  P[13][12] = associated_legendre_P_l13m12;
  P[13][13] = associated_legendre_P_l13m13;
  P[14][0] = associated_legendre_P_l14m0;
  P[14][1] = associated_legendre_P_l14m1;
  P[14][2] = associated_legendre_P_l14m2;
  P[14][3] = associated_legendre_P_l14m3;
  P[14][4] = associated_legendre_P_l14m4;
  P[14][5] = associated_legendre_P_l14m5;
  P[14][6] = associated_legendre_P_l14m6;
  P[14][7] = associated_legendre_P_l14m7;
  P[14][8] = associated_legendre_P_l14m8;
  P[14][9] = associated_legendre_P_l14m9;
  P[14][10] = associated_legendre_P_l14m10;
  P[14][11] = associated_legendre_P_l14m11;
  P[14][12] = associated_legendre_P_l14m12;
  P[14][13] = associated_legendre_P_l14m13;
  P[14][14] = associated_legendre_P_l14m14;
  P_[1][1] = associated_legendre_P_l1m_1;
  P_[2][1] = associated_legendre_P_l2m_1;
  P_[2][2] = associated_legendre_P_l2m_2;
  P_[3][1] = associated_legendre_P_l3m_1;
  P_[3][2] = associated_legendre_P_l3m_2;
  P_[3][3] = associated_legendre_P_l3m_3;
  P_[4][1] = associated_legendre_P_l4m_1;
  P_[4][2] = associated_legendre_P_l4m_2;
  P_[4][3] = associated_legendre_P_l4m_3;
  P_[4][4] = associated_legendre_P_l4m_4;
  P_[5][1] = associated_legendre_P_l5m_1;
  P_[5][2] = associated_legendre_P_l5m_2;
  P_[5][3] = associated_legendre_P_l5m_3;
  P_[5][4] = associated_legendre_P_l5m_4;
  P_[5][5] = associated_legendre_P_l5m_5;
  P_[6][1] = associated_legendre_P_l6m_1;
  P_[6][2] = associated_legendre_P_l6m_2;
  P_[6][3] = associated_legendre_P_l6m_3;
  P_[6][4] = associated_legendre_P_l6m_4;
  P_[6][5] = associated_legendre_P_l6m_5;
  P_[6][6] = associated_legendre_P_l6m_6;
  P_[7][1] = associated_legendre_P_l7m_1;
  P_[7][2] = associated_legendre_P_l7m_2;
  P_[7][3] = associated_legendre_P_l7m_3;
  P_[7][4] = associated_legendre_P_l7m_4;
  P_[7][5] = associated_legendre_P_l7m_5;
  P_[7][6] = associated_legendre_P_l7m_6;
  P_[7][7] = associated_legendre_P_l7m_7;
  P_[8][1] = associated_legendre_P_l8m_1;
  P_[8][2] = associated_legendre_P_l8m_2;
  P_[8][3] = associated_legendre_P_l8m_3;
  P_[8][4] = associated_legendre_P_l8m_4;
  P_[8][5] = associated_legendre_P_l8m_5;
  P_[8][6] = associated_legendre_P_l8m_6;
  P_[8][7] = associated_legendre_P_l8m_7;
  P_[8][8] = associated_legendre_P_l8m_8;
  P_[9][1] = associated_legendre_P_l9m_1;
  P_[9][2] = associated_legendre_P_l9m_2;
  P_[9][3] = associated_legendre_P_l9m_3;
  P_[9][4] = associated_legendre_P_l9m_4;
  P_[9][5] = associated_legendre_P_l9m_5;
  P_[9][6] = associated_legendre_P_l9m_6;
  P_[9][7] = associated_legendre_P_l9m_7;
  P_[9][8] = associated_legendre_P_l9m_8;
  P_[9][9] = associated_legendre_P_l9m_9;
  P_[10][1] = associated_legendre_P_l10m_1;
  P_[10][2] = associated_legendre_P_l10m_2;
  P_[10][3] = associated_legendre_P_l10m_3;
  P_[10][4] = associated_legendre_P_l10m_4;
  P_[10][5] = associated_legendre_P_l10m_5;
  P_[10][6] = associated_legendre_P_l10m_6;
  P_[10][7] = associated_legendre_P_l10m_7;
  P_[10][8] = associated_legendre_P_l10m_8;
  P_[10][9] = associated_legendre_P_l10m_9;
  P_[10][10] = associated_legendre_P_l10m_10;
  P_[11][1] = associated_legendre_P_l11m_1;
  P_[11][2] = associated_legendre_P_l11m_2;
  P_[11][3] = associated_legendre_P_l11m_3;
  P_[11][4] = associated_legendre_P_l11m_4;
  P_[11][5] = associated_legendre_P_l11m_5;
  P_[11][6] = associated_legendre_P_l11m_6;
  P_[11][7] = associated_legendre_P_l11m_7;
  P_[11][8] = associated_legendre_P_l11m_8;
  P_[11][9] = associated_legendre_P_l11m_9;
  P_[11][10] = associated_legendre_P_l11m_10;
  P_[11][11] = associated_legendre_P_l11m_11;
  P_[12][1] = associated_legendre_P_l12m_1;
  P_[12][2] = associated_legendre_P_l12m_2;
  P_[12][3] = associated_legendre_P_l12m_3;
  P_[12][4] = associated_legendre_P_l12m_4;
  P_[12][5] = associated_legendre_P_l12m_5;
  P_[12][6] = associated_legendre_P_l12m_6;
  P_[12][7] = associated_legendre_P_l12m_7;
  P_[12][8] = associated_legendre_P_l12m_8;
  P_[12][9] = associated_legendre_P_l12m_9;
  P_[12][10] = associated_legendre_P_l12m_10;
  P_[12][11] = associated_legendre_P_l12m_11;
  P_[12][12] = associated_legendre_P_l12m_12;
  P_[13][1] = associated_legendre_P_l13m_1;
  P_[13][2] = associated_legendre_P_l13m_2;
  P_[13][3] = associated_legendre_P_l13m_3;
  P_[13][4] = associated_legendre_P_l13m_4;
  P_[13][5] = associated_legendre_P_l13m_5;
  P_[13][6] = associated_legendre_P_l13m_6;
  P_[13][7] = associated_legendre_P_l13m_7;
  P_[13][8] = associated_legendre_P_l13m_8;
  P_[13][9] = associated_legendre_P_l13m_9;
  P_[13][10] = associated_legendre_P_l13m_10;
  P_[13][11] = associated_legendre_P_l13m_11;
  P_[13][12] = associated_legendre_P_l13m_12;
  P_[13][13] = associated_legendre_P_l13m_13;
  P_[14][1] = associated_legendre_P_l14m_1;
  P_[14][2] = associated_legendre_P_l14m_2;
  P_[14][3] = associated_legendre_P_l14m_3;
  P_[14][4] = associated_legendre_P_l14m_4;
  P_[14][5] = associated_legendre_P_l14m_5;
  P_[14][6] = associated_legendre_P_l14m_6;
  P_[14][7] = associated_legendre_P_l14m_7;
  P_[14][8] = associated_legendre_P_l14m_8;
  P_[14][9] = associated_legendre_P_l14m_9;
  P_[14][10] = associated_legendre_P_l14m_10;
  P_[14][11] = associated_legendre_P_l14m_11;
  P_[14][12] = associated_legendre_P_l14m_12;
  P_[14][13] = associated_legendre_P_l14m_13;
  P_[14][14] = associated_legendre_P_l14m_14;
}
/* P_{l}^{m}(x)=\left( -1\right) ^{m}\left( 1-x^{2}\right) ^{\frac {m} {2}}\frac {d^{m}P_{l}\left( x\right) } {dx^{m}}
// ->return value: P_{l}^{m}(x) */
double associated_legendre(const int l, const int m, const double x)
{
  if (x > 1. || x < -1.)
    Error0("x exceeds from [-1,1] interval.\n");
  if (l >= 15)
    Error0("l exceeds the maximum.\n"
            "To go higher number change lmax in 'associated_legendre.py'.\n");
  if (l < 0)
    Error0("l is negative.\n");
  if (abs(m) > l)
    return 0;
  if (m < 0)
    return P_[l][-m](x);
  return P[l][m](x);
}

/* P_{0}^{0} */
static double associated_legendre_P_l0m0(const double x)
{
  UNUSED(x);
  return 1;
}


/* P_{1}^{0} */
static double associated_legendre_P_l1m0(const double x)
{
  return x;
}


/* P_{1}^{1} */
static double associated_legendre_P_l1m1(const double x)
{
  return -sqrt(1 - pow(x, 2));
}


/* P_{2}^{0} */
static double associated_legendre_P_l2m0(const double x)
{
  return (3.0/2.0)*pow(x, 2) - 1.0/2.0;
}


/* P_{2}^{1} */
static double associated_legendre_P_l2m1(const double x)
{
  return -3*x*sqrt(1 - pow(x, 2));
}


/* P_{2}^{2} */
static double associated_legendre_P_l2m2(const double x)
{
  return 3 - 3*pow(x, 2);
}


/* P_{3}^{0} */
static double associated_legendre_P_l3m0(const double x)
{
  return (5.0/2.0)*pow(x, 3) - 3.0/2.0*x;
}


/* P_{3}^{1} */
static double associated_legendre_P_l3m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((15.0/2.0)*pow(x, 2) - 3.0/2.0);
}


/* P_{3}^{2} */
static double associated_legendre_P_l3m2(const double x)
{
  return 15*x*(1 - pow(x, 2));
}


/* P_{3}^{3} */
static double associated_legendre_P_l3m3(const double x)
{
  return -15*pow(1 - pow(x, 2), 3.0/2.0);
}


/* P_{4}^{0} */
static double associated_legendre_P_l4m0(const double x)
{
  return (35.0/8.0)*pow(x, 4) - 15.0/4.0*pow(x, 2) + 3.0/8.0;
}


/* P_{4}^{1} */
static double associated_legendre_P_l4m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((35.0/2.0)*pow(x, 3) - 15.0/2.0*x);
}


/* P_{4}^{2} */
static double associated_legendre_P_l4m2(const double x)
{
  return (1 - pow(x, 2))*((105.0/2.0)*pow(x, 2) - 15.0/2.0);
}


/* P_{4}^{3} */
static double associated_legendre_P_l4m3(const double x)
{
  return -105*x*pow(1 - pow(x, 2), 3.0/2.0);
}


/* P_{4}^{4} */
static double associated_legendre_P_l4m4(const double x)
{
  return 105*pow(1 - pow(x, 2), 2);
}


/* P_{5}^{0} */
static double associated_legendre_P_l5m0(const double x)
{
  return (63.0/8.0)*pow(x, 5) - 35.0/4.0*pow(x, 3) + (15.0/8.0)*x;
}


/* P_{5}^{1} */
static double associated_legendre_P_l5m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((315.0/8.0)*pow(x, 4) - 105.0/4.0*pow(x, 2) + 15.0/8.0);
}


/* P_{5}^{2} */
static double associated_legendre_P_l5m2(const double x)
{
  return (1 - pow(x, 2))*((315.0/2.0)*pow(x, 3) - 105.0/2.0*x);
}


/* P_{5}^{3} */
static double associated_legendre_P_l5m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((945.0/2.0)*pow(x, 2) - 105.0/2.0);
}


/* P_{5}^{4} */
static double associated_legendre_P_l5m4(const double x)
{
  return 945*x*pow(1 - pow(x, 2), 2);
}


/* P_{5}^{5} */
static double associated_legendre_P_l5m5(const double x)
{
  return -945*pow(1 - pow(x, 2), 5.0/2.0);
}


/* P_{6}^{0} */
static double associated_legendre_P_l6m0(const double x)
{
  return (231.0/16.0)*pow(x, 6) - 315.0/16.0*pow(x, 4) + (105.0/16.0)*pow(x, 2) - 5.0/16.0;
}


/* P_{6}^{1} */
static double associated_legendre_P_l6m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((693.0/8.0)*pow(x, 5) - 315.0/4.0*pow(x, 3) + (105.0/8.0)*x);
}


/* P_{6}^{2} */
static double associated_legendre_P_l6m2(const double x)
{
  return (1 - pow(x, 2))*((3465.0/8.0)*pow(x, 4) - 945.0/4.0*pow(x, 2) + 105.0/8.0);
}


/* P_{6}^{3} */
static double associated_legendre_P_l6m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((3465.0/2.0)*pow(x, 3) - 945.0/2.0*x);
}


/* P_{6}^{4} */
static double associated_legendre_P_l6m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((10395.0/2.0)*pow(x, 2) - 945.0/2.0);
}


/* P_{6}^{5} */
static double associated_legendre_P_l6m5(const double x)
{
  return -10395*x*pow(1 - pow(x, 2), 5.0/2.0);
}


/* P_{6}^{6} */
static double associated_legendre_P_l6m6(const double x)
{
  return 10395*pow(1 - pow(x, 2), 3);
}


/* P_{7}^{0} */
static double associated_legendre_P_l7m0(const double x)
{
  return (429.0/16.0)*pow(x, 7) - 693.0/16.0*pow(x, 5) + (315.0/16.0)*pow(x, 3) - 35.0/16.0*x;
}


/* P_{7}^{1} */
static double associated_legendre_P_l7m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((3003.0/16.0)*pow(x, 6) - 3465.0/16.0*pow(x, 4) + (945.0/16.0)*pow(x, 2) - 35.0/16.0);
}


/* P_{7}^{2} */
static double associated_legendre_P_l7m2(const double x)
{
  return (1 - pow(x, 2))*((9009.0/8.0)*pow(x, 5) - 3465.0/4.0*pow(x, 3) + (945.0/8.0)*x);
}


/* P_{7}^{3} */
static double associated_legendre_P_l7m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((45045.0/8.0)*pow(x, 4) - 10395.0/4.0*pow(x, 2) + 945.0/8.0);
}


/* P_{7}^{4} */
static double associated_legendre_P_l7m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((45045.0/2.0)*pow(x, 3) - 10395.0/2.0*x);
}


/* P_{7}^{5} */
static double associated_legendre_P_l7m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((135135.0/2.0)*pow(x, 2) - 10395.0/2.0);
}


/* P_{7}^{6} */
static double associated_legendre_P_l7m6(const double x)
{
  return 135135*x*pow(1 - pow(x, 2), 3);
}


/* P_{7}^{7} */
static double associated_legendre_P_l7m7(const double x)
{
  return -135135*pow(1 - pow(x, 2), 7.0/2.0);
}


/* P_{8}^{0} */
static double associated_legendre_P_l8m0(const double x)
{
  return (6435.0/128.0)*pow(x, 8) - 3003.0/32.0*pow(x, 6) + (3465.0/64.0)*pow(x, 4) - 315.0/32.0*pow(x, 2) + 35.0/128.0;
}


/* P_{8}^{1} */
static double associated_legendre_P_l8m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((6435.0/16.0)*pow(x, 7) - 9009.0/16.0*pow(x, 5) + (3465.0/16.0)*pow(x, 3) - 315.0/16.0*x);
}


/* P_{8}^{2} */
static double associated_legendre_P_l8m2(const double x)
{
  return (1 - pow(x, 2))*((45045.0/16.0)*pow(x, 6) - 45045.0/16.0*pow(x, 4) + (10395.0/16.0)*pow(x, 2) - 315.0/16.0);
}


/* P_{8}^{3} */
static double associated_legendre_P_l8m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((135135.0/8.0)*pow(x, 5) - 45045.0/4.0*pow(x, 3) + (10395.0/8.0)*x);
}


/* P_{8}^{4} */
static double associated_legendre_P_l8m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((675675.0/8.0)*pow(x, 4) - 135135.0/4.0*pow(x, 2) + 10395.0/8.0);
}


/* P_{8}^{5} */
static double associated_legendre_P_l8m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((675675.0/2.0)*pow(x, 3) - 135135.0/2.0*x);
}


/* P_{8}^{6} */
static double associated_legendre_P_l8m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((2027025.0/2.0)*pow(x, 2) - 135135.0/2.0);
}


/* P_{8}^{7} */
static double associated_legendre_P_l8m7(const double x)
{
  return -2027025*x*pow(1 - pow(x, 2), 7.0/2.0);
}


/* P_{8}^{8} */
static double associated_legendre_P_l8m8(const double x)
{
  return 2027025*pow(1 - pow(x, 2), 4);
}


/* P_{9}^{0} */
static double associated_legendre_P_l9m0(const double x)
{
  return (12155.0/128.0)*pow(x, 9) - 6435.0/32.0*pow(x, 7) + (9009.0/64.0)*pow(x, 5) - 1155.0/32.0*pow(x, 3) + (315.0/128.0)*x;
}


/* P_{9}^{1} */
static double associated_legendre_P_l9m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((109395.0/128.0)*pow(x, 8) - 45045.0/32.0*pow(x, 6) + (45045.0/64.0)*pow(x, 4) - 3465.0/32.0*pow(x, 2) + 315.0/128.0);
}


/* P_{9}^{2} */
static double associated_legendre_P_l9m2(const double x)
{
  return (1 - pow(x, 2))*((109395.0/16.0)*pow(x, 7) - 135135.0/16.0*pow(x, 5) + (45045.0/16.0)*pow(x, 3) - 3465.0/16.0*x);
}


/* P_{9}^{3} */
static double associated_legendre_P_l9m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((765765.0/16.0)*pow(x, 6) - 675675.0/16.0*pow(x, 4) + (135135.0/16.0)*pow(x, 2) - 3465.0/16.0);
}


/* P_{9}^{4} */
static double associated_legendre_P_l9m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((2297295.0/8.0)*pow(x, 5) - 675675.0/4.0*pow(x, 3) + (135135.0/8.0)*x);
}


/* P_{9}^{5} */
static double associated_legendre_P_l9m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((11486475.0/8.0)*pow(x, 4) - 2027025.0/4.0*pow(x, 2) + 135135.0/8.0);
}


/* P_{9}^{6} */
static double associated_legendre_P_l9m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((11486475.0/2.0)*pow(x, 3) - 2027025.0/2.0*x);
}


/* P_{9}^{7} */
static double associated_legendre_P_l9m7(const double x)
{
  return -pow(1 - pow(x, 2), 7.0/2.0)*((34459425.0/2.0)*pow(x, 2) - 2027025.0/2.0);
}


/* P_{9}^{8} */
static double associated_legendre_P_l9m8(const double x)
{
  return 34459425*x*pow(1 - pow(x, 2), 4);
}


/* P_{9}^{9} */
static double associated_legendre_P_l9m9(const double x)
{
  return -34459425*pow(1 - pow(x, 2), 9.0/2.0);
}


/* P_{10}^{0} */
static double associated_legendre_P_l10m0(const double x)
{
  return (46189.0/256.0)*pow(x, 10) - 109395.0/256.0*pow(x, 8) + (45045.0/128.0)*pow(x, 6) - 15015.0/128.0*pow(x, 4) + (3465.0/256.0)*pow(x, 2) - 63.0/256.0;
}


/* P_{10}^{1} */
static double associated_legendre_P_l10m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((230945.0/128.0)*pow(x, 9) - 109395.0/32.0*pow(x, 7) + (135135.0/64.0)*pow(x, 5) - 15015.0/32.0*pow(x, 3) + (3465.0/128.0)*x);
}


/* P_{10}^{2} */
static double associated_legendre_P_l10m2(const double x)
{
  return (1 - pow(x, 2))*((2078505.0/128.0)*pow(x, 8) - 765765.0/32.0*pow(x, 6) + (675675.0/64.0)*pow(x, 4) - 45045.0/32.0*pow(x, 2) + 3465.0/128.0);
}


/* P_{10}^{3} */
static double associated_legendre_P_l10m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((2078505.0/16.0)*pow(x, 7) - 2297295.0/16.0*pow(x, 5) + (675675.0/16.0)*pow(x, 3) - 45045.0/16.0*x);
}


/* P_{10}^{4} */
static double associated_legendre_P_l10m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((14549535.0/16.0)*pow(x, 6) - 11486475.0/16.0*pow(x, 4) + (2027025.0/16.0)*pow(x, 2) - 45045.0/16.0);
}


/* P_{10}^{5} */
static double associated_legendre_P_l10m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((43648605.0/8.0)*pow(x, 5) - 11486475.0/4.0*pow(x, 3) + (2027025.0/8.0)*x);
}


/* P_{10}^{6} */
static double associated_legendre_P_l10m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((218243025.0/8.0)*pow(x, 4) - 34459425.0/4.0*pow(x, 2) + 2027025.0/8.0);
}


/* P_{10}^{7} */
static double associated_legendre_P_l10m7(const double x)
{
  return -pow(1 - pow(x, 2), 7.0/2.0)*((218243025.0/2.0)*pow(x, 3) - 34459425.0/2.0*x);
}


/* P_{10}^{8} */
static double associated_legendre_P_l10m8(const double x)
{
  return pow(1 - pow(x, 2), 4)*((654729075.0/2.0)*pow(x, 2) - 34459425.0/2.0);
}


/* P_{10}^{9} */
static double associated_legendre_P_l10m9(const double x)
{
  return -654729075*x*pow(1 - pow(x, 2), 9.0/2.0);
}


/* P_{10}^{10} */
static double associated_legendre_P_l10m10(const double x)
{
  return 654729075*pow(1 - pow(x, 2), 5);
}


/* P_{11}^{0} */
static double associated_legendre_P_l11m0(const double x)
{
  return (88179.0/256.0)*pow(x, 11) - 230945.0/256.0*pow(x, 9) + (109395.0/128.0)*pow(x, 7) - 45045.0/128.0*pow(x, 5) + (15015.0/256.0)*pow(x, 3) - 693.0/256.0*x;
}


/* P_{11}^{1} */
static double associated_legendre_P_l11m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((969969.0/256.0)*pow(x, 10) - 2078505.0/256.0*pow(x, 8) + (765765.0/128.0)*pow(x, 6) - 225225.0/128.0*pow(x, 4) + (45045.0/256.0)*pow(x, 2) - 693.0/256.0);
}


/* P_{11}^{2} */
static double associated_legendre_P_l11m2(const double x)
{
  return (1 - pow(x, 2))*((4849845.0/128.0)*pow(x, 9) - 2078505.0/32.0*pow(x, 7) + (2297295.0/64.0)*pow(x, 5) - 225225.0/32.0*pow(x, 3) + (45045.0/128.0)*x);
}


/* P_{11}^{3} */
static double associated_legendre_P_l11m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((43648605.0/128.0)*pow(x, 8) - 14549535.0/32.0*pow(x, 6) + (11486475.0/64.0)*pow(x, 4) - 675675.0/32.0*pow(x, 2) + 45045.0/128.0);
}


/* P_{11}^{4} */
static double associated_legendre_P_l11m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((43648605.0/16.0)*pow(x, 7) - 43648605.0/16.0*pow(x, 5) + (11486475.0/16.0)*pow(x, 3) - 675675.0/16.0*x);
}


/* P_{11}^{5} */
static double associated_legendre_P_l11m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((305540235.0/16.0)*pow(x, 6) - 218243025.0/16.0*pow(x, 4) + (34459425.0/16.0)*pow(x, 2) - 675675.0/16.0);
}


/* P_{11}^{6} */
static double associated_legendre_P_l11m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((916620705.0/8.0)*pow(x, 5) - 218243025.0/4.0*pow(x, 3) + (34459425.0/8.0)*x);
}


/* P_{11}^{7} */
static double associated_legendre_P_l11m7(const double x)
{
  return -pow(1 - pow(x, 2), 7.0/2.0)*((4583103525.0/8.0)*pow(x, 4) - 654729075.0/4.0*pow(x, 2) + 34459425.0/8.0);
}


/* P_{11}^{8} */
static double associated_legendre_P_l11m8(const double x)
{
  return pow(1 - pow(x, 2), 4)*((4583103525.0/2.0)*pow(x, 3) - 654729075.0/2.0*x);
}


/* P_{11}^{9} */
static double associated_legendre_P_l11m9(const double x)
{
  return -pow(1 - pow(x, 2), 9.0/2.0)*((13749310575.0/2.0)*pow(x, 2) - 654729075.0/2.0);
}


/* P_{11}^{10} */
static double associated_legendre_P_l11m10(const double x)
{
  return 13749310575*x*pow(1 - pow(x, 2), 5);
}


/* P_{11}^{11} */
static double associated_legendre_P_l11m11(const double x)
{
  return -13749310575*pow(1 - pow(x, 2), 11.0/2.0);
}


/* P_{12}^{0} */
static double associated_legendre_P_l12m0(const double x)
{
  return (676039.0/1024.0)*pow(x, 12) - 969969.0/512.0*pow(x, 10) + (2078505.0/1024.0)*pow(x, 8) - 255255.0/256.0*pow(x, 6) + (225225.0/1024.0)*pow(x, 4) - 9009.0/512.0*pow(x, 2) + 231.0/1024.0;
}


/* P_{12}^{1} */
static double associated_legendre_P_l12m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((2028117.0/256.0)*pow(x, 11) - 4849845.0/256.0*pow(x, 9) + (2078505.0/128.0)*pow(x, 7) - 765765.0/128.0*pow(x, 5) + (225225.0/256.0)*pow(x, 3) - 9009.0/256.0*x);
}


/* P_{12}^{2} */
static double associated_legendre_P_l12m2(const double x)
{
  return (1 - pow(x, 2))*((22309287.0/256.0)*pow(x, 10) - 43648605.0/256.0*pow(x, 8) + (14549535.0/128.0)*pow(x, 6) - 3828825.0/128.0*pow(x, 4) + (675675.0/256.0)*pow(x, 2) - 9009.0/256.0);
}


/* P_{12}^{3} */
static double associated_legendre_P_l12m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((111546435.0/128.0)*pow(x, 9) - 43648605.0/32.0*pow(x, 7) + (43648605.0/64.0)*pow(x, 5) - 3828825.0/32.0*pow(x, 3) + (675675.0/128.0)*x);
}


/* P_{12}^{4} */
static double associated_legendre_P_l12m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((1003917915.0/128.0)*pow(x, 8) - 305540235.0/32.0*pow(x, 6) + (218243025.0/64.0)*pow(x, 4) - 11486475.0/32.0*pow(x, 2) + 675675.0/128.0);
}


/* P_{12}^{5} */
static double associated_legendre_P_l12m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((1003917915.0/16.0)*pow(x, 7) - 916620705.0/16.0*pow(x, 5) + (218243025.0/16.0)*pow(x, 3) - 11486475.0/16.0*x);
}


/* P_{12}^{6} */
static double associated_legendre_P_l12m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((7027425405.0/16.0)*pow(x, 6) - 4583103525.0/16.0*pow(x, 4) + (654729075.0/16.0)*pow(x, 2) - 11486475.0/16.0);
}


/* P_{12}^{7} */
static double associated_legendre_P_l12m7(const double x)
{
  return -pow(1 - pow(x, 2), 7.0/2.0)*((21082276215.0/8.0)*pow(x, 5) - 4583103525.0/4.0*pow(x, 3) + (654729075.0/8.0)*x);
}


/* P_{12}^{8} */
static double associated_legendre_P_l12m8(const double x)
{
  return pow(1 - pow(x, 2), 4)*((105411381075.0/8.0)*pow(x, 4) - 13749310575.0/4.0*pow(x, 2) + 654729075.0/8.0);
}


/* P_{12}^{9} */
static double associated_legendre_P_l12m9(const double x)
{
  return -pow(1 - pow(x, 2), 9.0/2.0)*((105411381075.0/2.0)*pow(x, 3) - 13749310575.0/2.0*x);
}


/* P_{12}^{10} */
static double associated_legendre_P_l12m10(const double x)
{
  return pow(1 - pow(x, 2), 5)*((316234143225.0/2.0)*pow(x, 2) - 13749310575.0/2.0);
}


/* P_{12}^{11} */
static double associated_legendre_P_l12m11(const double x)
{
  return -316234143225*x*pow(1 - pow(x, 2), 11.0/2.0);
}


/* P_{12}^{12} */
static double associated_legendre_P_l12m12(const double x)
{
  return 316234143225*pow(1 - pow(x, 2), 6);
}


/* P_{13}^{0} */
static double associated_legendre_P_l13m0(const double x)
{
  return (1300075.0/1024.0)*pow(x, 13) - 2028117.0/512.0*pow(x, 11) + (4849845.0/1024.0)*pow(x, 9) - 692835.0/256.0*pow(x, 7) + (765765.0/1024.0)*pow(x, 5) - 45045.0/512.0*pow(x, 3) + (3003.0/1024.0)*x;
}


/* P_{13}^{1} */
static double associated_legendre_P_l13m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((16900975.0/1024.0)*pow(x, 12) - 22309287.0/512.0*pow(x, 10) + (43648605.0/1024.0)*pow(x, 8) - 4849845.0/256.0*pow(x, 6) + (3828825.0/1024.0)*pow(x, 4) - 135135.0/512.0*pow(x, 2) + 3003.0/1024.0);
}


/* P_{13}^{2} */
static double associated_legendre_P_l13m2(const double x)
{
  return (1 - pow(x, 2))*((50702925.0/256.0)*pow(x, 11) - 111546435.0/256.0*pow(x, 9) + (43648605.0/128.0)*pow(x, 7) - 14549535.0/128.0*pow(x, 5) + (3828825.0/256.0)*pow(x, 3) - 135135.0/256.0*x);
}


/* P_{13}^{3} */
static double associated_legendre_P_l13m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((557732175.0/256.0)*pow(x, 10) - 1003917915.0/256.0*pow(x, 8) + (305540235.0/128.0)*pow(x, 6) - 72747675.0/128.0*pow(x, 4) + (11486475.0/256.0)*pow(x, 2) - 135135.0/256.0);
}


/* P_{13}^{4} */
static double associated_legendre_P_l13m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((2788660875.0/128.0)*pow(x, 9) - 1003917915.0/32.0*pow(x, 7) + (916620705.0/64.0)*pow(x, 5) - 72747675.0/32.0*pow(x, 3) + (11486475.0/128.0)*x);
}


/* P_{13}^{5} */
static double associated_legendre_P_l13m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((25097947875.0/128.0)*pow(x, 8) - 7027425405.0/32.0*pow(x, 6) + (4583103525.0/64.0)*pow(x, 4) - 218243025.0/32.0*pow(x, 2) + 11486475.0/128.0);
}


/* P_{13}^{6} */
static double associated_legendre_P_l13m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((25097947875.0/16.0)*pow(x, 7) - 21082276215.0/16.0*pow(x, 5) + (4583103525.0/16.0)*pow(x, 3) - 218243025.0/16.0*x);
}


/* P_{13}^{7} */
static double associated_legendre_P_l13m7(const double x)
{
  return -pow(1 - pow(x, 2), 7.0/2.0)*((175685635125.0/16.0)*pow(x, 6) - 105411381075.0/16.0*pow(x, 4) + (13749310575.0/16.0)*pow(x, 2) - 218243025.0/16.0);
}


/* P_{13}^{8} */
static double associated_legendre_P_l13m8(const double x)
{
  return pow(1 - pow(x, 2), 4)*((527056905375.0/8.0)*pow(x, 5) - 105411381075.0/4.0*pow(x, 3) + (13749310575.0/8.0)*x);
}


/* P_{13}^{9} */
static double associated_legendre_P_l13m9(const double x)
{
  return -pow(1 - pow(x, 2), 9.0/2.0)*((2635284526875.0/8.0)*pow(x, 4) - 316234143225.0/4.0*pow(x, 2) + 13749310575.0/8.0);
}


/* P_{13}^{10} */
static double associated_legendre_P_l13m10(const double x)
{
  return pow(1 - pow(x, 2), 5)*((2635284526875.0/2.0)*pow(x, 3) - 316234143225.0/2.0*x);
}


/* P_{13}^{11} */
static double associated_legendre_P_l13m11(const double x)
{
  return -pow(1 - pow(x, 2), 11.0/2.0)*((7905853580625.0/2.0)*pow(x, 2) - 316234143225.0/2.0);
}


/* P_{13}^{12} */
static double associated_legendre_P_l13m12(const double x)
{
  return 7905853580625*x*pow(1 - pow(x, 2), 6);
}


/* P_{13}^{13} */
static double associated_legendre_P_l13m13(const double x)
{
  return -7905853580625*pow(1 - pow(x, 2), 13.0/2.0);
}


/* P_{14}^{0} */
static double associated_legendre_P_l14m0(const double x)
{
  return (5014575.0/2048.0)*pow(x, 14) - 16900975.0/2048.0*pow(x, 12) + (22309287.0/2048.0)*pow(x, 10) - 14549535.0/2048.0*pow(x, 8) + (4849845.0/2048.0)*pow(x, 6) - 765765.0/2048.0*pow(x, 4) + (45045.0/2048.0)*pow(x, 2) - 429.0/2048.0;
}


/* P_{14}^{1} */
static double associated_legendre_P_l14m1(const double x)
{
  return -sqrt(1 - pow(x, 2))*((35102025.0/1024.0)*pow(x, 13) - 50702925.0/512.0*pow(x, 11) + (111546435.0/1024.0)*pow(x, 9) - 14549535.0/256.0*pow(x, 7) + (14549535.0/1024.0)*pow(x, 5) - 765765.0/512.0*pow(x, 3) + (45045.0/1024.0)*x);
}


/* P_{14}^{2} */
static double associated_legendre_P_l14m2(const double x)
{
  return (1 - pow(x, 2))*((456326325.0/1024.0)*pow(x, 12) - 557732175.0/512.0*pow(x, 10) + (1003917915.0/1024.0)*pow(x, 8) - 101846745.0/256.0*pow(x, 6) + (72747675.0/1024.0)*pow(x, 4) - 2297295.0/512.0*pow(x, 2) + 45045.0/1024.0);
}


/* P_{14}^{3} */
static double associated_legendre_P_l14m3(const double x)
{
  return -pow(1 - pow(x, 2), 3.0/2.0)*((1368978975.0/256.0)*pow(x, 11) - 2788660875.0/256.0*pow(x, 9) + (1003917915.0/128.0)*pow(x, 7) - 305540235.0/128.0*pow(x, 5) + (72747675.0/256.0)*pow(x, 3) - 2297295.0/256.0*x);
}


/* P_{14}^{4} */
static double associated_legendre_P_l14m4(const double x)
{
  return pow(1 - pow(x, 2), 2)*((15058768725.0/256.0)*pow(x, 10) - 25097947875.0/256.0*pow(x, 8) + (7027425405.0/128.0)*pow(x, 6) - 1527701175.0/128.0*pow(x, 4) + (218243025.0/256.0)*pow(x, 2) - 2297295.0/256.0);
}


/* P_{14}^{5} */
static double associated_legendre_P_l14m5(const double x)
{
  return -pow(1 - pow(x, 2), 5.0/2.0)*((75293843625.0/128.0)*pow(x, 9) - 25097947875.0/32.0*pow(x, 7) + (21082276215.0/64.0)*pow(x, 5) - 1527701175.0/32.0*pow(x, 3) + (218243025.0/128.0)*x);
}


/* P_{14}^{6} */
static double associated_legendre_P_l14m6(const double x)
{
  return pow(1 - pow(x, 2), 3)*((677644592625.0/128.0)*pow(x, 8) - 175685635125.0/32.0*pow(x, 6) + (105411381075.0/64.0)*pow(x, 4) - 4583103525.0/32.0*pow(x, 2) + 218243025.0/128.0);
}


/* P_{14}^{7} */
static double associated_legendre_P_l14m7(const double x)
{
  return -pow(1 - pow(x, 2), 7.0/2.0)*((677644592625.0/16.0)*pow(x, 7) - 527056905375.0/16.0*pow(x, 5) + (105411381075.0/16.0)*pow(x, 3) - 4583103525.0/16.0*x);
}


/* P_{14}^{8} */
static double associated_legendre_P_l14m8(const double x)
{
  return pow(1 - pow(x, 2), 4)*((4743512148375.0/16.0)*pow(x, 6) - 2635284526875.0/16.0*pow(x, 4) + (316234143225.0/16.0)*pow(x, 2) - 4583103525.0/16.0);
}


/* P_{14}^{9} */
static double associated_legendre_P_l14m9(const double x)
{
  return -pow(1 - pow(x, 2), 9.0/2.0)*((14230536445125.0/8.0)*pow(x, 5) - 2635284526875.0/4.0*pow(x, 3) + (316234143225.0/8.0)*x);
}


/* P_{14}^{10} */
static double associated_legendre_P_l14m10(const double x)
{
  return pow(1 - pow(x, 2), 5)*((71152682225625.0/8.0)*pow(x, 4) - 7905853580625.0/4.0*pow(x, 2) + 316234143225.0/8.0);
}


/* P_{14}^{11} */
static double associated_legendre_P_l14m11(const double x)
{
  return -pow(1 - pow(x, 2), 11.0/2.0)*((71152682225625.0/2.0)*pow(x, 3) - 7905853580625.0/2.0*x);
}


/* P_{14}^{12} */
static double associated_legendre_P_l14m12(const double x)
{
  return pow(1 - pow(x, 2), 6)*((213458046676875.0/2.0)*pow(x, 2) - 7905853580625.0/2.0);
}


/* P_{14}^{13} */
static double associated_legendre_P_l14m13(const double x)
{
  return -213458046676875*x*pow(1 - pow(x, 2), 13.0/2.0);
}


/* P_{14}^{14} */
static double associated_legendre_P_l14m14(const double x)
{
  return 213458046676875*pow(1 - pow(x, 2), 7);
}


/* P_{1}^{-1} */
static double associated_legendre_P_l1m_1(const double x)
{
  return (1.0/2.0)*sqrt(1 - pow(x, 2));
}


/* P_{2}^{-1} */
static double associated_legendre_P_l2m_1(const double x)
{
  return (1.0/2.0)*x*sqrt(1 - pow(x, 2));
}


/* P_{2}^{-2} */
static double associated_legendre_P_l2m_2(const double x)
{
  return 1.0/8.0 - 1.0/8.0*pow(x, 2);
}


/* P_{3}^{-1} */
static double associated_legendre_P_l3m_1(const double x)
{
  return (1.0/12.0)*sqrt(1 - pow(x, 2))*((15.0/2.0)*pow(x, 2) - 3.0/2.0);
}


/* P_{3}^{-2} */
static double associated_legendre_P_l3m_2(const double x)
{
  return (1.0/8.0)*x*(1 - pow(x, 2));
}


/* P_{3}^{-3} */
static double associated_legendre_P_l3m_3(const double x)
{
  return (1.0/48.0)*pow(1 - pow(x, 2), 3.0/2.0);
}


/* P_{4}^{-1} */
static double associated_legendre_P_l4m_1(const double x)
{
  return (1.0/20.0)*sqrt(1 - pow(x, 2))*((35.0/2.0)*pow(x, 3) - 15.0/2.0*x);
}


/* P_{4}^{-2} */
static double associated_legendre_P_l4m_2(const double x)
{
  return (1.0/360.0)*(1 - pow(x, 2))*((105.0/2.0)*pow(x, 2) - 15.0/2.0);
}


/* P_{4}^{-3} */
static double associated_legendre_P_l4m_3(const double x)
{
  return (1.0/48.0)*x*pow(1 - pow(x, 2), 3.0/2.0);
}


/* P_{4}^{-4} */
static double associated_legendre_P_l4m_4(const double x)
{
  return (1.0/384.0)*pow(1 - pow(x, 2), 2);
}


/* P_{5}^{-1} */
static double associated_legendre_P_l5m_1(const double x)
{
  return (1.0/30.0)*sqrt(1 - pow(x, 2))*((315.0/8.0)*pow(x, 4) - 105.0/4.0*pow(x, 2) + 15.0/8.0);
}


/* P_{5}^{-2} */
static double associated_legendre_P_l5m_2(const double x)
{
  return (1.0/840.0)*(1 - pow(x, 2))*((315.0/2.0)*pow(x, 3) - 105.0/2.0*x);
}


/* P_{5}^{-3} */
static double associated_legendre_P_l5m_3(const double x)
{
  return (1.0/20160.0)*pow(1 - pow(x, 2), 3.0/2.0)*((945.0/2.0)*pow(x, 2) - 105.0/2.0);
}


/* P_{5}^{-4} */
static double associated_legendre_P_l5m_4(const double x)
{
  return (1.0/384.0)*x*pow(1 - pow(x, 2), 2);
}


/* P_{5}^{-5} */
static double associated_legendre_P_l5m_5(const double x)
{
  return (1.0/3840.0)*pow(1 - pow(x, 2), 5.0/2.0);
}


/* P_{6}^{-1} */
static double associated_legendre_P_l6m_1(const double x)
{
  return (1.0/42.0)*sqrt(1 - pow(x, 2))*((693.0/8.0)*pow(x, 5) - 315.0/4.0*pow(x, 3) + (105.0/8.0)*x);
}


/* P_{6}^{-2} */
static double associated_legendre_P_l6m_2(const double x)
{
  return (1.0/1680.0)*(1 - pow(x, 2))*((3465.0/8.0)*pow(x, 4) - 945.0/4.0*pow(x, 2) + 105.0/8.0);
}


/* P_{6}^{-3} */
static double associated_legendre_P_l6m_3(const double x)
{
  return (1.0/60480.0)*pow(1 - pow(x, 2), 3.0/2.0)*((3465.0/2.0)*pow(x, 3) - 945.0/2.0*x);
}


/* P_{6}^{-4} */
static double associated_legendre_P_l6m_4(const double x)
{
  return (1.0/1814400.0)*pow(1 - pow(x, 2), 2)*((10395.0/2.0)*pow(x, 2) - 945.0/2.0);
}


/* P_{6}^{-5} */
static double associated_legendre_P_l6m_5(const double x)
{
  return (1.0/3840.0)*x*pow(1 - pow(x, 2), 5.0/2.0);
}


/* P_{6}^{-6} */
static double associated_legendre_P_l6m_6(const double x)
{
  return (1.0/46080.0)*pow(1 - pow(x, 2), 3);
}


/* P_{7}^{-1} */
static double associated_legendre_P_l7m_1(const double x)
{
  return (1.0/56.0)*sqrt(1 - pow(x, 2))*((3003.0/16.0)*pow(x, 6) - 3465.0/16.0*pow(x, 4) + (945.0/16.0)*pow(x, 2) - 35.0/16.0);
}


/* P_{7}^{-2} */
static double associated_legendre_P_l7m_2(const double x)
{
  return (1.0/3024.0)*(1 - pow(x, 2))*((9009.0/8.0)*pow(x, 5) - 3465.0/4.0*pow(x, 3) + (945.0/8.0)*x);
}


/* P_{7}^{-3} */
static double associated_legendre_P_l7m_3(const double x)
{
  return (1.0/151200.0)*pow(1 - pow(x, 2), 3.0/2.0)*((45045.0/8.0)*pow(x, 4) - 10395.0/4.0*pow(x, 2) + 945.0/8.0);
}


/* P_{7}^{-4} */
static double associated_legendre_P_l7m_4(const double x)
{
  return (1.0/6652800.0)*pow(1 - pow(x, 2), 2)*((45045.0/2.0)*pow(x, 3) - 10395.0/2.0*x);
}


/* P_{7}^{-5} */
static double associated_legendre_P_l7m_5(const double x)
{
  return (1.0/239500800.0)*pow(1 - pow(x, 2), 5.0/2.0)*((135135.0/2.0)*pow(x, 2) - 10395.0/2.0);
}


/* P_{7}^{-6} */
static double associated_legendre_P_l7m_6(const double x)
{
  return (1.0/46080.0)*x*pow(1 - pow(x, 2), 3);
}


/* P_{7}^{-7} */
static double associated_legendre_P_l7m_7(const double x)
{
  return (1.0/645120.0)*pow(1 - pow(x, 2), 7.0/2.0);
}


/* P_{8}^{-1} */
static double associated_legendre_P_l8m_1(const double x)
{
  return (1.0/72.0)*sqrt(1 - pow(x, 2))*((6435.0/16.0)*pow(x, 7) - 9009.0/16.0*pow(x, 5) + (3465.0/16.0)*pow(x, 3) - 315.0/16.0*x);
}


/* P_{8}^{-2} */
static double associated_legendre_P_l8m_2(const double x)
{
  return (1.0/5040.0)*(1 - pow(x, 2))*((45045.0/16.0)*pow(x, 6) - 45045.0/16.0*pow(x, 4) + (10395.0/16.0)*pow(x, 2) - 315.0/16.0);
}


/* P_{8}^{-3} */
static double associated_legendre_P_l8m_3(const double x)
{
  return (1.0/332640.0)*pow(1 - pow(x, 2), 3.0/2.0)*((135135.0/8.0)*pow(x, 5) - 45045.0/4.0*pow(x, 3) + (10395.0/8.0)*x);
}


/* P_{8}^{-4} */
static double associated_legendre_P_l8m_4(const double x)
{
  return (1.0/19958400.0)*pow(1 - pow(x, 2), 2)*((675675.0/8.0)*pow(x, 4) - 135135.0/4.0*pow(x, 2) + 10395.0/8.0);
}


/* P_{8}^{-5} */
static double associated_legendre_P_l8m_5(const double x)
{
  return (1.0/1037836800.0)*pow(1 - pow(x, 2), 5.0/2.0)*((675675.0/2.0)*pow(x, 3) - 135135.0/2.0*x);
}


/* P_{8}^{-6} */
static double associated_legendre_P_l8m_6(const double x)
{
  return (1.0/43589145600.0)*pow(1 - pow(x, 2), 3)*((2027025.0/2.0)*pow(x, 2) - 135135.0/2.0);
}


/* P_{8}^{-7} */
static double associated_legendre_P_l8m_7(const double x)
{
  return (1.0/645120.0)*x*pow(1 - pow(x, 2), 7.0/2.0);
}


/* P_{8}^{-8} */
static double associated_legendre_P_l8m_8(const double x)
{
  return (1.0/10321920.0)*pow(1 - pow(x, 2), 4);
}


/* P_{9}^{-1} */
static double associated_legendre_P_l9m_1(const double x)
{
  return (1.0/90.0)*sqrt(1 - pow(x, 2))*((109395.0/128.0)*pow(x, 8) - 45045.0/32.0*pow(x, 6) + (45045.0/64.0)*pow(x, 4) - 3465.0/32.0*pow(x, 2) + 315.0/128.0);
}


/* P_{9}^{-2} */
static double associated_legendre_P_l9m_2(const double x)
{
  return (1.0/7920.0)*(1 - pow(x, 2))*((109395.0/16.0)*pow(x, 7) - 135135.0/16.0*pow(x, 5) + (45045.0/16.0)*pow(x, 3) - 3465.0/16.0*x);
}


/* P_{9}^{-3} */
static double associated_legendre_P_l9m_3(const double x)
{
  return (1.0/665280.0)*pow(1 - pow(x, 2), 3.0/2.0)*((765765.0/16.0)*pow(x, 6) - 675675.0/16.0*pow(x, 4) + (135135.0/16.0)*pow(x, 2) - 3465.0/16.0);
}


/* P_{9}^{-4} */
static double associated_legendre_P_l9m_4(const double x)
{
  return (1.0/51891840.0)*pow(1 - pow(x, 2), 2)*((2297295.0/8.0)*pow(x, 5) - 675675.0/4.0*pow(x, 3) + (135135.0/8.0)*x);
}


/* P_{9}^{-5} */
static double associated_legendre_P_l9m_5(const double x)
{
  return (1.0/3632428800.0)*pow(1 - pow(x, 2), 5.0/2.0)*((11486475.0/8.0)*pow(x, 4) - 2027025.0/4.0*pow(x, 2) + 135135.0/8.0);
}


/* P_{9}^{-6} */
static double associated_legendre_P_l9m_6(const double x)
{
  return (1.0/217945728000.0)*pow(1 - pow(x, 2), 3)*((11486475.0/2.0)*pow(x, 3) - 2027025.0/2.0*x);
}


/* P_{9}^{-7} */
static double associated_legendre_P_l9m_7(const double x)
{
  return (1.0/10461394944000.0)*pow(1 - pow(x, 2), 7.0/2.0)*((34459425.0/2.0)*pow(x, 2) - 2027025.0/2.0);
}


/* P_{9}^{-8} */
static double associated_legendre_P_l9m_8(const double x)
{
  return (1.0/10321920.0)*x*pow(1 - pow(x, 2), 4);
}


/* P_{9}^{-9} */
static double associated_legendre_P_l9m_9(const double x)
{
  return (1.0/185794560.0)*pow(1 - pow(x, 2), 9.0/2.0);
}


/* P_{10}^{-1} */
static double associated_legendre_P_l10m_1(const double x)
{
  return (1.0/110.0)*sqrt(1 - pow(x, 2))*((230945.0/128.0)*pow(x, 9) - 109395.0/32.0*pow(x, 7) + (135135.0/64.0)*pow(x, 5) - 15015.0/32.0*pow(x, 3) + (3465.0/128.0)*x);
}


/* P_{10}^{-2} */
static double associated_legendre_P_l10m_2(const double x)
{
  return (1.0/11880.0)*(1 - pow(x, 2))*((2078505.0/128.0)*pow(x, 8) - 765765.0/32.0*pow(x, 6) + (675675.0/64.0)*pow(x, 4) - 45045.0/32.0*pow(x, 2) + 3465.0/128.0);
}


/* P_{10}^{-3} */
static double associated_legendre_P_l10m_3(const double x)
{
  return (1.0/1235520.0)*pow(1 - pow(x, 2), 3.0/2.0)*((2078505.0/16.0)*pow(x, 7) - 2297295.0/16.0*pow(x, 5) + (675675.0/16.0)*pow(x, 3) - 45045.0/16.0*x);
}


/* P_{10}^{-4} */
static double associated_legendre_P_l10m_4(const double x)
{
  return (1.0/121080960.0)*pow(1 - pow(x, 2), 2)*((14549535.0/16.0)*pow(x, 6) - 11486475.0/16.0*pow(x, 4) + (2027025.0/16.0)*pow(x, 2) - 45045.0/16.0);
}


/* P_{10}^{-5} */
static double associated_legendre_P_l10m_5(const double x)
{
  return (1.0/10897286400.0)*pow(1 - pow(x, 2), 5.0/2.0)*((43648605.0/8.0)*pow(x, 5) - 11486475.0/4.0*pow(x, 3) + (2027025.0/8.0)*x);
}


/* P_{10}^{-6} */
static double associated_legendre_P_l10m_6(const double x)
{
  return (1.0/871782912000.0)*pow(1 - pow(x, 2), 3)*((218243025.0/8.0)*pow(x, 4) - 34459425.0/4.0*pow(x, 2) + 2027025.0/8.0);
}


/* P_{10}^{-7} */
static double associated_legendre_P_l10m_7(const double x)
{
  return (1.0/59281238016000.0)*pow(1 - pow(x, 2), 7.0/2.0)*((218243025.0/2.0)*pow(x, 3) - 34459425.0/2.0*x);
}


/* P_{10}^{-8} */
static double associated_legendre_P_l10m_8(const double x)
{
  return (1.0/3201186852864000.0)*pow(1 - pow(x, 2), 4)*((654729075.0/2.0)*pow(x, 2) - 34459425.0/2.0);
}


/* P_{10}^{-9} */
static double associated_legendre_P_l10m_9(const double x)
{
  return (1.0/185794560.0)*x*pow(1 - pow(x, 2), 9.0/2.0);
}


/* P_{10}^{-10} */
static double associated_legendre_P_l10m_10(const double x)
{
  return (1.0/3715891200.0)*pow(1 - pow(x, 2), 5);
}


/* P_{11}^{-1} */
static double associated_legendre_P_l11m_1(const double x)
{
  return (1.0/132.0)*sqrt(1 - pow(x, 2))*((969969.0/256.0)*pow(x, 10) - 2078505.0/256.0*pow(x, 8) + (765765.0/128.0)*pow(x, 6) - 225225.0/128.0*pow(x, 4) + (45045.0/256.0)*pow(x, 2) - 693.0/256.0);
}


/* P_{11}^{-2} */
static double associated_legendre_P_l11m_2(const double x)
{
  return (1.0/17160.0)*(1 - pow(x, 2))*((4849845.0/128.0)*pow(x, 9) - 2078505.0/32.0*pow(x, 7) + (2297295.0/64.0)*pow(x, 5) - 225225.0/32.0*pow(x, 3) + (45045.0/128.0)*x);
}


/* P_{11}^{-3} */
static double associated_legendre_P_l11m_3(const double x)
{
  return (1.0/2162160.0)*pow(1 - pow(x, 2), 3.0/2.0)*((43648605.0/128.0)*pow(x, 8) - 14549535.0/32.0*pow(x, 6) + (11486475.0/64.0)*pow(x, 4) - 675675.0/32.0*pow(x, 2) + 45045.0/128.0);
}


/* P_{11}^{-4} */
static double associated_legendre_P_l11m_4(const double x)
{
  return (1.0/259459200.0)*pow(1 - pow(x, 2), 2)*((43648605.0/16.0)*pow(x, 7) - 43648605.0/16.0*pow(x, 5) + (11486475.0/16.0)*pow(x, 3) - 675675.0/16.0*x);
}


/* P_{11}^{-5} */
static double associated_legendre_P_l11m_5(const double x)
{
  return (1.0/29059430400.0)*pow(1 - pow(x, 2), 5.0/2.0)*((305540235.0/16.0)*pow(x, 6) - 218243025.0/16.0*pow(x, 4) + (34459425.0/16.0)*pow(x, 2) - 675675.0/16.0);
}


/* P_{11}^{-6} */
static double associated_legendre_P_l11m_6(const double x)
{
  return (1.0/2964061900800.0)*pow(1 - pow(x, 2), 3)*((916620705.0/8.0)*pow(x, 5) - 218243025.0/4.0*pow(x, 3) + (34459425.0/8.0)*x);
}


/* P_{11}^{-7} */
static double associated_legendre_P_l11m_7(const double x)
{
  return (1.0/266765571072000.0)*pow(1 - pow(x, 2), 7.0/2.0)*((4583103525.0/8.0)*pow(x, 4) - 654729075.0/4.0*pow(x, 2) + 34459425.0/8.0);
}


/* P_{11}^{-8} */
static double associated_legendre_P_l11m_8(const double x)
{
  return (1.0/20274183401472000.0)*pow(1 - pow(x, 2), 4)*((4583103525.0/2.0)*pow(x, 3) - 654729075.0/2.0*x);
}


/* P_{11}^{-9} */
static double associated_legendre_P_l11m_9(const double x)
{
  return (1.0/1216451004088320000.0)*pow(1 - pow(x, 2), 9.0/2.0)*((13749310575.0/2.0)*pow(x, 2) - 654729075.0/2.0);
}


/* P_{11}^{-10} */
static double associated_legendre_P_l11m_10(const double x)
{
  return (1.0/3715891200.0)*x*pow(1 - pow(x, 2), 5);
}


/* P_{11}^{-11} */
static double associated_legendre_P_l11m_11(const double x)
{
  return (1.0/81749606400.0)*pow(1 - pow(x, 2), 11.0/2.0);
}


/* P_{12}^{-1} */
static double associated_legendre_P_l12m_1(const double x)
{
  return (1.0/156.0)*sqrt(1 - pow(x, 2))*((2028117.0/256.0)*pow(x, 11) - 4849845.0/256.0*pow(x, 9) + (2078505.0/128.0)*pow(x, 7) - 765765.0/128.0*pow(x, 5) + (225225.0/256.0)*pow(x, 3) - 9009.0/256.0*x);
}


/* P_{12}^{-2} */
static double associated_legendre_P_l12m_2(const double x)
{
  return (1.0/24024.0)*(1 - pow(x, 2))*((22309287.0/256.0)*pow(x, 10) - 43648605.0/256.0*pow(x, 8) + (14549535.0/128.0)*pow(x, 6) - 3828825.0/128.0*pow(x, 4) + (675675.0/256.0)*pow(x, 2) - 9009.0/256.0);
}


/* P_{12}^{-3} */
static double associated_legendre_P_l12m_3(const double x)
{
  return (1.0/3603600.0)*pow(1 - pow(x, 2), 3.0/2.0)*((111546435.0/128.0)*pow(x, 9) - 43648605.0/32.0*pow(x, 7) + (43648605.0/64.0)*pow(x, 5) - 3828825.0/32.0*pow(x, 3) + (675675.0/128.0)*x);
}


/* P_{12}^{-4} */
static double associated_legendre_P_l12m_4(const double x)
{
  return (1.0/518918400.0)*pow(1 - pow(x, 2), 2)*((1003917915.0/128.0)*pow(x, 8) - 305540235.0/32.0*pow(x, 6) + (218243025.0/64.0)*pow(x, 4) - 11486475.0/32.0*pow(x, 2) + 675675.0/128.0);
}


/* P_{12}^{-5} */
static double associated_legendre_P_l12m_5(const double x)
{
  return (1.0/70572902400.0)*pow(1 - pow(x, 2), 5.0/2.0)*((1003917915.0/16.0)*pow(x, 7) - 916620705.0/16.0*pow(x, 5) + (218243025.0/16.0)*pow(x, 3) - 11486475.0/16.0*x);
}


/* P_{12}^{-6} */
static double associated_legendre_P_l12m_6(const double x)
{
  return (1.0/8892185702400.0)*pow(1 - pow(x, 2), 3)*((7027425405.0/16.0)*pow(x, 6) - 4583103525.0/16.0*pow(x, 4) + (654729075.0/16.0)*pow(x, 2) - 11486475.0/16.0);
}


/* P_{12}^{-7} */
static double associated_legendre_P_l12m_7(const double x)
{
  return (1.0/1013709170073600.0)*pow(1 - pow(x, 2), 7.0/2.0)*((21082276215.0/8.0)*pow(x, 5) - 4583103525.0/4.0*pow(x, 3) + (654729075.0/8.0)*x);
}


/* P_{12}^{-8} */
static double associated_legendre_P_l12m_8(const double x)
{
  return (1.0/101370917007360000.0)*pow(1 - pow(x, 2), 4)*((105411381075.0/8.0)*pow(x, 4) - 13749310575.0/4.0*pow(x, 2) + 654729075.0/8.0);
}


/* P_{12}^{-9} */
static double associated_legendre_P_l12m_9(const double x)
{
  return (1.0/8515157028618240000.0)*pow(1 - pow(x, 2), 9.0/2.0)*((105411381075.0/2.0)*pow(x, 3) - 13749310575.0/2.0*x);
}


/* P_{12}^{-10} */
static double associated_legendre_P_l12m_10(const double x)
{
  return (1.0/562000363888803840000.0)*pow(1 - pow(x, 2), 5)*((316234143225.0/2.0)*pow(x, 2) - 13749310575.0/2.0);
}


/* P_{12}^{-11} */
static double associated_legendre_P_l12m_11(const double x)
{
  return (1.0/81749606400.0)*x*pow(1 - pow(x, 2), 11.0/2.0);
}


/* P_{12}^{-12} */
static double associated_legendre_P_l12m_12(const double x)
{
  return (1.0/1961990553600.0)*pow(1 - pow(x, 2), 6);
}


/* P_{13}^{-1} */
static double associated_legendre_P_l13m_1(const double x)
{
  return (1.0/182.0)*sqrt(1 - pow(x, 2))*((16900975.0/1024.0)*pow(x, 12) - 22309287.0/512.0*pow(x, 10) + (43648605.0/1024.0)*pow(x, 8) - 4849845.0/256.0*pow(x, 6) + (3828825.0/1024.0)*pow(x, 4) - 135135.0/512.0*pow(x, 2) + 3003.0/1024.0);
}


/* P_{13}^{-2} */
static double associated_legendre_P_l13m_2(const double x)
{
  return (1.0/32760.0)*(1 - pow(x, 2))*((50702925.0/256.0)*pow(x, 11) - 111546435.0/256.0*pow(x, 9) + (43648605.0/128.0)*pow(x, 7) - 14549535.0/128.0*pow(x, 5) + (3828825.0/256.0)*pow(x, 3) - 135135.0/256.0*x);
}


/* P_{13}^{-3} */
static double associated_legendre_P_l13m_3(const double x)
{
  return (1.0/5765760.0)*pow(1 - pow(x, 2), 3.0/2.0)*((557732175.0/256.0)*pow(x, 10) - 1003917915.0/256.0*pow(x, 8) + (305540235.0/128.0)*pow(x, 6) - 72747675.0/128.0*pow(x, 4) + (11486475.0/256.0)*pow(x, 2) - 135135.0/256.0);
}


/* P_{13}^{-4} */
static double associated_legendre_P_l13m_4(const double x)
{
  return (1.0/980179200.0)*pow(1 - pow(x, 2), 2)*((2788660875.0/128.0)*pow(x, 9) - 1003917915.0/32.0*pow(x, 7) + (916620705.0/64.0)*pow(x, 5) - 72747675.0/32.0*pow(x, 3) + (11486475.0/128.0)*x);
}


/* P_{13}^{-5} */
static double associated_legendre_P_l13m_5(const double x)
{
  return (1.0/158789030400.0)*pow(1 - pow(x, 2), 5.0/2.0)*((25097947875.0/128.0)*pow(x, 8) - 7027425405.0/32.0*pow(x, 6) + (4583103525.0/64.0)*pow(x, 4) - 218243025.0/32.0*pow(x, 2) + 11486475.0/128.0);
}


/* P_{13}^{-6} */
static double associated_legendre_P_l13m_6(const double x)
{
  return (1.0/24135932620800.0)*pow(1 - pow(x, 2), 3)*((25097947875.0/16.0)*pow(x, 7) - 21082276215.0/16.0*pow(x, 5) + (4583103525.0/16.0)*pow(x, 3) - 218243025.0/16.0*x);
}


/* P_{13}^{-7} */
static double associated_legendre_P_l13m_7(const double x)
{
  return (1.0/3379030566912000.0)*pow(1 - pow(x, 2), 7.0/2.0)*((175685635125.0/16.0)*pow(x, 6) - 105411381075.0/16.0*pow(x, 4) + (13749310575.0/16.0)*pow(x, 2) - 218243025.0/16.0);
}


/* P_{13}^{-8} */
static double associated_legendre_P_l13m_8(const double x)
{
  return (1.0/425757851430912000.0)*pow(1 - pow(x, 2), 4)*((527056905375.0/8.0)*pow(x, 5) - 105411381075.0/4.0*pow(x, 3) + (13749310575.0/8.0)*x);
}


/* P_{13}^{-9} */
static double associated_legendre_P_l13m_9(const double x)
{
  return (1.0/46833363657400320000.0)*pow(1 - pow(x, 2), 9.0/2.0)*((2635284526875.0/8.0)*pow(x, 4) - 316234143225.0/4.0*pow(x, 2) + 13749310575.0/8.0);
}


/* P_{13}^{-10} */
static double associated_legendre_P_l13m_10(const double x)
{
  return (1.0/4308669456480829440000.0)*pow(1 - pow(x, 2), 5)*((2635284526875.0/2.0)*pow(x, 3) - 316234143225.0/2.0*x);
}


/* P_{13}^{-11} */
static double associated_legendre_P_l13m_11(const double x)
{
  return (1.0/310224200866619719680000.0)*pow(1 - pow(x, 2), 11.0/2.0)*((7905853580625.0/2.0)*pow(x, 2) - 316234143225.0/2.0);
}


/* P_{13}^{-12} */
static double associated_legendre_P_l13m_12(const double x)
{
  return (1.0/1961990553600.0)*x*pow(1 - pow(x, 2), 6);
}


/* P_{13}^{-13} */
static double associated_legendre_P_l13m_13(const double x)
{
  return (1.0/51011754393600.0)*pow(1 - pow(x, 2), 13.0/2.0);
}


/* P_{14}^{-1} */
static double associated_legendre_P_l14m_1(const double x)
{
  return (1.0/210.0)*sqrt(1 - pow(x, 2))*((35102025.0/1024.0)*pow(x, 13) - 50702925.0/512.0*pow(x, 11) + (111546435.0/1024.0)*pow(x, 9) - 14549535.0/256.0*pow(x, 7) + (14549535.0/1024.0)*pow(x, 5) - 765765.0/512.0*pow(x, 3) + (45045.0/1024.0)*x);
}


/* P_{14}^{-2} */
static double associated_legendre_P_l14m_2(const double x)
{
  return (1.0/43680.0)*(1 - pow(x, 2))*((456326325.0/1024.0)*pow(x, 12) - 557732175.0/512.0*pow(x, 10) + (1003917915.0/1024.0)*pow(x, 8) - 101846745.0/256.0*pow(x, 6) + (72747675.0/1024.0)*pow(x, 4) - 2297295.0/512.0*pow(x, 2) + 45045.0/1024.0);
}


/* P_{14}^{-3} */
static double associated_legendre_P_l14m_3(const double x)
{
  return (1.0/8910720.0)*pow(1 - pow(x, 2), 3.0/2.0)*((1368978975.0/256.0)*pow(x, 11) - 2788660875.0/256.0*pow(x, 9) + (1003917915.0/128.0)*pow(x, 7) - 305540235.0/128.0*pow(x, 5) + (72747675.0/256.0)*pow(x, 3) - 2297295.0/256.0*x);
}


/* P_{14}^{-4} */
static double associated_legendre_P_l14m_4(const double x)
{
  return (1.0/1764322560.0)*pow(1 - pow(x, 2), 2)*((15058768725.0/256.0)*pow(x, 10) - 25097947875.0/256.0*pow(x, 8) + (7027425405.0/128.0)*pow(x, 6) - 1527701175.0/128.0*pow(x, 4) + (218243025.0/256.0)*pow(x, 2) - 2297295.0/256.0);
}


/* P_{14}^{-5} */
static double associated_legendre_P_l14m_5(const double x)
{
  return (1.0/335221286400.0)*pow(1 - pow(x, 2), 5.0/2.0)*((75293843625.0/128.0)*pow(x, 9) - 25097947875.0/32.0*pow(x, 7) + (21082276215.0/64.0)*pow(x, 5) - 1527701175.0/32.0*pow(x, 3) + (218243025.0/128.0)*x);
}


/* P_{14}^{-6} */
static double associated_legendre_P_l14m_6(const double x)
{
  return (1.0/60339831552000.0)*pow(1 - pow(x, 2), 3)*((677644592625.0/128.0)*pow(x, 8) - 175685635125.0/32.0*pow(x, 6) + (105411381075.0/64.0)*pow(x, 4) - 4583103525.0/32.0*pow(x, 2) + 218243025.0/128.0);
}


/* P_{14}^{-7} */
static double associated_legendre_P_l14m_7(const double x)
{
  return (1.0/10137091700736000.0)*pow(1 - pow(x, 2), 7.0/2.0)*((677644592625.0/16.0)*pow(x, 7) - 527056905375.0/16.0*pow(x, 5) + (105411381075.0/16.0)*pow(x, 3) - 4583103525.0/16.0*x);
}


/* P_{14}^{-8} */
static double associated_legendre_P_l14m_8(const double x)
{
  return (1.0/1561112121913344000.0)*pow(1 - pow(x, 2), 4)*((4743512148375.0/16.0)*pow(x, 6) - 2635284526875.0/16.0*pow(x, 4) + (316234143225.0/16.0)*pow(x, 2) - 4583103525.0/16.0);
}


/* P_{14}^{-9} */
static double associated_legendre_P_l14m_9(const double x)
{
  return (1.0/215433472824041472000.0)*pow(1 - pow(x, 2), 9.0/2.0)*((14230536445125.0/8.0)*pow(x, 5) - 2635284526875.0/4.0*pow(x, 3) + (316234143225.0/8.0)*x);
}


/* P_{14}^{-10} */
static double associated_legendre_P_l14m_10(const double x)
{
  return (1.0/25852016738884976640000.0)*pow(1 - pow(x, 2), 5)*((71152682225625.0/8.0)*pow(x, 4) - 7905853580625.0/4.0*pow(x, 2) + 316234143225.0/8.0);
}


/* P_{14}^{-11} */
static double associated_legendre_P_l14m_11(const double x)
{
  return (1.0/2585201673888497664000000.0)*pow(1 - pow(x, 2), 11.0/2.0)*((71152682225625.0/2.0)*pow(x, 3) - 7905853580625.0/2.0*x);
}


/* P_{14}^{-12} */
static double associated_legendre_P_l14m_12(const double x)
{
  return (1.0/201645730563302817792000000.0)*pow(1 - pow(x, 2), 6)*((213458046676875.0/2.0)*pow(x, 2) - 7905853580625.0/2.0);
}


/* P_{14}^{-13} */
static double associated_legendre_P_l14m_13(const double x)
{
  return (1.0/51011754393600.0)*x*pow(1 - pow(x, 2), 13.0/2.0);
}


/* P_{14}^{-14} */
static double associated_legendre_P_l14m_14(const double x)
{
  return (1.0/1428329123020800.0)*pow(1 - pow(x, 2), 7);
}

