 /* direct output of this file to spherical_harmonic_Ylm.c */
#include "core_lib.h"
#include "maths_analytic_lib.h"
#include "maths_general_lib.h"
#include <complex.h>

#define LMAX_ARRAY_SIZE_YLM 15
typedef double complex fYlm_T (const double theta, const double phi);
fYlm_T *Y[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];
fYlm_T *Y_[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];
fYlm_T *dY_dphi[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];
fYlm_T *dY_dphi_[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];
fYlm_T *dY_dtheta[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];
fYlm_T *dY_dtheta_[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];

void init_Ylm(void);
void init_dYlm_dphi(void);
void init_dYlm_dtheta(void);
double complex Ylm(const unsigned l, const int m, const double theta, const double phi);
double complex dYlm_dphi(const unsigned l, const int m, const double theta, const double phi);
double complex dYlm_dtheta(const unsigned l, const int m, const double theta, const double phi);
static double complex Ylm_l0m0(const double theta, const double phi);
static double complex Ylm_l1m0(const double theta, const double phi);
static double complex Ylm_l1m1(const double theta, const double phi);
static double complex Ylm_l2m0(const double theta, const double phi);
static double complex Ylm_l2m1(const double theta, const double phi);
static double complex Ylm_l2m2(const double theta, const double phi);
static double complex Ylm_l3m0(const double theta, const double phi);
static double complex Ylm_l3m1(const double theta, const double phi);
static double complex Ylm_l3m2(const double theta, const double phi);
static double complex Ylm_l3m3(const double theta, const double phi);
static double complex Ylm_l4m0(const double theta, const double phi);
static double complex Ylm_l4m1(const double theta, const double phi);
static double complex Ylm_l4m2(const double theta, const double phi);
static double complex Ylm_l4m3(const double theta, const double phi);
static double complex Ylm_l4m4(const double theta, const double phi);
static double complex Ylm_l5m0(const double theta, const double phi);
static double complex Ylm_l5m1(const double theta, const double phi);
static double complex Ylm_l5m2(const double theta, const double phi);
static double complex Ylm_l5m3(const double theta, const double phi);
static double complex Ylm_l5m4(const double theta, const double phi);
static double complex Ylm_l5m5(const double theta, const double phi);
static double complex Ylm_l6m0(const double theta, const double phi);
static double complex Ylm_l6m1(const double theta, const double phi);
static double complex Ylm_l6m2(const double theta, const double phi);
static double complex Ylm_l6m3(const double theta, const double phi);
static double complex Ylm_l6m4(const double theta, const double phi);
static double complex Ylm_l6m5(const double theta, const double phi);
static double complex Ylm_l6m6(const double theta, const double phi);
static double complex Ylm_l7m0(const double theta, const double phi);
static double complex Ylm_l7m1(const double theta, const double phi);
static double complex Ylm_l7m2(const double theta, const double phi);
static double complex Ylm_l7m3(const double theta, const double phi);
static double complex Ylm_l7m4(const double theta, const double phi);
static double complex Ylm_l7m5(const double theta, const double phi);
static double complex Ylm_l7m6(const double theta, const double phi);
static double complex Ylm_l7m7(const double theta, const double phi);
static double complex Ylm_l8m0(const double theta, const double phi);
static double complex Ylm_l8m1(const double theta, const double phi);
static double complex Ylm_l8m2(const double theta, const double phi);
static double complex Ylm_l8m3(const double theta, const double phi);
static double complex Ylm_l8m4(const double theta, const double phi);
static double complex Ylm_l8m5(const double theta, const double phi);
static double complex Ylm_l8m6(const double theta, const double phi);
static double complex Ylm_l8m7(const double theta, const double phi);
static double complex Ylm_l8m8(const double theta, const double phi);
static double complex Ylm_l9m0(const double theta, const double phi);
static double complex Ylm_l9m1(const double theta, const double phi);
static double complex Ylm_l9m2(const double theta, const double phi);
static double complex Ylm_l9m3(const double theta, const double phi);
static double complex Ylm_l9m4(const double theta, const double phi);
static double complex Ylm_l9m5(const double theta, const double phi);
static double complex Ylm_l9m6(const double theta, const double phi);
static double complex Ylm_l9m7(const double theta, const double phi);
static double complex Ylm_l9m8(const double theta, const double phi);
static double complex Ylm_l9m9(const double theta, const double phi);
static double complex Ylm_l10m0(const double theta, const double phi);
static double complex Ylm_l10m1(const double theta, const double phi);
static double complex Ylm_l10m2(const double theta, const double phi);
static double complex Ylm_l10m3(const double theta, const double phi);
static double complex Ylm_l10m4(const double theta, const double phi);
static double complex Ylm_l10m5(const double theta, const double phi);
static double complex Ylm_l10m6(const double theta, const double phi);
static double complex Ylm_l10m7(const double theta, const double phi);
static double complex Ylm_l10m8(const double theta, const double phi);
static double complex Ylm_l10m9(const double theta, const double phi);
static double complex Ylm_l10m10(const double theta, const double phi);
static double complex Ylm_l11m0(const double theta, const double phi);
static double complex Ylm_l11m1(const double theta, const double phi);
static double complex Ylm_l11m2(const double theta, const double phi);
static double complex Ylm_l11m3(const double theta, const double phi);
static double complex Ylm_l11m4(const double theta, const double phi);
static double complex Ylm_l11m5(const double theta, const double phi);
static double complex Ylm_l11m6(const double theta, const double phi);
static double complex Ylm_l11m7(const double theta, const double phi);
static double complex Ylm_l11m8(const double theta, const double phi);
static double complex Ylm_l11m9(const double theta, const double phi);
static double complex Ylm_l11m10(const double theta, const double phi);
static double complex Ylm_l11m11(const double theta, const double phi);
static double complex Ylm_l12m0(const double theta, const double phi);
static double complex Ylm_l12m1(const double theta, const double phi);
static double complex Ylm_l12m2(const double theta, const double phi);
static double complex Ylm_l12m3(const double theta, const double phi);
static double complex Ylm_l12m4(const double theta, const double phi);
static double complex Ylm_l12m5(const double theta, const double phi);
static double complex Ylm_l12m6(const double theta, const double phi);
static double complex Ylm_l12m7(const double theta, const double phi);
static double complex Ylm_l12m8(const double theta, const double phi);
static double complex Ylm_l12m9(const double theta, const double phi);
static double complex Ylm_l12m10(const double theta, const double phi);
static double complex Ylm_l12m11(const double theta, const double phi);
static double complex Ylm_l12m12(const double theta, const double phi);
static double complex Ylm_l13m0(const double theta, const double phi);
static double complex Ylm_l13m1(const double theta, const double phi);
static double complex Ylm_l13m2(const double theta, const double phi);
static double complex Ylm_l13m3(const double theta, const double phi);
static double complex Ylm_l13m4(const double theta, const double phi);
static double complex Ylm_l13m5(const double theta, const double phi);
static double complex Ylm_l13m6(const double theta, const double phi);
static double complex Ylm_l13m7(const double theta, const double phi);
static double complex Ylm_l13m8(const double theta, const double phi);
static double complex Ylm_l13m9(const double theta, const double phi);
static double complex Ylm_l13m10(const double theta, const double phi);
static double complex Ylm_l13m11(const double theta, const double phi);
static double complex Ylm_l13m12(const double theta, const double phi);
static double complex Ylm_l13m13(const double theta, const double phi);
static double complex Ylm_l14m0(const double theta, const double phi);
static double complex Ylm_l14m1(const double theta, const double phi);
static double complex Ylm_l14m2(const double theta, const double phi);
static double complex Ylm_l14m3(const double theta, const double phi);
static double complex Ylm_l14m4(const double theta, const double phi);
static double complex Ylm_l14m5(const double theta, const double phi);
static double complex Ylm_l14m6(const double theta, const double phi);
static double complex Ylm_l14m7(const double theta, const double phi);
static double complex Ylm_l14m8(const double theta, const double phi);
static double complex Ylm_l14m9(const double theta, const double phi);
static double complex Ylm_l14m10(const double theta, const double phi);
static double complex Ylm_l14m11(const double theta, const double phi);
static double complex Ylm_l14m12(const double theta, const double phi);
static double complex Ylm_l14m13(const double theta, const double phi);
static double complex Ylm_l14m14(const double theta, const double phi);
static double complex Ylm_l1m_1(const double theta, const double phi);
static double complex Ylm_l2m_1(const double theta, const double phi);
static double complex Ylm_l2m_2(const double theta, const double phi);
static double complex Ylm_l3m_1(const double theta, const double phi);
static double complex Ylm_l3m_2(const double theta, const double phi);
static double complex Ylm_l3m_3(const double theta, const double phi);
static double complex Ylm_l4m_1(const double theta, const double phi);
static double complex Ylm_l4m_2(const double theta, const double phi);
static double complex Ylm_l4m_3(const double theta, const double phi);
static double complex Ylm_l4m_4(const double theta, const double phi);
static double complex Ylm_l5m_1(const double theta, const double phi);
static double complex Ylm_l5m_2(const double theta, const double phi);
static double complex Ylm_l5m_3(const double theta, const double phi);
static double complex Ylm_l5m_4(const double theta, const double phi);
static double complex Ylm_l5m_5(const double theta, const double phi);
static double complex Ylm_l6m_1(const double theta, const double phi);
static double complex Ylm_l6m_2(const double theta, const double phi);
static double complex Ylm_l6m_3(const double theta, const double phi);
static double complex Ylm_l6m_4(const double theta, const double phi);
static double complex Ylm_l6m_5(const double theta, const double phi);
static double complex Ylm_l6m_6(const double theta, const double phi);
static double complex Ylm_l7m_1(const double theta, const double phi);
static double complex Ylm_l7m_2(const double theta, const double phi);
static double complex Ylm_l7m_3(const double theta, const double phi);
static double complex Ylm_l7m_4(const double theta, const double phi);
static double complex Ylm_l7m_5(const double theta, const double phi);
static double complex Ylm_l7m_6(const double theta, const double phi);
static double complex Ylm_l7m_7(const double theta, const double phi);
static double complex Ylm_l8m_1(const double theta, const double phi);
static double complex Ylm_l8m_2(const double theta, const double phi);
static double complex Ylm_l8m_3(const double theta, const double phi);
static double complex Ylm_l8m_4(const double theta, const double phi);
static double complex Ylm_l8m_5(const double theta, const double phi);
static double complex Ylm_l8m_6(const double theta, const double phi);
static double complex Ylm_l8m_7(const double theta, const double phi);
static double complex Ylm_l8m_8(const double theta, const double phi);
static double complex Ylm_l9m_1(const double theta, const double phi);
static double complex Ylm_l9m_2(const double theta, const double phi);
static double complex Ylm_l9m_3(const double theta, const double phi);
static double complex Ylm_l9m_4(const double theta, const double phi);
static double complex Ylm_l9m_5(const double theta, const double phi);
static double complex Ylm_l9m_6(const double theta, const double phi);
static double complex Ylm_l9m_7(const double theta, const double phi);
static double complex Ylm_l9m_8(const double theta, const double phi);
static double complex Ylm_l9m_9(const double theta, const double phi);
static double complex Ylm_l10m_1(const double theta, const double phi);
static double complex Ylm_l10m_2(const double theta, const double phi);
static double complex Ylm_l10m_3(const double theta, const double phi);
static double complex Ylm_l10m_4(const double theta, const double phi);
static double complex Ylm_l10m_5(const double theta, const double phi);
static double complex Ylm_l10m_6(const double theta, const double phi);
static double complex Ylm_l10m_7(const double theta, const double phi);
static double complex Ylm_l10m_8(const double theta, const double phi);
static double complex Ylm_l10m_9(const double theta, const double phi);
static double complex Ylm_l10m_10(const double theta, const double phi);
static double complex Ylm_l11m_1(const double theta, const double phi);
static double complex Ylm_l11m_2(const double theta, const double phi);
static double complex Ylm_l11m_3(const double theta, const double phi);
static double complex Ylm_l11m_4(const double theta, const double phi);
static double complex Ylm_l11m_5(const double theta, const double phi);
static double complex Ylm_l11m_6(const double theta, const double phi);
static double complex Ylm_l11m_7(const double theta, const double phi);
static double complex Ylm_l11m_8(const double theta, const double phi);
static double complex Ylm_l11m_9(const double theta, const double phi);
static double complex Ylm_l11m_10(const double theta, const double phi);
static double complex Ylm_l11m_11(const double theta, const double phi);
static double complex Ylm_l12m_1(const double theta, const double phi);
static double complex Ylm_l12m_2(const double theta, const double phi);
static double complex Ylm_l12m_3(const double theta, const double phi);
static double complex Ylm_l12m_4(const double theta, const double phi);
static double complex Ylm_l12m_5(const double theta, const double phi);
static double complex Ylm_l12m_6(const double theta, const double phi);
static double complex Ylm_l12m_7(const double theta, const double phi);
static double complex Ylm_l12m_8(const double theta, const double phi);
static double complex Ylm_l12m_9(const double theta, const double phi);
static double complex Ylm_l12m_10(const double theta, const double phi);
static double complex Ylm_l12m_11(const double theta, const double phi);
static double complex Ylm_l12m_12(const double theta, const double phi);
static double complex Ylm_l13m_1(const double theta, const double phi);
static double complex Ylm_l13m_2(const double theta, const double phi);
static double complex Ylm_l13m_3(const double theta, const double phi);
static double complex Ylm_l13m_4(const double theta, const double phi);
static double complex Ylm_l13m_5(const double theta, const double phi);
static double complex Ylm_l13m_6(const double theta, const double phi);
static double complex Ylm_l13m_7(const double theta, const double phi);
static double complex Ylm_l13m_8(const double theta, const double phi);
static double complex Ylm_l13m_9(const double theta, const double phi);
static double complex Ylm_l13m_10(const double theta, const double phi);
static double complex Ylm_l13m_11(const double theta, const double phi);
static double complex Ylm_l13m_12(const double theta, const double phi);
static double complex Ylm_l13m_13(const double theta, const double phi);
static double complex Ylm_l14m_1(const double theta, const double phi);
static double complex Ylm_l14m_2(const double theta, const double phi);
static double complex Ylm_l14m_3(const double theta, const double phi);
static double complex Ylm_l14m_4(const double theta, const double phi);
static double complex Ylm_l14m_5(const double theta, const double phi);
static double complex Ylm_l14m_6(const double theta, const double phi);
static double complex Ylm_l14m_7(const double theta, const double phi);
static double complex Ylm_l14m_8(const double theta, const double phi);
static double complex Ylm_l14m_9(const double theta, const double phi);
static double complex Ylm_l14m_10(const double theta, const double phi);
static double complex Ylm_l14m_11(const double theta, const double phi);
static double complex Ylm_l14m_12(const double theta, const double phi);
static double complex Ylm_l14m_13(const double theta, const double phi);
static double complex Ylm_l14m_14(const double theta, const double phi);
static double complex dYlm_dphi_l0m0(const double theta, const double phi);
static double complex dYlm_dphi_l1m0(const double theta, const double phi);
static double complex dYlm_dphi_l1m1(const double theta, const double phi);
static double complex dYlm_dphi_l2m0(const double theta, const double phi);
static double complex dYlm_dphi_l2m1(const double theta, const double phi);
static double complex dYlm_dphi_l2m2(const double theta, const double phi);
static double complex dYlm_dphi_l3m0(const double theta, const double phi);
static double complex dYlm_dphi_l3m1(const double theta, const double phi);
static double complex dYlm_dphi_l3m2(const double theta, const double phi);
static double complex dYlm_dphi_l3m3(const double theta, const double phi);
static double complex dYlm_dphi_l4m0(const double theta, const double phi);
static double complex dYlm_dphi_l4m1(const double theta, const double phi);
static double complex dYlm_dphi_l4m2(const double theta, const double phi);
static double complex dYlm_dphi_l4m3(const double theta, const double phi);
static double complex dYlm_dphi_l4m4(const double theta, const double phi);
static double complex dYlm_dphi_l5m0(const double theta, const double phi);
static double complex dYlm_dphi_l5m1(const double theta, const double phi);
static double complex dYlm_dphi_l5m2(const double theta, const double phi);
static double complex dYlm_dphi_l5m3(const double theta, const double phi);
static double complex dYlm_dphi_l5m4(const double theta, const double phi);
static double complex dYlm_dphi_l5m5(const double theta, const double phi);
static double complex dYlm_dphi_l6m0(const double theta, const double phi);
static double complex dYlm_dphi_l6m1(const double theta, const double phi);
static double complex dYlm_dphi_l6m2(const double theta, const double phi);
static double complex dYlm_dphi_l6m3(const double theta, const double phi);
static double complex dYlm_dphi_l6m4(const double theta, const double phi);
static double complex dYlm_dphi_l6m5(const double theta, const double phi);
static double complex dYlm_dphi_l6m6(const double theta, const double phi);
static double complex dYlm_dphi_l7m0(const double theta, const double phi);
static double complex dYlm_dphi_l7m1(const double theta, const double phi);
static double complex dYlm_dphi_l7m2(const double theta, const double phi);
static double complex dYlm_dphi_l7m3(const double theta, const double phi);
static double complex dYlm_dphi_l7m4(const double theta, const double phi);
static double complex dYlm_dphi_l7m5(const double theta, const double phi);
static double complex dYlm_dphi_l7m6(const double theta, const double phi);
static double complex dYlm_dphi_l7m7(const double theta, const double phi);
static double complex dYlm_dphi_l8m0(const double theta, const double phi);
static double complex dYlm_dphi_l8m1(const double theta, const double phi);
static double complex dYlm_dphi_l8m2(const double theta, const double phi);
static double complex dYlm_dphi_l8m3(const double theta, const double phi);
static double complex dYlm_dphi_l8m4(const double theta, const double phi);
static double complex dYlm_dphi_l8m5(const double theta, const double phi);
static double complex dYlm_dphi_l8m6(const double theta, const double phi);
static double complex dYlm_dphi_l8m7(const double theta, const double phi);
static double complex dYlm_dphi_l8m8(const double theta, const double phi);
static double complex dYlm_dphi_l9m0(const double theta, const double phi);
static double complex dYlm_dphi_l9m1(const double theta, const double phi);
static double complex dYlm_dphi_l9m2(const double theta, const double phi);
static double complex dYlm_dphi_l9m3(const double theta, const double phi);
static double complex dYlm_dphi_l9m4(const double theta, const double phi);
static double complex dYlm_dphi_l9m5(const double theta, const double phi);
static double complex dYlm_dphi_l9m6(const double theta, const double phi);
static double complex dYlm_dphi_l9m7(const double theta, const double phi);
static double complex dYlm_dphi_l9m8(const double theta, const double phi);
static double complex dYlm_dphi_l9m9(const double theta, const double phi);
static double complex dYlm_dphi_l10m0(const double theta, const double phi);
static double complex dYlm_dphi_l10m1(const double theta, const double phi);
static double complex dYlm_dphi_l10m2(const double theta, const double phi);
static double complex dYlm_dphi_l10m3(const double theta, const double phi);
static double complex dYlm_dphi_l10m4(const double theta, const double phi);
static double complex dYlm_dphi_l10m5(const double theta, const double phi);
static double complex dYlm_dphi_l10m6(const double theta, const double phi);
static double complex dYlm_dphi_l10m7(const double theta, const double phi);
static double complex dYlm_dphi_l10m8(const double theta, const double phi);
static double complex dYlm_dphi_l10m9(const double theta, const double phi);
static double complex dYlm_dphi_l10m10(const double theta, const double phi);
static double complex dYlm_dphi_l11m0(const double theta, const double phi);
static double complex dYlm_dphi_l11m1(const double theta, const double phi);
static double complex dYlm_dphi_l11m2(const double theta, const double phi);
static double complex dYlm_dphi_l11m3(const double theta, const double phi);
static double complex dYlm_dphi_l11m4(const double theta, const double phi);
static double complex dYlm_dphi_l11m5(const double theta, const double phi);
static double complex dYlm_dphi_l11m6(const double theta, const double phi);
static double complex dYlm_dphi_l11m7(const double theta, const double phi);
static double complex dYlm_dphi_l11m8(const double theta, const double phi);
static double complex dYlm_dphi_l11m9(const double theta, const double phi);
static double complex dYlm_dphi_l11m10(const double theta, const double phi);
static double complex dYlm_dphi_l11m11(const double theta, const double phi);
static double complex dYlm_dphi_l12m0(const double theta, const double phi);
static double complex dYlm_dphi_l12m1(const double theta, const double phi);
static double complex dYlm_dphi_l12m2(const double theta, const double phi);
static double complex dYlm_dphi_l12m3(const double theta, const double phi);
static double complex dYlm_dphi_l12m4(const double theta, const double phi);
static double complex dYlm_dphi_l12m5(const double theta, const double phi);
static double complex dYlm_dphi_l12m6(const double theta, const double phi);
static double complex dYlm_dphi_l12m7(const double theta, const double phi);
static double complex dYlm_dphi_l12m8(const double theta, const double phi);
static double complex dYlm_dphi_l12m9(const double theta, const double phi);
static double complex dYlm_dphi_l12m10(const double theta, const double phi);
static double complex dYlm_dphi_l12m11(const double theta, const double phi);
static double complex dYlm_dphi_l12m12(const double theta, const double phi);
static double complex dYlm_dphi_l13m0(const double theta, const double phi);
static double complex dYlm_dphi_l13m1(const double theta, const double phi);
static double complex dYlm_dphi_l13m2(const double theta, const double phi);
static double complex dYlm_dphi_l13m3(const double theta, const double phi);
static double complex dYlm_dphi_l13m4(const double theta, const double phi);
static double complex dYlm_dphi_l13m5(const double theta, const double phi);
static double complex dYlm_dphi_l13m6(const double theta, const double phi);
static double complex dYlm_dphi_l13m7(const double theta, const double phi);
static double complex dYlm_dphi_l13m8(const double theta, const double phi);
static double complex dYlm_dphi_l13m9(const double theta, const double phi);
static double complex dYlm_dphi_l13m10(const double theta, const double phi);
static double complex dYlm_dphi_l13m11(const double theta, const double phi);
static double complex dYlm_dphi_l13m12(const double theta, const double phi);
static double complex dYlm_dphi_l13m13(const double theta, const double phi);
static double complex dYlm_dphi_l14m0(const double theta, const double phi);
static double complex dYlm_dphi_l14m1(const double theta, const double phi);
static double complex dYlm_dphi_l14m2(const double theta, const double phi);
static double complex dYlm_dphi_l14m3(const double theta, const double phi);
static double complex dYlm_dphi_l14m4(const double theta, const double phi);
static double complex dYlm_dphi_l14m5(const double theta, const double phi);
static double complex dYlm_dphi_l14m6(const double theta, const double phi);
static double complex dYlm_dphi_l14m7(const double theta, const double phi);
static double complex dYlm_dphi_l14m8(const double theta, const double phi);
static double complex dYlm_dphi_l14m9(const double theta, const double phi);
static double complex dYlm_dphi_l14m10(const double theta, const double phi);
static double complex dYlm_dphi_l14m11(const double theta, const double phi);
static double complex dYlm_dphi_l14m12(const double theta, const double phi);
static double complex dYlm_dphi_l14m13(const double theta, const double phi);
static double complex dYlm_dphi_l14m14(const double theta, const double phi);
static double complex dYlm_dphi_l1m_1(const double theta, const double phi);
static double complex dYlm_dphi_l2m_1(const double theta, const double phi);
static double complex dYlm_dphi_l2m_2(const double theta, const double phi);
static double complex dYlm_dphi_l3m_1(const double theta, const double phi);
static double complex dYlm_dphi_l3m_2(const double theta, const double phi);
static double complex dYlm_dphi_l3m_3(const double theta, const double phi);
static double complex dYlm_dphi_l4m_1(const double theta, const double phi);
static double complex dYlm_dphi_l4m_2(const double theta, const double phi);
static double complex dYlm_dphi_l4m_3(const double theta, const double phi);
static double complex dYlm_dphi_l4m_4(const double theta, const double phi);
static double complex dYlm_dphi_l5m_1(const double theta, const double phi);
static double complex dYlm_dphi_l5m_2(const double theta, const double phi);
static double complex dYlm_dphi_l5m_3(const double theta, const double phi);
static double complex dYlm_dphi_l5m_4(const double theta, const double phi);
static double complex dYlm_dphi_l5m_5(const double theta, const double phi);
static double complex dYlm_dphi_l6m_1(const double theta, const double phi);
static double complex dYlm_dphi_l6m_2(const double theta, const double phi);
static double complex dYlm_dphi_l6m_3(const double theta, const double phi);
static double complex dYlm_dphi_l6m_4(const double theta, const double phi);
static double complex dYlm_dphi_l6m_5(const double theta, const double phi);
static double complex dYlm_dphi_l6m_6(const double theta, const double phi);
static double complex dYlm_dphi_l7m_1(const double theta, const double phi);
static double complex dYlm_dphi_l7m_2(const double theta, const double phi);
static double complex dYlm_dphi_l7m_3(const double theta, const double phi);
static double complex dYlm_dphi_l7m_4(const double theta, const double phi);
static double complex dYlm_dphi_l7m_5(const double theta, const double phi);
static double complex dYlm_dphi_l7m_6(const double theta, const double phi);
static double complex dYlm_dphi_l7m_7(const double theta, const double phi);
static double complex dYlm_dphi_l8m_1(const double theta, const double phi);
static double complex dYlm_dphi_l8m_2(const double theta, const double phi);
static double complex dYlm_dphi_l8m_3(const double theta, const double phi);
static double complex dYlm_dphi_l8m_4(const double theta, const double phi);
static double complex dYlm_dphi_l8m_5(const double theta, const double phi);
static double complex dYlm_dphi_l8m_6(const double theta, const double phi);
static double complex dYlm_dphi_l8m_7(const double theta, const double phi);
static double complex dYlm_dphi_l8m_8(const double theta, const double phi);
static double complex dYlm_dphi_l9m_1(const double theta, const double phi);
static double complex dYlm_dphi_l9m_2(const double theta, const double phi);
static double complex dYlm_dphi_l9m_3(const double theta, const double phi);
static double complex dYlm_dphi_l9m_4(const double theta, const double phi);
static double complex dYlm_dphi_l9m_5(const double theta, const double phi);
static double complex dYlm_dphi_l9m_6(const double theta, const double phi);
static double complex dYlm_dphi_l9m_7(const double theta, const double phi);
static double complex dYlm_dphi_l9m_8(const double theta, const double phi);
static double complex dYlm_dphi_l9m_9(const double theta, const double phi);
static double complex dYlm_dphi_l10m_1(const double theta, const double phi);
static double complex dYlm_dphi_l10m_2(const double theta, const double phi);
static double complex dYlm_dphi_l10m_3(const double theta, const double phi);
static double complex dYlm_dphi_l10m_4(const double theta, const double phi);
static double complex dYlm_dphi_l10m_5(const double theta, const double phi);
static double complex dYlm_dphi_l10m_6(const double theta, const double phi);
static double complex dYlm_dphi_l10m_7(const double theta, const double phi);
static double complex dYlm_dphi_l10m_8(const double theta, const double phi);
static double complex dYlm_dphi_l10m_9(const double theta, const double phi);
static double complex dYlm_dphi_l10m_10(const double theta, const double phi);
static double complex dYlm_dphi_l11m_1(const double theta, const double phi);
static double complex dYlm_dphi_l11m_2(const double theta, const double phi);
static double complex dYlm_dphi_l11m_3(const double theta, const double phi);
static double complex dYlm_dphi_l11m_4(const double theta, const double phi);
static double complex dYlm_dphi_l11m_5(const double theta, const double phi);
static double complex dYlm_dphi_l11m_6(const double theta, const double phi);
static double complex dYlm_dphi_l11m_7(const double theta, const double phi);
static double complex dYlm_dphi_l11m_8(const double theta, const double phi);
static double complex dYlm_dphi_l11m_9(const double theta, const double phi);
static double complex dYlm_dphi_l11m_10(const double theta, const double phi);
static double complex dYlm_dphi_l11m_11(const double theta, const double phi);
static double complex dYlm_dphi_l12m_1(const double theta, const double phi);
static double complex dYlm_dphi_l12m_2(const double theta, const double phi);
static double complex dYlm_dphi_l12m_3(const double theta, const double phi);
static double complex dYlm_dphi_l12m_4(const double theta, const double phi);
static double complex dYlm_dphi_l12m_5(const double theta, const double phi);
static double complex dYlm_dphi_l12m_6(const double theta, const double phi);
static double complex dYlm_dphi_l12m_7(const double theta, const double phi);
static double complex dYlm_dphi_l12m_8(const double theta, const double phi);
static double complex dYlm_dphi_l12m_9(const double theta, const double phi);
static double complex dYlm_dphi_l12m_10(const double theta, const double phi);
static double complex dYlm_dphi_l12m_11(const double theta, const double phi);
static double complex dYlm_dphi_l12m_12(const double theta, const double phi);
static double complex dYlm_dphi_l13m_1(const double theta, const double phi);
static double complex dYlm_dphi_l13m_2(const double theta, const double phi);
static double complex dYlm_dphi_l13m_3(const double theta, const double phi);
static double complex dYlm_dphi_l13m_4(const double theta, const double phi);
static double complex dYlm_dphi_l13m_5(const double theta, const double phi);
static double complex dYlm_dphi_l13m_6(const double theta, const double phi);
static double complex dYlm_dphi_l13m_7(const double theta, const double phi);
static double complex dYlm_dphi_l13m_8(const double theta, const double phi);
static double complex dYlm_dphi_l13m_9(const double theta, const double phi);
static double complex dYlm_dphi_l13m_10(const double theta, const double phi);
static double complex dYlm_dphi_l13m_11(const double theta, const double phi);
static double complex dYlm_dphi_l13m_12(const double theta, const double phi);
static double complex dYlm_dphi_l13m_13(const double theta, const double phi);
static double complex dYlm_dphi_l14m_1(const double theta, const double phi);
static double complex dYlm_dphi_l14m_2(const double theta, const double phi);
static double complex dYlm_dphi_l14m_3(const double theta, const double phi);
static double complex dYlm_dphi_l14m_4(const double theta, const double phi);
static double complex dYlm_dphi_l14m_5(const double theta, const double phi);
static double complex dYlm_dphi_l14m_6(const double theta, const double phi);
static double complex dYlm_dphi_l14m_7(const double theta, const double phi);
static double complex dYlm_dphi_l14m_8(const double theta, const double phi);
static double complex dYlm_dphi_l14m_9(const double theta, const double phi);
static double complex dYlm_dphi_l14m_10(const double theta, const double phi);
static double complex dYlm_dphi_l14m_11(const double theta, const double phi);
static double complex dYlm_dphi_l14m_12(const double theta, const double phi);
static double complex dYlm_dphi_l14m_13(const double theta, const double phi);
static double complex dYlm_dphi_l14m_14(const double theta, const double phi);
static double complex dYlm_dtheta_l0m0(const double theta, const double phi);
static double complex dYlm_dtheta_l1m0(const double theta, const double phi);
static double complex dYlm_dtheta_l1m1(const double theta, const double phi);
static double complex dYlm_dtheta_l2m0(const double theta, const double phi);
static double complex dYlm_dtheta_l2m1(const double theta, const double phi);
static double complex dYlm_dtheta_l2m2(const double theta, const double phi);
static double complex dYlm_dtheta_l3m0(const double theta, const double phi);
static double complex dYlm_dtheta_l3m1(const double theta, const double phi);
static double complex dYlm_dtheta_l3m2(const double theta, const double phi);
static double complex dYlm_dtheta_l3m3(const double theta, const double phi);
static double complex dYlm_dtheta_l4m0(const double theta, const double phi);
static double complex dYlm_dtheta_l4m1(const double theta, const double phi);
static double complex dYlm_dtheta_l4m2(const double theta, const double phi);
static double complex dYlm_dtheta_l4m3(const double theta, const double phi);
static double complex dYlm_dtheta_l4m4(const double theta, const double phi);
static double complex dYlm_dtheta_l5m0(const double theta, const double phi);
static double complex dYlm_dtheta_l5m1(const double theta, const double phi);
static double complex dYlm_dtheta_l5m2(const double theta, const double phi);
static double complex dYlm_dtheta_l5m3(const double theta, const double phi);
static double complex dYlm_dtheta_l5m4(const double theta, const double phi);
static double complex dYlm_dtheta_l5m5(const double theta, const double phi);
static double complex dYlm_dtheta_l6m0(const double theta, const double phi);
static double complex dYlm_dtheta_l6m1(const double theta, const double phi);
static double complex dYlm_dtheta_l6m2(const double theta, const double phi);
static double complex dYlm_dtheta_l6m3(const double theta, const double phi);
static double complex dYlm_dtheta_l6m4(const double theta, const double phi);
static double complex dYlm_dtheta_l6m5(const double theta, const double phi);
static double complex dYlm_dtheta_l6m6(const double theta, const double phi);
static double complex dYlm_dtheta_l7m0(const double theta, const double phi);
static double complex dYlm_dtheta_l7m1(const double theta, const double phi);
static double complex dYlm_dtheta_l7m2(const double theta, const double phi);
static double complex dYlm_dtheta_l7m3(const double theta, const double phi);
static double complex dYlm_dtheta_l7m4(const double theta, const double phi);
static double complex dYlm_dtheta_l7m5(const double theta, const double phi);
static double complex dYlm_dtheta_l7m6(const double theta, const double phi);
static double complex dYlm_dtheta_l7m7(const double theta, const double phi);
static double complex dYlm_dtheta_l8m0(const double theta, const double phi);
static double complex dYlm_dtheta_l8m1(const double theta, const double phi);
static double complex dYlm_dtheta_l8m2(const double theta, const double phi);
static double complex dYlm_dtheta_l8m3(const double theta, const double phi);
static double complex dYlm_dtheta_l8m4(const double theta, const double phi);
static double complex dYlm_dtheta_l8m5(const double theta, const double phi);
static double complex dYlm_dtheta_l8m6(const double theta, const double phi);
static double complex dYlm_dtheta_l8m7(const double theta, const double phi);
static double complex dYlm_dtheta_l8m8(const double theta, const double phi);
static double complex dYlm_dtheta_l9m0(const double theta, const double phi);
static double complex dYlm_dtheta_l9m1(const double theta, const double phi);
static double complex dYlm_dtheta_l9m2(const double theta, const double phi);
static double complex dYlm_dtheta_l9m3(const double theta, const double phi);
static double complex dYlm_dtheta_l9m4(const double theta, const double phi);
static double complex dYlm_dtheta_l9m5(const double theta, const double phi);
static double complex dYlm_dtheta_l9m6(const double theta, const double phi);
static double complex dYlm_dtheta_l9m7(const double theta, const double phi);
static double complex dYlm_dtheta_l9m8(const double theta, const double phi);
static double complex dYlm_dtheta_l9m9(const double theta, const double phi);
static double complex dYlm_dtheta_l10m0(const double theta, const double phi);
static double complex dYlm_dtheta_l10m1(const double theta, const double phi);
static double complex dYlm_dtheta_l10m2(const double theta, const double phi);
static double complex dYlm_dtheta_l10m3(const double theta, const double phi);
static double complex dYlm_dtheta_l10m4(const double theta, const double phi);
static double complex dYlm_dtheta_l10m5(const double theta, const double phi);
static double complex dYlm_dtheta_l10m6(const double theta, const double phi);
static double complex dYlm_dtheta_l10m7(const double theta, const double phi);
static double complex dYlm_dtheta_l10m8(const double theta, const double phi);
static double complex dYlm_dtheta_l10m9(const double theta, const double phi);
static double complex dYlm_dtheta_l10m10(const double theta, const double phi);
static double complex dYlm_dtheta_l11m0(const double theta, const double phi);
static double complex dYlm_dtheta_l11m1(const double theta, const double phi);
static double complex dYlm_dtheta_l11m2(const double theta, const double phi);
static double complex dYlm_dtheta_l11m3(const double theta, const double phi);
static double complex dYlm_dtheta_l11m4(const double theta, const double phi);
static double complex dYlm_dtheta_l11m5(const double theta, const double phi);
static double complex dYlm_dtheta_l11m6(const double theta, const double phi);
static double complex dYlm_dtheta_l11m7(const double theta, const double phi);
static double complex dYlm_dtheta_l11m8(const double theta, const double phi);
static double complex dYlm_dtheta_l11m9(const double theta, const double phi);
static double complex dYlm_dtheta_l11m10(const double theta, const double phi);
static double complex dYlm_dtheta_l11m11(const double theta, const double phi);
static double complex dYlm_dtheta_l12m0(const double theta, const double phi);
static double complex dYlm_dtheta_l12m1(const double theta, const double phi);
static double complex dYlm_dtheta_l12m2(const double theta, const double phi);
static double complex dYlm_dtheta_l12m3(const double theta, const double phi);
static double complex dYlm_dtheta_l12m4(const double theta, const double phi);
static double complex dYlm_dtheta_l12m5(const double theta, const double phi);
static double complex dYlm_dtheta_l12m6(const double theta, const double phi);
static double complex dYlm_dtheta_l12m7(const double theta, const double phi);
static double complex dYlm_dtheta_l12m8(const double theta, const double phi);
static double complex dYlm_dtheta_l12m9(const double theta, const double phi);
static double complex dYlm_dtheta_l12m10(const double theta, const double phi);
static double complex dYlm_dtheta_l12m11(const double theta, const double phi);
static double complex dYlm_dtheta_l12m12(const double theta, const double phi);
static double complex dYlm_dtheta_l13m0(const double theta, const double phi);
static double complex dYlm_dtheta_l13m1(const double theta, const double phi);
static double complex dYlm_dtheta_l13m2(const double theta, const double phi);
static double complex dYlm_dtheta_l13m3(const double theta, const double phi);
static double complex dYlm_dtheta_l13m4(const double theta, const double phi);
static double complex dYlm_dtheta_l13m5(const double theta, const double phi);
static double complex dYlm_dtheta_l13m6(const double theta, const double phi);
static double complex dYlm_dtheta_l13m7(const double theta, const double phi);
static double complex dYlm_dtheta_l13m8(const double theta, const double phi);
static double complex dYlm_dtheta_l13m9(const double theta, const double phi);
static double complex dYlm_dtheta_l13m10(const double theta, const double phi);
static double complex dYlm_dtheta_l13m11(const double theta, const double phi);
static double complex dYlm_dtheta_l13m12(const double theta, const double phi);
static double complex dYlm_dtheta_l13m13(const double theta, const double phi);
static double complex dYlm_dtheta_l14m0(const double theta, const double phi);
static double complex dYlm_dtheta_l14m1(const double theta, const double phi);
static double complex dYlm_dtheta_l14m2(const double theta, const double phi);
static double complex dYlm_dtheta_l14m3(const double theta, const double phi);
static double complex dYlm_dtheta_l14m4(const double theta, const double phi);
static double complex dYlm_dtheta_l14m5(const double theta, const double phi);
static double complex dYlm_dtheta_l14m6(const double theta, const double phi);
static double complex dYlm_dtheta_l14m7(const double theta, const double phi);
static double complex dYlm_dtheta_l14m8(const double theta, const double phi);
static double complex dYlm_dtheta_l14m9(const double theta, const double phi);
static double complex dYlm_dtheta_l14m10(const double theta, const double phi);
static double complex dYlm_dtheta_l14m11(const double theta, const double phi);
static double complex dYlm_dtheta_l14m12(const double theta, const double phi);
static double complex dYlm_dtheta_l14m13(const double theta, const double phi);
static double complex dYlm_dtheta_l14m14(const double theta, const double phi);
static double complex dYlm_dtheta_l1m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l2m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l2m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l3m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l3m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l3m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l4m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l4m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l4m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l4m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l5m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l5m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l5m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l5m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l5m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l6m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l6m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l6m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l6m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l6m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l6m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l7m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l8m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l9m_9(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_9(const double theta, const double phi);
static double complex dYlm_dtheta_l10m_10(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_9(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_10(const double theta, const double phi);
static double complex dYlm_dtheta_l11m_11(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_9(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_10(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_11(const double theta, const double phi);
static double complex dYlm_dtheta_l12m_12(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_9(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_10(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_11(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_12(const double theta, const double phi);
static double complex dYlm_dtheta_l13m_13(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_1(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_2(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_3(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_4(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_5(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_6(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_7(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_8(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_9(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_10(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_11(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_12(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_13(const double theta, const double phi);
static double complex dYlm_dtheta_l14m_14(const double theta, const double phi);



/* initializing Y_{l}^{m}(x) */
void init_Ylm(void)
{
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
}
/* Y_n^m(\theta, \varphi) := \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}} exp(i m \varphi)\mathrm{P}_n^m\left(\cos(\theta)\right) 
// ->return value: Y_{l}^{m}(x) */
double complex Ylm(const unsigned l, int m, const double theta, const double phi)
{
  if (theta > M_PI || theta < 0)
    Error0("theta exceeds from [0,pi] interval.\n");
  if (phi > 2*M_PI || phi < 0)
    Error0("phi exceeds from [0,2*pi] interval.\n");
  if (l >= 15)
    Error0("l exceeds the maximum.\n");
  if ((unsigned)abs(m) > l)
    return 0;
  if (m < 0)
    return Y_[l][-m](theta,phi);
  return Y[l][m](theta,phi);
}



/* initializing dYlm_dphi */
void init_dYlm_dphi(void)
{
  dY_dphi[0][0] = dYlm_dphi_l0m0;
  dY_dphi[1][0] = dYlm_dphi_l1m0;
  dY_dphi[1][1] = dYlm_dphi_l1m1;
  dY_dphi[2][0] = dYlm_dphi_l2m0;
  dY_dphi[2][1] = dYlm_dphi_l2m1;
  dY_dphi[2][2] = dYlm_dphi_l2m2;
  dY_dphi[3][0] = dYlm_dphi_l3m0;
  dY_dphi[3][1] = dYlm_dphi_l3m1;
  dY_dphi[3][2] = dYlm_dphi_l3m2;
  dY_dphi[3][3] = dYlm_dphi_l3m3;
  dY_dphi[4][0] = dYlm_dphi_l4m0;
  dY_dphi[4][1] = dYlm_dphi_l4m1;
  dY_dphi[4][2] = dYlm_dphi_l4m2;
  dY_dphi[4][3] = dYlm_dphi_l4m3;
  dY_dphi[4][4] = dYlm_dphi_l4m4;
  dY_dphi[5][0] = dYlm_dphi_l5m0;
  dY_dphi[5][1] = dYlm_dphi_l5m1;
  dY_dphi[5][2] = dYlm_dphi_l5m2;
  dY_dphi[5][3] = dYlm_dphi_l5m3;
  dY_dphi[5][4] = dYlm_dphi_l5m4;
  dY_dphi[5][5] = dYlm_dphi_l5m5;
  dY_dphi[6][0] = dYlm_dphi_l6m0;
  dY_dphi[6][1] = dYlm_dphi_l6m1;
  dY_dphi[6][2] = dYlm_dphi_l6m2;
  dY_dphi[6][3] = dYlm_dphi_l6m3;
  dY_dphi[6][4] = dYlm_dphi_l6m4;
  dY_dphi[6][5] = dYlm_dphi_l6m5;
  dY_dphi[6][6] = dYlm_dphi_l6m6;
  dY_dphi[7][0] = dYlm_dphi_l7m0;
  dY_dphi[7][1] = dYlm_dphi_l7m1;
  dY_dphi[7][2] = dYlm_dphi_l7m2;
  dY_dphi[7][3] = dYlm_dphi_l7m3;
  dY_dphi[7][4] = dYlm_dphi_l7m4;
  dY_dphi[7][5] = dYlm_dphi_l7m5;
  dY_dphi[7][6] = dYlm_dphi_l7m6;
  dY_dphi[7][7] = dYlm_dphi_l7m7;
  dY_dphi[8][0] = dYlm_dphi_l8m0;
  dY_dphi[8][1] = dYlm_dphi_l8m1;
  dY_dphi[8][2] = dYlm_dphi_l8m2;
  dY_dphi[8][3] = dYlm_dphi_l8m3;
  dY_dphi[8][4] = dYlm_dphi_l8m4;
  dY_dphi[8][5] = dYlm_dphi_l8m5;
  dY_dphi[8][6] = dYlm_dphi_l8m6;
  dY_dphi[8][7] = dYlm_dphi_l8m7;
  dY_dphi[8][8] = dYlm_dphi_l8m8;
  dY_dphi[9][0] = dYlm_dphi_l9m0;
  dY_dphi[9][1] = dYlm_dphi_l9m1;
  dY_dphi[9][2] = dYlm_dphi_l9m2;
  dY_dphi[9][3] = dYlm_dphi_l9m3;
  dY_dphi[9][4] = dYlm_dphi_l9m4;
  dY_dphi[9][5] = dYlm_dphi_l9m5;
  dY_dphi[9][6] = dYlm_dphi_l9m6;
  dY_dphi[9][7] = dYlm_dphi_l9m7;
  dY_dphi[9][8] = dYlm_dphi_l9m8;
  dY_dphi[9][9] = dYlm_dphi_l9m9;
  dY_dphi[10][0] = dYlm_dphi_l10m0;
  dY_dphi[10][1] = dYlm_dphi_l10m1;
  dY_dphi[10][2] = dYlm_dphi_l10m2;
  dY_dphi[10][3] = dYlm_dphi_l10m3;
  dY_dphi[10][4] = dYlm_dphi_l10m4;
  dY_dphi[10][5] = dYlm_dphi_l10m5;
  dY_dphi[10][6] = dYlm_dphi_l10m6;
  dY_dphi[10][7] = dYlm_dphi_l10m7;
  dY_dphi[10][8] = dYlm_dphi_l10m8;
  dY_dphi[10][9] = dYlm_dphi_l10m9;
  dY_dphi[10][10] = dYlm_dphi_l10m10;
  dY_dphi[11][0] = dYlm_dphi_l11m0;
  dY_dphi[11][1] = dYlm_dphi_l11m1;
  dY_dphi[11][2] = dYlm_dphi_l11m2;
  dY_dphi[11][3] = dYlm_dphi_l11m3;
  dY_dphi[11][4] = dYlm_dphi_l11m4;
  dY_dphi[11][5] = dYlm_dphi_l11m5;
  dY_dphi[11][6] = dYlm_dphi_l11m6;
  dY_dphi[11][7] = dYlm_dphi_l11m7;
  dY_dphi[11][8] = dYlm_dphi_l11m8;
  dY_dphi[11][9] = dYlm_dphi_l11m9;
  dY_dphi[11][10] = dYlm_dphi_l11m10;
  dY_dphi[11][11] = dYlm_dphi_l11m11;
  dY_dphi[12][0] = dYlm_dphi_l12m0;
  dY_dphi[12][1] = dYlm_dphi_l12m1;
  dY_dphi[12][2] = dYlm_dphi_l12m2;
  dY_dphi[12][3] = dYlm_dphi_l12m3;
  dY_dphi[12][4] = dYlm_dphi_l12m4;
  dY_dphi[12][5] = dYlm_dphi_l12m5;
  dY_dphi[12][6] = dYlm_dphi_l12m6;
  dY_dphi[12][7] = dYlm_dphi_l12m7;
  dY_dphi[12][8] = dYlm_dphi_l12m8;
  dY_dphi[12][9] = dYlm_dphi_l12m9;
  dY_dphi[12][10] = dYlm_dphi_l12m10;
  dY_dphi[12][11] = dYlm_dphi_l12m11;
  dY_dphi[12][12] = dYlm_dphi_l12m12;
  dY_dphi[13][0] = dYlm_dphi_l13m0;
  dY_dphi[13][1] = dYlm_dphi_l13m1;
  dY_dphi[13][2] = dYlm_dphi_l13m2;
  dY_dphi[13][3] = dYlm_dphi_l13m3;
  dY_dphi[13][4] = dYlm_dphi_l13m4;
  dY_dphi[13][5] = dYlm_dphi_l13m5;
  dY_dphi[13][6] = dYlm_dphi_l13m6;
  dY_dphi[13][7] = dYlm_dphi_l13m7;
  dY_dphi[13][8] = dYlm_dphi_l13m8;
  dY_dphi[13][9] = dYlm_dphi_l13m9;
  dY_dphi[13][10] = dYlm_dphi_l13m10;
  dY_dphi[13][11] = dYlm_dphi_l13m11;
  dY_dphi[13][12] = dYlm_dphi_l13m12;
  dY_dphi[13][13] = dYlm_dphi_l13m13;
  dY_dphi[14][0] = dYlm_dphi_l14m0;
  dY_dphi[14][1] = dYlm_dphi_l14m1;
  dY_dphi[14][2] = dYlm_dphi_l14m2;
  dY_dphi[14][3] = dYlm_dphi_l14m3;
  dY_dphi[14][4] = dYlm_dphi_l14m4;
  dY_dphi[14][5] = dYlm_dphi_l14m5;
  dY_dphi[14][6] = dYlm_dphi_l14m6;
  dY_dphi[14][7] = dYlm_dphi_l14m7;
  dY_dphi[14][8] = dYlm_dphi_l14m8;
  dY_dphi[14][9] = dYlm_dphi_l14m9;
  dY_dphi[14][10] = dYlm_dphi_l14m10;
  dY_dphi[14][11] = dYlm_dphi_l14m11;
  dY_dphi[14][12] = dYlm_dphi_l14m12;
  dY_dphi[14][13] = dYlm_dphi_l14m13;
  dY_dphi[14][14] = dYlm_dphi_l14m14;
  dY_dphi_[1][1] = dYlm_dphi_l1m_1;
  dY_dphi_[2][1] = dYlm_dphi_l2m_1;
  dY_dphi_[2][2] = dYlm_dphi_l2m_2;
  dY_dphi_[3][1] = dYlm_dphi_l3m_1;
  dY_dphi_[3][2] = dYlm_dphi_l3m_2;
  dY_dphi_[3][3] = dYlm_dphi_l3m_3;
  dY_dphi_[4][1] = dYlm_dphi_l4m_1;
  dY_dphi_[4][2] = dYlm_dphi_l4m_2;
  dY_dphi_[4][3] = dYlm_dphi_l4m_3;
  dY_dphi_[4][4] = dYlm_dphi_l4m_4;
  dY_dphi_[5][1] = dYlm_dphi_l5m_1;
  dY_dphi_[5][2] = dYlm_dphi_l5m_2;
  dY_dphi_[5][3] = dYlm_dphi_l5m_3;
  dY_dphi_[5][4] = dYlm_dphi_l5m_4;
  dY_dphi_[5][5] = dYlm_dphi_l5m_5;
  dY_dphi_[6][1] = dYlm_dphi_l6m_1;
  dY_dphi_[6][2] = dYlm_dphi_l6m_2;
  dY_dphi_[6][3] = dYlm_dphi_l6m_3;
  dY_dphi_[6][4] = dYlm_dphi_l6m_4;
  dY_dphi_[6][5] = dYlm_dphi_l6m_5;
  dY_dphi_[6][6] = dYlm_dphi_l6m_6;
  dY_dphi_[7][1] = dYlm_dphi_l7m_1;
  dY_dphi_[7][2] = dYlm_dphi_l7m_2;
  dY_dphi_[7][3] = dYlm_dphi_l7m_3;
  dY_dphi_[7][4] = dYlm_dphi_l7m_4;
  dY_dphi_[7][5] = dYlm_dphi_l7m_5;
  dY_dphi_[7][6] = dYlm_dphi_l7m_6;
  dY_dphi_[7][7] = dYlm_dphi_l7m_7;
  dY_dphi_[8][1] = dYlm_dphi_l8m_1;
  dY_dphi_[8][2] = dYlm_dphi_l8m_2;
  dY_dphi_[8][3] = dYlm_dphi_l8m_3;
  dY_dphi_[8][4] = dYlm_dphi_l8m_4;
  dY_dphi_[8][5] = dYlm_dphi_l8m_5;
  dY_dphi_[8][6] = dYlm_dphi_l8m_6;
  dY_dphi_[8][7] = dYlm_dphi_l8m_7;
  dY_dphi_[8][8] = dYlm_dphi_l8m_8;
  dY_dphi_[9][1] = dYlm_dphi_l9m_1;
  dY_dphi_[9][2] = dYlm_dphi_l9m_2;
  dY_dphi_[9][3] = dYlm_dphi_l9m_3;
  dY_dphi_[9][4] = dYlm_dphi_l9m_4;
  dY_dphi_[9][5] = dYlm_dphi_l9m_5;
  dY_dphi_[9][6] = dYlm_dphi_l9m_6;
  dY_dphi_[9][7] = dYlm_dphi_l9m_7;
  dY_dphi_[9][8] = dYlm_dphi_l9m_8;
  dY_dphi_[9][9] = dYlm_dphi_l9m_9;
  dY_dphi_[10][1] = dYlm_dphi_l10m_1;
  dY_dphi_[10][2] = dYlm_dphi_l10m_2;
  dY_dphi_[10][3] = dYlm_dphi_l10m_3;
  dY_dphi_[10][4] = dYlm_dphi_l10m_4;
  dY_dphi_[10][5] = dYlm_dphi_l10m_5;
  dY_dphi_[10][6] = dYlm_dphi_l10m_6;
  dY_dphi_[10][7] = dYlm_dphi_l10m_7;
  dY_dphi_[10][8] = dYlm_dphi_l10m_8;
  dY_dphi_[10][9] = dYlm_dphi_l10m_9;
  dY_dphi_[10][10] = dYlm_dphi_l10m_10;
  dY_dphi_[11][1] = dYlm_dphi_l11m_1;
  dY_dphi_[11][2] = dYlm_dphi_l11m_2;
  dY_dphi_[11][3] = dYlm_dphi_l11m_3;
  dY_dphi_[11][4] = dYlm_dphi_l11m_4;
  dY_dphi_[11][5] = dYlm_dphi_l11m_5;
  dY_dphi_[11][6] = dYlm_dphi_l11m_6;
  dY_dphi_[11][7] = dYlm_dphi_l11m_7;
  dY_dphi_[11][8] = dYlm_dphi_l11m_8;
  dY_dphi_[11][9] = dYlm_dphi_l11m_9;
  dY_dphi_[11][10] = dYlm_dphi_l11m_10;
  dY_dphi_[11][11] = dYlm_dphi_l11m_11;
  dY_dphi_[12][1] = dYlm_dphi_l12m_1;
  dY_dphi_[12][2] = dYlm_dphi_l12m_2;
  dY_dphi_[12][3] = dYlm_dphi_l12m_3;
  dY_dphi_[12][4] = dYlm_dphi_l12m_4;
  dY_dphi_[12][5] = dYlm_dphi_l12m_5;
  dY_dphi_[12][6] = dYlm_dphi_l12m_6;
  dY_dphi_[12][7] = dYlm_dphi_l12m_7;
  dY_dphi_[12][8] = dYlm_dphi_l12m_8;
  dY_dphi_[12][9] = dYlm_dphi_l12m_9;
  dY_dphi_[12][10] = dYlm_dphi_l12m_10;
  dY_dphi_[12][11] = dYlm_dphi_l12m_11;
  dY_dphi_[12][12] = dYlm_dphi_l12m_12;
  dY_dphi_[13][1] = dYlm_dphi_l13m_1;
  dY_dphi_[13][2] = dYlm_dphi_l13m_2;
  dY_dphi_[13][3] = dYlm_dphi_l13m_3;
  dY_dphi_[13][4] = dYlm_dphi_l13m_4;
  dY_dphi_[13][5] = dYlm_dphi_l13m_5;
  dY_dphi_[13][6] = dYlm_dphi_l13m_6;
  dY_dphi_[13][7] = dYlm_dphi_l13m_7;
  dY_dphi_[13][8] = dYlm_dphi_l13m_8;
  dY_dphi_[13][9] = dYlm_dphi_l13m_9;
  dY_dphi_[13][10] = dYlm_dphi_l13m_10;
  dY_dphi_[13][11] = dYlm_dphi_l13m_11;
  dY_dphi_[13][12] = dYlm_dphi_l13m_12;
  dY_dphi_[13][13] = dYlm_dphi_l13m_13;
  dY_dphi_[14][1] = dYlm_dphi_l14m_1;
  dY_dphi_[14][2] = dYlm_dphi_l14m_2;
  dY_dphi_[14][3] = dYlm_dphi_l14m_3;
  dY_dphi_[14][4] = dYlm_dphi_l14m_4;
  dY_dphi_[14][5] = dYlm_dphi_l14m_5;
  dY_dphi_[14][6] = dYlm_dphi_l14m_6;
  dY_dphi_[14][7] = dYlm_dphi_l14m_7;
  dY_dphi_[14][8] = dYlm_dphi_l14m_8;
  dY_dphi_[14][9] = dYlm_dphi_l14m_9;
  dY_dphi_[14][10] = dYlm_dphi_l14m_10;
  dY_dphi_[14][11] = dYlm_dphi_l14m_11;
  dY_dphi_[14][12] = dYlm_dphi_l14m_12;
  dY_dphi_[14][13] = dYlm_dphi_l14m_13;
  dY_dphi_[14][14] = dYlm_dphi_l14m_14;
}
/* ->return value: dY_{l}^{m}(x)/dphi */
double complex dYlm_dphi(const unsigned l, int m, const double theta, const double phi)
{
  if (theta > M_PI || theta < 0)
    Error0("theta exceeds from [0,pi] interval.\n");
  if (phi > 2*M_PI || phi < 0)
    Error0("phi exceeds from [0,2*pi] interval.\n");
  if (l >= 15)
    Error0("l exceeds the maximum.\n");
  if ((unsigned)abs(m) > l)
    return 0;
  if (m < 0)
    return dY_dphi_[l][-m](theta,phi);
  return dY_dphi[l][m](theta,phi);
}



/* initialize dYlm_dtheta */
void init_dYlm_dtheta(void)
{
  dY_dtheta[0][0] = dYlm_dtheta_l0m0;
  dY_dtheta[1][0] = dYlm_dtheta_l1m0;
  dY_dtheta[1][1] = dYlm_dtheta_l1m1;
  dY_dtheta[2][0] = dYlm_dtheta_l2m0;
  dY_dtheta[2][1] = dYlm_dtheta_l2m1;
  dY_dtheta[2][2] = dYlm_dtheta_l2m2;
  dY_dtheta[3][0] = dYlm_dtheta_l3m0;
  dY_dtheta[3][1] = dYlm_dtheta_l3m1;
  dY_dtheta[3][2] = dYlm_dtheta_l3m2;
  dY_dtheta[3][3] = dYlm_dtheta_l3m3;
  dY_dtheta[4][0] = dYlm_dtheta_l4m0;
  dY_dtheta[4][1] = dYlm_dtheta_l4m1;
  dY_dtheta[4][2] = dYlm_dtheta_l4m2;
  dY_dtheta[4][3] = dYlm_dtheta_l4m3;
  dY_dtheta[4][4] = dYlm_dtheta_l4m4;
  dY_dtheta[5][0] = dYlm_dtheta_l5m0;
  dY_dtheta[5][1] = dYlm_dtheta_l5m1;
  dY_dtheta[5][2] = dYlm_dtheta_l5m2;
  dY_dtheta[5][3] = dYlm_dtheta_l5m3;
  dY_dtheta[5][4] = dYlm_dtheta_l5m4;
  dY_dtheta[5][5] = dYlm_dtheta_l5m5;
  dY_dtheta[6][0] = dYlm_dtheta_l6m0;
  dY_dtheta[6][1] = dYlm_dtheta_l6m1;
  dY_dtheta[6][2] = dYlm_dtheta_l6m2;
  dY_dtheta[6][3] = dYlm_dtheta_l6m3;
  dY_dtheta[6][4] = dYlm_dtheta_l6m4;
  dY_dtheta[6][5] = dYlm_dtheta_l6m5;
  dY_dtheta[6][6] = dYlm_dtheta_l6m6;
  dY_dtheta[7][0] = dYlm_dtheta_l7m0;
  dY_dtheta[7][1] = dYlm_dtheta_l7m1;
  dY_dtheta[7][2] = dYlm_dtheta_l7m2;
  dY_dtheta[7][3] = dYlm_dtheta_l7m3;
  dY_dtheta[7][4] = dYlm_dtheta_l7m4;
  dY_dtheta[7][5] = dYlm_dtheta_l7m5;
  dY_dtheta[7][6] = dYlm_dtheta_l7m6;
  dY_dtheta[7][7] = dYlm_dtheta_l7m7;
  dY_dtheta[8][0] = dYlm_dtheta_l8m0;
  dY_dtheta[8][1] = dYlm_dtheta_l8m1;
  dY_dtheta[8][2] = dYlm_dtheta_l8m2;
  dY_dtheta[8][3] = dYlm_dtheta_l8m3;
  dY_dtheta[8][4] = dYlm_dtheta_l8m4;
  dY_dtheta[8][5] = dYlm_dtheta_l8m5;
  dY_dtheta[8][6] = dYlm_dtheta_l8m6;
  dY_dtheta[8][7] = dYlm_dtheta_l8m7;
  dY_dtheta[8][8] = dYlm_dtheta_l8m8;
  dY_dtheta[9][0] = dYlm_dtheta_l9m0;
  dY_dtheta[9][1] = dYlm_dtheta_l9m1;
  dY_dtheta[9][2] = dYlm_dtheta_l9m2;
  dY_dtheta[9][3] = dYlm_dtheta_l9m3;
  dY_dtheta[9][4] = dYlm_dtheta_l9m4;
  dY_dtheta[9][5] = dYlm_dtheta_l9m5;
  dY_dtheta[9][6] = dYlm_dtheta_l9m6;
  dY_dtheta[9][7] = dYlm_dtheta_l9m7;
  dY_dtheta[9][8] = dYlm_dtheta_l9m8;
  dY_dtheta[9][9] = dYlm_dtheta_l9m9;
  dY_dtheta[10][0] = dYlm_dtheta_l10m0;
  dY_dtheta[10][1] = dYlm_dtheta_l10m1;
  dY_dtheta[10][2] = dYlm_dtheta_l10m2;
  dY_dtheta[10][3] = dYlm_dtheta_l10m3;
  dY_dtheta[10][4] = dYlm_dtheta_l10m4;
  dY_dtheta[10][5] = dYlm_dtheta_l10m5;
  dY_dtheta[10][6] = dYlm_dtheta_l10m6;
  dY_dtheta[10][7] = dYlm_dtheta_l10m7;
  dY_dtheta[10][8] = dYlm_dtheta_l10m8;
  dY_dtheta[10][9] = dYlm_dtheta_l10m9;
  dY_dtheta[10][10] = dYlm_dtheta_l10m10;
  dY_dtheta[11][0] = dYlm_dtheta_l11m0;
  dY_dtheta[11][1] = dYlm_dtheta_l11m1;
  dY_dtheta[11][2] = dYlm_dtheta_l11m2;
  dY_dtheta[11][3] = dYlm_dtheta_l11m3;
  dY_dtheta[11][4] = dYlm_dtheta_l11m4;
  dY_dtheta[11][5] = dYlm_dtheta_l11m5;
  dY_dtheta[11][6] = dYlm_dtheta_l11m6;
  dY_dtheta[11][7] = dYlm_dtheta_l11m7;
  dY_dtheta[11][8] = dYlm_dtheta_l11m8;
  dY_dtheta[11][9] = dYlm_dtheta_l11m9;
  dY_dtheta[11][10] = dYlm_dtheta_l11m10;
  dY_dtheta[11][11] = dYlm_dtheta_l11m11;
  dY_dtheta[12][0] = dYlm_dtheta_l12m0;
  dY_dtheta[12][1] = dYlm_dtheta_l12m1;
  dY_dtheta[12][2] = dYlm_dtheta_l12m2;
  dY_dtheta[12][3] = dYlm_dtheta_l12m3;
  dY_dtheta[12][4] = dYlm_dtheta_l12m4;
  dY_dtheta[12][5] = dYlm_dtheta_l12m5;
  dY_dtheta[12][6] = dYlm_dtheta_l12m6;
  dY_dtheta[12][7] = dYlm_dtheta_l12m7;
  dY_dtheta[12][8] = dYlm_dtheta_l12m8;
  dY_dtheta[12][9] = dYlm_dtheta_l12m9;
  dY_dtheta[12][10] = dYlm_dtheta_l12m10;
  dY_dtheta[12][11] = dYlm_dtheta_l12m11;
  dY_dtheta[12][12] = dYlm_dtheta_l12m12;
  dY_dtheta[13][0] = dYlm_dtheta_l13m0;
  dY_dtheta[13][1] = dYlm_dtheta_l13m1;
  dY_dtheta[13][2] = dYlm_dtheta_l13m2;
  dY_dtheta[13][3] = dYlm_dtheta_l13m3;
  dY_dtheta[13][4] = dYlm_dtheta_l13m4;
  dY_dtheta[13][5] = dYlm_dtheta_l13m5;
  dY_dtheta[13][6] = dYlm_dtheta_l13m6;
  dY_dtheta[13][7] = dYlm_dtheta_l13m7;
  dY_dtheta[13][8] = dYlm_dtheta_l13m8;
  dY_dtheta[13][9] = dYlm_dtheta_l13m9;
  dY_dtheta[13][10] = dYlm_dtheta_l13m10;
  dY_dtheta[13][11] = dYlm_dtheta_l13m11;
  dY_dtheta[13][12] = dYlm_dtheta_l13m12;
  dY_dtheta[13][13] = dYlm_dtheta_l13m13;
  dY_dtheta[14][0] = dYlm_dtheta_l14m0;
  dY_dtheta[14][1] = dYlm_dtheta_l14m1;
  dY_dtheta[14][2] = dYlm_dtheta_l14m2;
  dY_dtheta[14][3] = dYlm_dtheta_l14m3;
  dY_dtheta[14][4] = dYlm_dtheta_l14m4;
  dY_dtheta[14][5] = dYlm_dtheta_l14m5;
  dY_dtheta[14][6] = dYlm_dtheta_l14m6;
  dY_dtheta[14][7] = dYlm_dtheta_l14m7;
  dY_dtheta[14][8] = dYlm_dtheta_l14m8;
  dY_dtheta[14][9] = dYlm_dtheta_l14m9;
  dY_dtheta[14][10] = dYlm_dtheta_l14m10;
  dY_dtheta[14][11] = dYlm_dtheta_l14m11;
  dY_dtheta[14][12] = dYlm_dtheta_l14m12;
  dY_dtheta[14][13] = dYlm_dtheta_l14m13;
  dY_dtheta[14][14] = dYlm_dtheta_l14m14;
  dY_dtheta_[1][1] = dYlm_dtheta_l1m_1;
  dY_dtheta_[2][1] = dYlm_dtheta_l2m_1;
  dY_dtheta_[2][2] = dYlm_dtheta_l2m_2;
  dY_dtheta_[3][1] = dYlm_dtheta_l3m_1;
  dY_dtheta_[3][2] = dYlm_dtheta_l3m_2;
  dY_dtheta_[3][3] = dYlm_dtheta_l3m_3;
  dY_dtheta_[4][1] = dYlm_dtheta_l4m_1;
  dY_dtheta_[4][2] = dYlm_dtheta_l4m_2;
  dY_dtheta_[4][3] = dYlm_dtheta_l4m_3;
  dY_dtheta_[4][4] = dYlm_dtheta_l4m_4;
  dY_dtheta_[5][1] = dYlm_dtheta_l5m_1;
  dY_dtheta_[5][2] = dYlm_dtheta_l5m_2;
  dY_dtheta_[5][3] = dYlm_dtheta_l5m_3;
  dY_dtheta_[5][4] = dYlm_dtheta_l5m_4;
  dY_dtheta_[5][5] = dYlm_dtheta_l5m_5;
  dY_dtheta_[6][1] = dYlm_dtheta_l6m_1;
  dY_dtheta_[6][2] = dYlm_dtheta_l6m_2;
  dY_dtheta_[6][3] = dYlm_dtheta_l6m_3;
  dY_dtheta_[6][4] = dYlm_dtheta_l6m_4;
  dY_dtheta_[6][5] = dYlm_dtheta_l6m_5;
  dY_dtheta_[6][6] = dYlm_dtheta_l6m_6;
  dY_dtheta_[7][1] = dYlm_dtheta_l7m_1;
  dY_dtheta_[7][2] = dYlm_dtheta_l7m_2;
  dY_dtheta_[7][3] = dYlm_dtheta_l7m_3;
  dY_dtheta_[7][4] = dYlm_dtheta_l7m_4;
  dY_dtheta_[7][5] = dYlm_dtheta_l7m_5;
  dY_dtheta_[7][6] = dYlm_dtheta_l7m_6;
  dY_dtheta_[7][7] = dYlm_dtheta_l7m_7;
  dY_dtheta_[8][1] = dYlm_dtheta_l8m_1;
  dY_dtheta_[8][2] = dYlm_dtheta_l8m_2;
  dY_dtheta_[8][3] = dYlm_dtheta_l8m_3;
  dY_dtheta_[8][4] = dYlm_dtheta_l8m_4;
  dY_dtheta_[8][5] = dYlm_dtheta_l8m_5;
  dY_dtheta_[8][6] = dYlm_dtheta_l8m_6;
  dY_dtheta_[8][7] = dYlm_dtheta_l8m_7;
  dY_dtheta_[8][8] = dYlm_dtheta_l8m_8;
  dY_dtheta_[9][1] = dYlm_dtheta_l9m_1;
  dY_dtheta_[9][2] = dYlm_dtheta_l9m_2;
  dY_dtheta_[9][3] = dYlm_dtheta_l9m_3;
  dY_dtheta_[9][4] = dYlm_dtheta_l9m_4;
  dY_dtheta_[9][5] = dYlm_dtheta_l9m_5;
  dY_dtheta_[9][6] = dYlm_dtheta_l9m_6;
  dY_dtheta_[9][7] = dYlm_dtheta_l9m_7;
  dY_dtheta_[9][8] = dYlm_dtheta_l9m_8;
  dY_dtheta_[9][9] = dYlm_dtheta_l9m_9;
  dY_dtheta_[10][1] = dYlm_dtheta_l10m_1;
  dY_dtheta_[10][2] = dYlm_dtheta_l10m_2;
  dY_dtheta_[10][3] = dYlm_dtheta_l10m_3;
  dY_dtheta_[10][4] = dYlm_dtheta_l10m_4;
  dY_dtheta_[10][5] = dYlm_dtheta_l10m_5;
  dY_dtheta_[10][6] = dYlm_dtheta_l10m_6;
  dY_dtheta_[10][7] = dYlm_dtheta_l10m_7;
  dY_dtheta_[10][8] = dYlm_dtheta_l10m_8;
  dY_dtheta_[10][9] = dYlm_dtheta_l10m_9;
  dY_dtheta_[10][10] = dYlm_dtheta_l10m_10;
  dY_dtheta_[11][1] = dYlm_dtheta_l11m_1;
  dY_dtheta_[11][2] = dYlm_dtheta_l11m_2;
  dY_dtheta_[11][3] = dYlm_dtheta_l11m_3;
  dY_dtheta_[11][4] = dYlm_dtheta_l11m_4;
  dY_dtheta_[11][5] = dYlm_dtheta_l11m_5;
  dY_dtheta_[11][6] = dYlm_dtheta_l11m_6;
  dY_dtheta_[11][7] = dYlm_dtheta_l11m_7;
  dY_dtheta_[11][8] = dYlm_dtheta_l11m_8;
  dY_dtheta_[11][9] = dYlm_dtheta_l11m_9;
  dY_dtheta_[11][10] = dYlm_dtheta_l11m_10;
  dY_dtheta_[11][11] = dYlm_dtheta_l11m_11;
  dY_dtheta_[12][1] = dYlm_dtheta_l12m_1;
  dY_dtheta_[12][2] = dYlm_dtheta_l12m_2;
  dY_dtheta_[12][3] = dYlm_dtheta_l12m_3;
  dY_dtheta_[12][4] = dYlm_dtheta_l12m_4;
  dY_dtheta_[12][5] = dYlm_dtheta_l12m_5;
  dY_dtheta_[12][6] = dYlm_dtheta_l12m_6;
  dY_dtheta_[12][7] = dYlm_dtheta_l12m_7;
  dY_dtheta_[12][8] = dYlm_dtheta_l12m_8;
  dY_dtheta_[12][9] = dYlm_dtheta_l12m_9;
  dY_dtheta_[12][10] = dYlm_dtheta_l12m_10;
  dY_dtheta_[12][11] = dYlm_dtheta_l12m_11;
  dY_dtheta_[12][12] = dYlm_dtheta_l12m_12;
  dY_dtheta_[13][1] = dYlm_dtheta_l13m_1;
  dY_dtheta_[13][2] = dYlm_dtheta_l13m_2;
  dY_dtheta_[13][3] = dYlm_dtheta_l13m_3;
  dY_dtheta_[13][4] = dYlm_dtheta_l13m_4;
  dY_dtheta_[13][5] = dYlm_dtheta_l13m_5;
  dY_dtheta_[13][6] = dYlm_dtheta_l13m_6;
  dY_dtheta_[13][7] = dYlm_dtheta_l13m_7;
  dY_dtheta_[13][8] = dYlm_dtheta_l13m_8;
  dY_dtheta_[13][9] = dYlm_dtheta_l13m_9;
  dY_dtheta_[13][10] = dYlm_dtheta_l13m_10;
  dY_dtheta_[13][11] = dYlm_dtheta_l13m_11;
  dY_dtheta_[13][12] = dYlm_dtheta_l13m_12;
  dY_dtheta_[13][13] = dYlm_dtheta_l13m_13;
  dY_dtheta_[14][1] = dYlm_dtheta_l14m_1;
  dY_dtheta_[14][2] = dYlm_dtheta_l14m_2;
  dY_dtheta_[14][3] = dYlm_dtheta_l14m_3;
  dY_dtheta_[14][4] = dYlm_dtheta_l14m_4;
  dY_dtheta_[14][5] = dYlm_dtheta_l14m_5;
  dY_dtheta_[14][6] = dYlm_dtheta_l14m_6;
  dY_dtheta_[14][7] = dYlm_dtheta_l14m_7;
  dY_dtheta_[14][8] = dYlm_dtheta_l14m_8;
  dY_dtheta_[14][9] = dYlm_dtheta_l14m_9;
  dY_dtheta_[14][10] = dYlm_dtheta_l14m_10;
  dY_dtheta_[14][11] = dYlm_dtheta_l14m_11;
  dY_dtheta_[14][12] = dYlm_dtheta_l14m_12;
  dY_dtheta_[14][13] = dYlm_dtheta_l14m_13;
  dY_dtheta_[14][14] = dYlm_dtheta_l14m_14;
}
/* ->return value: dY_{l}^{m}(x)/dtheta */
double complex dYlm_dtheta(const unsigned l, int m, const double theta, const double phi)
{
  if (theta > M_PI || theta < 0)
    Error0("theta exceeds from [0,pi] interval.\n");
  if (phi > 2*M_PI || phi < 0)
    Error0("phi exceeds from [0,2*pi] interval.\n");
  if (l >= 15)
    Error0("l exceeds the maximum.\n");
  if ((unsigned)abs(m) > l)
    return 0;
  if (m < 0)
    return dY_dtheta_[l][-m](theta,phi);
  return dY_dtheta[l][m](theta,phi);
}

/* Y_{0}^{0} */
static double complex Ylm_l0m0(const double theta, const double phi)
{
  UNUSED(theta);
  UNUSED(phi);
  return (1.0/2.0)/sqrt(M_PI);
}


/* Y_{1}^{0} */
static double complex Ylm_l1m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/2.0)*sqrt(3)*cos(theta)/sqrt(M_PI);
}


/* Y_{1}^{1} */
static double complex Ylm_l1m1(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(6)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{2}^{0} */
static double complex Ylm_l2m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/4.0)*sqrt(5)*(3*pow(cos(theta), 2) - 1)/sqrt(M_PI);
}


/* Y_{2}^{1} */
static double complex Ylm_l2m1(const double theta, const double phi)
{
  return -1.0/8.0*sqrt(30)*cexp(I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* Y_{2}^{2} */
static double complex Ylm_l2m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{3}^{0} */
static double complex Ylm_l3m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/4.0)*sqrt(7)*(5*pow(cos(theta), 2) - 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{3}^{1} */
static double complex Ylm_l3m1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*(1 - 5*pow(cos(theta), 2))*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{3}^{2} */
static double complex Ylm_l3m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(210)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{3}^{3} */
static double complex Ylm_l3m3(const double theta, const double phi)
{
  return -1.0/8.0*sqrt(35)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{4}^{0} */
static double complex Ylm_l4m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (3.0/16.0)*(35*pow(cos(theta), 4) - 30*pow(cos(theta), 2) + 3)/sqrt(M_PI);
}


/* Y_{4}^{1} */
static double complex Ylm_l4m1(const double theta, const double phi)
{
  return -sqrt(5)*((3.0/32.0)*sin(2*theta) + (21.0/64.0)*sin(4*theta))*cexp(I*phi)/sqrt(M_PI);
}


/* Y_{4}^{2} */
static double complex Ylm_l4m2(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(10)*(6 - 7*pow(sin(theta), 2))*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{4}^{3} */
static double complex Ylm_l4m3(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(35)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{4}^{4} */
static double complex Ylm_l4m4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(70)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{5}^{0} */
static double complex Ylm_l5m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/16.0)*sqrt(11)*(63*pow(cos(theta), 4) - 70*pow(cos(theta), 2) + 15)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{1} */
static double complex Ylm_l5m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*(-21*pow(cos(theta), 4) + 14*pow(cos(theta), 2) - 1)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{5}^{2} */
static double complex Ylm_l5m2(const double theta, const double phi)
{
  return (1.0/16.0)*sqrt(2310)*(2 - 3*pow(sin(theta), 2))*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{3} */
static double complex Ylm_l5m3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(385)*(1 - 9*pow(cos(theta), 2))*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{5}^{4} */
static double complex Ylm_l5m4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(770)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{5} */
static double complex Ylm_l5m5(const double theta, const double phi)
{
  return -3.0/32.0*sqrt(77)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{6}^{0} */
static double complex Ylm_l6m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/32.0)*sqrt(13)*(231*pow(cos(theta), 6) - 315*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)/sqrt(M_PI);
}


/* Y_{6}^{1} */
static double complex Ylm_l6m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*(-33*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 5)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{2} */
static double complex Ylm_l6m2(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(1365)*(33*pow(sin(theta), 4) - 48*pow(sin(theta), 2) + 16)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{6}^{3} */
static double complex Ylm_l6m3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*(3 - 11*pow(cos(theta), 2))*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{4} */
static double complex Ylm_l6m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(182)*(11*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{6}^{5} */
static double complex Ylm_l6m5(const double theta, const double phi)
{
  return -3.0/32.0*sqrt(1001)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{6} */
static double complex Ylm_l6m6(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(3003)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{7}^{0} */
static double complex Ylm_l7m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/32.0)*sqrt(15)*(429*pow(cos(theta), 6) - 693*pow(cos(theta), 4) + 315*pow(cos(theta), 2) - 35)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{1} */
static double complex Ylm_l7m1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*(-429*pow(cos(theta), 6) + 495*pow(cos(theta), 4) - 135*pow(cos(theta), 2) + 5)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{7}^{2} */
static double complex Ylm_l7m2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(35)*(143*pow(sin(theta), 4) - 176*pow(sin(theta), 2) + 48)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{3} */
static double complex Ylm_l7m3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(70)*(-143*pow(cos(theta), 4) + 66*pow(cos(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{7}^{4} */
static double complex Ylm_l7m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(13*pow(cos(theta), 2) - 3)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{5} */
static double complex Ylm_l7m5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(770)*(1 - 13*pow(cos(theta), 2))*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{7}^{6} */
static double complex Ylm_l7m6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(5005)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{7} */
static double complex Ylm_l7m7(const double theta, const double phi)
{
  return -3.0/128.0*sqrt(1430)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{8}^{0} */
static double complex Ylm_l8m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/256.0)*sqrt(17)*(6435*pow(cos(theta), 8) - 12012*pow(cos(theta), 6) + 6930*pow(cos(theta), 4) - 1260*pow(cos(theta), 2) + 35)/sqrt(M_PI);
}


/* Y_{8}^{1} */
static double complex Ylm_l8m1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*(-715*pow(cos(theta), 6) + 1001*pow(cos(theta), 4) - 385*pow(cos(theta), 2) + 35)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{2} */
static double complex Ylm_l8m2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(595)*(-143*pow(sin(theta), 4) + 253*pow(sin(theta), 2) + 143*pow(cos(theta), 6) - 111)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{8}^{3} */
static double complex Ylm_l8m3(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(39270)*(-39*pow(cos(theta), 4) + 26*pow(cos(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{4} */
static double complex Ylm_l8m4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2618)*(65*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{8}^{5} */
static double complex Ylm_l8m5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34034)*(1 - 5*pow(cos(theta), 2))*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{6} */
static double complex Ylm_l8m6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(7293)*(15*pow(cos(theta), 2) - 1)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{8}^{7} */
static double complex Ylm_l8m7(const double theta, const double phi)
{
  return -3.0/128.0*sqrt(24310)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{8} */
static double complex Ylm_l8m8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(24310)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{9}^{0} */
static double complex Ylm_l9m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/256.0)*sqrt(19)*(12155*pow(cos(theta), 8) - 25740*pow(cos(theta), 6) + 18018*pow(cos(theta), 4) - 4620*pow(cos(theta), 2) + 315)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{1} */
static double complex Ylm_l9m1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*(-2431*pow(cos(theta), 8) + 4004*pow(cos(theta), 6) - 2002*pow(cos(theta), 4) + 308*pow(cos(theta), 2) - 7)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{9}^{2} */
static double complex Ylm_l9m2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1045)*(-273*pow(sin(theta), 4) + 455*pow(sin(theta), 2) + 221*pow(cos(theta), 6) - 189)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{3} */
static double complex Ylm_l9m3(const double theta, const double phi)
{
  return (1.0/256.0)*sqrt(21945)*(-221*pow(cos(theta), 6) + 195*pow(cos(theta), 4) - 39*pow(cos(theta), 2) + 1)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{9}^{4} */
static double complex Ylm_l9m4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(190190)*(17*pow(cos(theta), 4) - 10*pow(cos(theta), 2) + 1)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{5} */
static double complex Ylm_l9m5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2717)*(-85*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 1)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{9}^{6} */
static double complex Ylm_l9m6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(40755)*(17*pow(cos(theta), 2) - 3)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{7} */
static double complex Ylm_l9m7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(13585)*(1 - 17*pow(cos(theta), 2))*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{9}^{8} */
static double complex Ylm_l9m8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(461890)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{9} */
static double complex Ylm_l9m9(const double theta, const double phi)
{
  return -1.0/512.0*sqrt(230945)*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{10}^{0} */
static double complex Ylm_l10m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/512.0)*sqrt(21)*(46189*pow(cos(theta), 10) - 109395*pow(cos(theta), 8) + 90090*pow(cos(theta), 6) - 30030*pow(cos(theta), 4) + 3465*pow(cos(theta), 2) - 63)/sqrt(M_PI);
}


/* Y_{10}^{1} */
static double complex Ylm_l10m1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*(-4199*pow(cos(theta), 8) + 7956*pow(cos(theta), 6) - 4914*pow(cos(theta), 4) + 1092*pow(cos(theta), 2) - 63)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{2} */
static double complex Ylm_l10m2(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(770)*(2730*pow(sin(theta), 4) - 5096*pow(sin(theta), 2) + 4199*pow(cos(theta), 8) - 6188*pow(cos(theta), 6) + 2373)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{10}^{3} */
static double complex Ylm_l10m3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(5005)*(-323*pow(cos(theta), 6) + 357*pow(cos(theta), 4) - 105*pow(cos(theta), 2) + 7)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{4} */
static double complex Ylm_l10m4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(10010)*(323*pow(cos(theta), 6) - 255*pow(cos(theta), 4) + 45*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{10}^{5} */
static double complex Ylm_l10m5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(1001)*(-323*pow(cos(theta), 4) + 170*pow(cos(theta), 2) - 15)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{6} */
static double complex Ylm_l10m6(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(5005)*(323*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 3)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{10}^{7} */
static double complex Ylm_l10m7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(85085)*(3 - 19*pow(cos(theta), 2))*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{8} */
static double complex Ylm_l10m8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(510510)*(19*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{10}^{9} */
static double complex Ylm_l10m9(const double theta, const double phi)
{
  return -1.0/512.0*sqrt(4849845)*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{10} */
static double complex Ylm_l10m10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(969969)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{11}^{0} */
static double complex Ylm_l11m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/512.0)*sqrt(23)*(88179*pow(cos(theta), 10) - 230945*pow(cos(theta), 8) + 218790*pow(cos(theta), 6) - 90090*pow(cos(theta), 4) + 15015*pow(cos(theta), 2) - 693)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{1} */
static double complex Ylm_l11m1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*(-29393*pow(cos(theta), 10) + 62985*pow(cos(theta), 8) - 46410*pow(cos(theta), 6) + 13650*pow(cos(theta), 4) - 1365*pow(cos(theta), 2) + 21)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{11}^{2} */
static double complex Ylm_l11m2(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(98670)*(2142*pow(sin(theta), 4) - 3864*pow(sin(theta), 2) + 2261*pow(cos(theta), 8) - 3876*pow(cos(theta), 6) + 1743)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{3} */
static double complex Ylm_l11m3(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(345345)*(-969*pow(cos(theta), 8) + 1292*pow(cos(theta), 6) - 510*pow(cos(theta), 4) + 60*pow(cos(theta), 2) - 1)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{11}^{4} */
static double complex Ylm_l11m4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(46046)*(323*pow(cos(theta), 6) - 323*pow(cos(theta), 4) + 85*pow(cos(theta), 2) - 5)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{5} */
static double complex Ylm_l11m5(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(6578)*(-2261*pow(cos(theta), 6) + 1615*pow(cos(theta), 4) - 255*pow(cos(theta), 2) + 5)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{11}^{6} */
static double complex Ylm_l11m6(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(167739)*(399*pow(cos(theta), 4) - 190*pow(cos(theta), 2) + 15)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{7} */
static double complex Ylm_l11m7(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(1677390)*(-133*pow(cos(theta), 4) + 38*pow(cos(theta), 2) - 1)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{11}^{8} */
static double complex Ylm_l11m8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(31870410)*(7*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{9} */
static double complex Ylm_l11m9(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(2124694)*(1 - 21*pow(cos(theta), 2))*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{11}^{10} */
static double complex Ylm_l11m10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(22309287)*cexp(10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{11} */
static double complex Ylm_l11m11(const double theta, const double phi)
{
  return -1.0/2048.0*sqrt(4056234)*cexp(11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{12}^{0} */
static double complex Ylm_l12m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (5.0/2048.0)*(676039*pow(cos(theta), 12) - 1939938*pow(cos(theta), 10) + 2078505*pow(cos(theta), 8) - 1021020*pow(cos(theta), 6) + 225225*pow(cos(theta), 4) - 18018*pow(cos(theta), 2) + 231)/sqrt(M_PI);
}


/* Y_{12}^{1} */
static double complex Ylm_l12m1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*(-52003*pow(cos(theta), 10) + 124355*pow(cos(theta), 8) - 106590*pow(cos(theta), 6) + 39270*pow(cos(theta), 4) - 5775*pow(cos(theta), 2) + 231)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{2} */
static double complex Ylm_l12m2(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(6006)*(-2550*pow(sin(theta), 4) + 4875*pow(sin(theta), 2) + 7429*pow(cos(theta), 10) - 14535*pow(cos(theta), 8) + 9690*pow(cos(theta), 6) - 2328)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{12}^{3} */
static double complex Ylm_l12m3(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(1001)*(-7429*pow(cos(theta), 8) + 11628*pow(cos(theta), 6) - 5814*pow(cos(theta), 4) + 1020*pow(cos(theta), 2) - 45)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{4} */
static double complex Ylm_l12m4(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(1001)*(7429*pow(cos(theta), 8) - 9044*pow(cos(theta), 6) + 3230*pow(cos(theta), 4) - 340*pow(cos(theta), 2) + 5)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{12}^{5} */
static double complex Ylm_l12m5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(34034)*(-437*pow(cos(theta), 6) + 399*pow(cos(theta), 4) - 95*pow(cos(theta), 2) + 5)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{6} */
static double complex Ylm_l12m6(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(2431)*(3059*pow(cos(theta), 6) - 1995*pow(cos(theta), 4) + 285*pow(cos(theta), 2) - 5)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{12}^{7} */
static double complex Ylm_l12m7(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(277134)*(-161*pow(cos(theta), 4) + 70*pow(cos(theta), 2) - 5)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{8} */
static double complex Ylm_l12m8(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(277134)*(161*pow(cos(theta), 4) - 42*pow(cos(theta), 2) + 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{12}^{9} */
static double complex Ylm_l12m9(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(646646)*(3 - 23*pow(cos(theta), 2))*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{10} */
static double complex Ylm_l12m10(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(88179)*(23*pow(cos(theta), 2) - 1)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{12}^{11} */
static double complex Ylm_l12m11(const double theta, const double phi)
{
  return -5.0/2048.0*sqrt(4056234)*cexp(11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{12} */
static double complex Ylm_l12m12(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(676039)*cexp(12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{13}^{0} */
static double complex Ylm_l13m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (3.0/2048.0)*sqrt(3)*(1300075*pow(cos(theta), 12) - 4056234*pow(cos(theta), 10) + 4849845*pow(cos(theta), 8) - 2771340*pow(cos(theta), 6) + 765765*pow(cos(theta), 4) - 90090*pow(cos(theta), 2) + 3003)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{1} */
static double complex Ylm_l13m1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*(-185725*pow(cos(theta), 12) + 490314*pow(cos(theta), 10) - 479655*pow(cos(theta), 8) + 213180*pow(cos(theta), 6) - 42075*pow(cos(theta), 4) + 2970*pow(cos(theta), 2) - 33)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{13}^{2} */
static double complex Ylm_l13m2(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2730)*(-21318*pow(sin(theta), 4) + 39831*pow(sin(theta), 2) + 37145*pow(cos(theta), 10) - 81719*pow(cos(theta), 8) + 63954*pow(cos(theta), 6) - 18612)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{3} */
static double complex Ylm_l13m3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(30030)*(-37145*pow(cos(theta), 10) + 66861*pow(cos(theta), 8) - 40698*pow(cos(theta), 6) + 9690*pow(cos(theta), 4) - 765*pow(cos(theta), 2) + 9)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{13}^{4} */
static double complex Ylm_l13m4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(51051)*(10925*pow(cos(theta), 8) - 15732*pow(cos(theta), 6) + 7182*pow(cos(theta), 4) - 1140*pow(cos(theta), 2) + 45)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{5} */
static double complex Ylm_l13m5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(102102)*(-10925*pow(cos(theta), 8) + 12236*pow(cos(theta), 6) - 3990*pow(cos(theta), 4) + 380*pow(cos(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{13}^{6} */
static double complex Ylm_l13m6(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(969969)*(575*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{7} */
static double complex Ylm_l13m7(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(692835)*(-805*pow(cos(theta), 6) + 483*pow(cos(theta), 4) - 63*pow(cos(theta), 2) + 1)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{13}^{8} */
static double complex Ylm_l13m8(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(9699690)*(115*pow(cos(theta), 4) - 46*pow(cos(theta), 2) + 3)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{9} */
static double complex Ylm_l13m9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(88179)*(-575*pow(cos(theta), 4) + 138*pow(cos(theta), 2) - 3)*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{13}^{10} */
static double complex Ylm_l13m10(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2028117)*(25*pow(cos(theta), 2) - 3)*cexp(10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{11} */
static double complex Ylm_l13m11(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(4056234)*(1 - 25*pow(cos(theta), 2))*cexp(11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{13}^{12} */
static double complex Ylm_l13m12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(2028117)*cexp(12*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{13} */
static double complex Ylm_l13m13(const double theta, const double phi)
{
  return -15.0/8192.0*sqrt(312018)*cexp(13*I*phi)*pow(sin(theta), 13)/sqrt(M_PI);
}


/* Y_{14}^{0} */
static double complex Ylm_l14m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (1.0/4096.0)*sqrt(29)*(5014575*pow(cos(theta), 14) - 16900975*pow(cos(theta), 12) + 22309287*pow(cos(theta), 10) - 14549535*pow(cos(theta), 8) + 4849845*pow(cos(theta), 6) - 765765*pow(cos(theta), 4) + 45045*pow(cos(theta), 2) - 429)/sqrt(M_PI);
}


/* Y_{14}^{1} */
static double complex Ylm_l14m1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*(-334305*pow(cos(theta), 12) + 965770*pow(cos(theta), 10) - 1062347*pow(cos(theta), 8) + 554268*pow(cos(theta), 6) - 138567*pow(cos(theta), 4) + 14586*pow(cos(theta), 2) - 429)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{2} */
static double complex Ylm_l14m2(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(79170)*(-334305*pow(cos(theta), 14) + 1151495*pow(cos(theta), 12) - 1552661*pow(cos(theta), 10) + 1033923*pow(cos(theta), 8) - 351747*pow(cos(theta), 6) + 56661*pow(cos(theta), 4) - 3399*pow(cos(theta), 2) + 33)*cexp(2*I*phi)/sqrt(M_PI);
}


/* Y_{14}^{3} */
static double complex Ylm_l14m3(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(448630)*(-58995*pow(cos(theta), 10) + 120175*pow(cos(theta), 8) - 86526*pow(cos(theta), 6) + 26334*pow(cos(theta), 4) - 3135*pow(cos(theta), 2) + 99)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{4} */
static double complex Ylm_l14m4(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(2467465)*(6555*pow(cos(theta), 10) - 10925*pow(cos(theta), 8) + 6118*pow(cos(theta), 6) - 1330*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{14}^{5} */
static double complex Ylm_l14m5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(18752734)*(-1725*pow(cos(theta), 8) + 2300*pow(cos(theta), 6) - 966*pow(cos(theta), 4) + 140*pow(cos(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{6} */
static double complex Ylm_l14m6(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(93763670)*(3105*pow(cos(theta), 8) - 3220*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 84*pow(cos(theta), 2) + 1)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{14}^{7} */
static double complex Ylm_l14m7(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(20092215)*(-1035*pow(cos(theta), 6) + 805*pow(cos(theta), 4) - 161*pow(cos(theta), 2) + 7)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{8} */
static double complex Ylm_l14m8(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(25571910)*(1035*pow(cos(theta), 6) - 575*pow(cos(theta), 4) + 69*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{14}^{9} */
static double complex Ylm_l14m9(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(98025655)*(-135*pow(cos(theta), 4) + 50*pow(cos(theta), 2) - 3)*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{10} */
static double complex Ylm_l14m10(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(117630786)*(225*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 1)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{14}^{11} */
static double complex Ylm_l14m11(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(1 - 9*pow(cos(theta), 2))*cexp(11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{12} */
static double complex Ylm_l14m12(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(1508087)*(27*pow(cos(theta), 2) - 1)*cexp(12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{14}^{13} */
static double complex Ylm_l14m13(const double theta, const double phi)
{
  return -15.0/8192.0*sqrt(9048522)*cexp(13*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{14} */
static double complex Ylm_l14m14(const double theta, const double phi)
{
  return (15.0/16384.0)*sqrt(1292646)*cexp(14*I*phi)*pow(sin(theta), 14)/sqrt(M_PI);
}


/* Y_{1}^{-1} */
static double complex Ylm_l1m_1(const double theta, const double phi)
{
  return (1.0/4.0)*sqrt(6)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{2}^{-1} */
static double complex Ylm_l2m_1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(-I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* Y_{2}^{-2} */
static double complex Ylm_l2m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{3}^{-1} */
static double complex Ylm_l3m_1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*(5*pow(cos(theta), 2) - 1)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{3}^{-2} */
static double complex Ylm_l3m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(210)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{3}^{-3} */
static double complex Ylm_l3m_3(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(35)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{4}^{-1} */
static double complex Ylm_l4m_1(const double theta, const double phi)
{
  return sqrt(5)*((3.0/32.0)*sin(2*theta) + (21.0/64.0)*sin(4*theta))*cexp(-I*phi)/sqrt(M_PI);
}


/* Y_{4}^{-2} */
static double complex Ylm_l4m_2(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(10)*(6 - 7*pow(sin(theta), 2))*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{4}^{-3} */
static double complex Ylm_l4m_3(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(35)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{4}^{-4} */
static double complex Ylm_l4m_4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(70)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{5}^{-1} */
static double complex Ylm_l5m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*(21*pow(cos(theta), 4) - 14*pow(cos(theta), 2) + 1)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{5}^{-2} */
static double complex Ylm_l5m_2(const double theta, const double phi)
{
  return (1.0/16.0)*sqrt(2310)*(2 - 3*pow(sin(theta), 2))*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{-3} */
static double complex Ylm_l5m_3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(385)*(9*pow(cos(theta), 2) - 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{5}^{-4} */
static double complex Ylm_l5m_4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(770)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{5}^{-5} */
static double complex Ylm_l5m_5(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(77)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{6}^{-1} */
static double complex Ylm_l6m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*(33*pow(cos(theta), 4) - 30*pow(cos(theta), 2) + 5)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{-2} */
static double complex Ylm_l6m_2(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(1365)*(33*pow(sin(theta), 4) - 48*pow(sin(theta), 2) + 16)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{6}^{-3} */
static double complex Ylm_l6m_3(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*(11*pow(cos(theta), 2) - 3)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{-4} */
static double complex Ylm_l6m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(182)*(11*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{6}^{-5} */
static double complex Ylm_l6m_5(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1001)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{6}^{-6} */
static double complex Ylm_l6m_6(const double theta, const double phi)
{
  return (1.0/64.0)*sqrt(3003)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{7}^{-1} */
static double complex Ylm_l7m_1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*(429*pow(cos(theta), 6) - 495*pow(cos(theta), 4) + 135*pow(cos(theta), 2) - 5)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{7}^{-2} */
static double complex Ylm_l7m_2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(35)*(143*pow(sin(theta), 4) - 176*pow(sin(theta), 2) + 48)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{-3} */
static double complex Ylm_l7m_3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(70)*(143*pow(cos(theta), 4) - 66*pow(cos(theta), 2) + 3)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{7}^{-4} */
static double complex Ylm_l7m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(13*pow(cos(theta), 2) - 3)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{-5} */
static double complex Ylm_l7m_5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(770)*(13*pow(cos(theta), 2) - 1)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{7}^{-6} */
static double complex Ylm_l7m_6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(5005)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{7}^{-7} */
static double complex Ylm_l7m_7(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1430)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{8}^{-1} */
static double complex Ylm_l8m_1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*(715*pow(cos(theta), 6) - 1001*pow(cos(theta), 4) + 385*pow(cos(theta), 2) - 35)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-2} */
static double complex Ylm_l8m_2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(595)*(-143*pow(sin(theta), 4) + 253*pow(sin(theta), 2) + 143*pow(cos(theta), 6) - 111)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{8}^{-3} */
static double complex Ylm_l8m_3(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(39270)*(39*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 3)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-4} */
static double complex Ylm_l8m_4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2618)*(65*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{8}^{-5} */
static double complex Ylm_l8m_5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34034)*(5*pow(cos(theta), 2) - 1)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-6} */
static double complex Ylm_l8m_6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(7293)*(15*pow(cos(theta), 2) - 1)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{8}^{-7} */
static double complex Ylm_l8m_7(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(24310)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{8}^{-8} */
static double complex Ylm_l8m_8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(24310)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{9}^{-1} */
static double complex Ylm_l9m_1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*(2431*pow(cos(theta), 8) - 4004*pow(cos(theta), 6) + 2002*pow(cos(theta), 4) - 308*pow(cos(theta), 2) + 7)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{9}^{-2} */
static double complex Ylm_l9m_2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1045)*(-273*pow(sin(theta), 4) + 455*pow(sin(theta), 2) + 221*pow(cos(theta), 6) - 189)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-3} */
static double complex Ylm_l9m_3(const double theta, const double phi)
{
  return (1.0/256.0)*sqrt(21945)*(221*pow(cos(theta), 6) - 195*pow(cos(theta), 4) + 39*pow(cos(theta), 2) - 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{9}^{-4} */
static double complex Ylm_l9m_4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(190190)*(17*pow(cos(theta), 4) - 10*pow(cos(theta), 2) + 1)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-5} */
static double complex Ylm_l9m_5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(2717)*(85*pow(cos(theta), 4) - 30*pow(cos(theta), 2) + 1)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{9}^{-6} */
static double complex Ylm_l9m_6(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(40755)*(17*pow(cos(theta), 2) - 3)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-7} */
static double complex Ylm_l9m_7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(13585)*(17*pow(cos(theta), 2) - 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{9}^{-8} */
static double complex Ylm_l9m_8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(461890)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{9}^{-9} */
static double complex Ylm_l9m_9(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(230945)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{10}^{-1} */
static double complex Ylm_l10m_1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*(4199*pow(cos(theta), 8) - 7956*pow(cos(theta), 6) + 4914*pow(cos(theta), 4) - 1092*pow(cos(theta), 2) + 63)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-2} */
static double complex Ylm_l10m_2(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(770)*(2730*pow(sin(theta), 4) - 5096*pow(sin(theta), 2) + 4199*pow(cos(theta), 8) - 6188*pow(cos(theta), 6) + 2373)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{10}^{-3} */
static double complex Ylm_l10m_3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(5005)*(323*pow(cos(theta), 6) - 357*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 7)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-4} */
static double complex Ylm_l10m_4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(10010)*(323*pow(cos(theta), 6) - 255*pow(cos(theta), 4) + 45*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{10}^{-5} */
static double complex Ylm_l10m_5(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(1001)*(323*pow(cos(theta), 4) - 170*pow(cos(theta), 2) + 15)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-6} */
static double complex Ylm_l10m_6(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(5005)*(323*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 3)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{10}^{-7} */
static double complex Ylm_l10m_7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(85085)*(19*pow(cos(theta), 2) - 3)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-8} */
static double complex Ylm_l10m_8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(510510)*(19*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{10}^{-9} */
static double complex Ylm_l10m_9(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(4849845)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{10}^{-10} */
static double complex Ylm_l10m_10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(969969)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{11}^{-1} */
static double complex Ylm_l11m_1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*(29393*pow(cos(theta), 10) - 62985*pow(cos(theta), 8) + 46410*pow(cos(theta), 6) - 13650*pow(cos(theta), 4) + 1365*pow(cos(theta), 2) - 21)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{11}^{-2} */
static double complex Ylm_l11m_2(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(98670)*(2142*pow(sin(theta), 4) - 3864*pow(sin(theta), 2) + 2261*pow(cos(theta), 8) - 3876*pow(cos(theta), 6) + 1743)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-3} */
static double complex Ylm_l11m_3(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(345345)*(969*pow(cos(theta), 8) - 1292*pow(cos(theta), 6) + 510*pow(cos(theta), 4) - 60*pow(cos(theta), 2) + 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{11}^{-4} */
static double complex Ylm_l11m_4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(46046)*(323*pow(cos(theta), 6) - 323*pow(cos(theta), 4) + 85*pow(cos(theta), 2) - 5)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-5} */
static double complex Ylm_l11m_5(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(6578)*(2261*pow(cos(theta), 6) - 1615*pow(cos(theta), 4) + 255*pow(cos(theta), 2) - 5)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{11}^{-6} */
static double complex Ylm_l11m_6(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(167739)*(399*pow(cos(theta), 4) - 190*pow(cos(theta), 2) + 15)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-7} */
static double complex Ylm_l11m_7(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(1677390)*(133*pow(cos(theta), 4) - 38*pow(cos(theta), 2) + 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{11}^{-8} */
static double complex Ylm_l11m_8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(31870410)*(7*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-9} */
static double complex Ylm_l11m_9(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(2124694)*(21*pow(cos(theta), 2) - 1)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{11}^{-10} */
static double complex Ylm_l11m_10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(22309287)*cexp(-10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{11}^{-11} */
static double complex Ylm_l11m_11(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(4056234)*cexp(-11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{12}^{-1} */
static double complex Ylm_l12m_1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*(52003*pow(cos(theta), 10) - 124355*pow(cos(theta), 8) + 106590*pow(cos(theta), 6) - 39270*pow(cos(theta), 4) + 5775*pow(cos(theta), 2) - 231)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-2} */
static double complex Ylm_l12m_2(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(6006)*(-2550*pow(sin(theta), 4) + 4875*pow(sin(theta), 2) + 7429*pow(cos(theta), 10) - 14535*pow(cos(theta), 8) + 9690*pow(cos(theta), 6) - 2328)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* Y_{12}^{-3} */
static double complex Ylm_l12m_3(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(1001)*(7429*pow(cos(theta), 8) - 11628*pow(cos(theta), 6) + 5814*pow(cos(theta), 4) - 1020*pow(cos(theta), 2) + 45)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-4} */
static double complex Ylm_l12m_4(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(1001)*(7429*pow(cos(theta), 8) - 9044*pow(cos(theta), 6) + 3230*pow(cos(theta), 4) - 340*pow(cos(theta), 2) + 5)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{12}^{-5} */
static double complex Ylm_l12m_5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(34034)*(437*pow(cos(theta), 6) - 399*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 5)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-6} */
static double complex Ylm_l12m_6(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(2431)*(3059*pow(cos(theta), 6) - 1995*pow(cos(theta), 4) + 285*pow(cos(theta), 2) - 5)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{12}^{-7} */
static double complex Ylm_l12m_7(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(277134)*(161*pow(cos(theta), 4) - 70*pow(cos(theta), 2) + 5)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-8} */
static double complex Ylm_l12m_8(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(277134)*(161*pow(cos(theta), 4) - 42*pow(cos(theta), 2) + 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{12}^{-9} */
static double complex Ylm_l12m_9(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(646646)*(23*pow(cos(theta), 2) - 3)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-10} */
static double complex Ylm_l12m_10(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(88179)*(23*pow(cos(theta), 2) - 1)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{12}^{-11} */
static double complex Ylm_l12m_11(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(4056234)*cexp(-11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{12}^{-12} */
static double complex Ylm_l12m_12(const double theta, const double phi)
{
  return (5.0/4096.0)*sqrt(676039)*cexp(-12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{13}^{-1} */
static double complex Ylm_l13m_1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*(185725*pow(cos(theta), 12) - 490314*pow(cos(theta), 10) + 479655*pow(cos(theta), 8) - 213180*pow(cos(theta), 6) + 42075*pow(cos(theta), 4) - 2970*pow(cos(theta), 2) + 33)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* Y_{13}^{-2} */
static double complex Ylm_l13m_2(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2730)*(-21318*pow(sin(theta), 4) + 39831*pow(sin(theta), 2) + 37145*pow(cos(theta), 10) - 81719*pow(cos(theta), 8) + 63954*pow(cos(theta), 6) - 18612)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-3} */
static double complex Ylm_l13m_3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(30030)*(37145*pow(cos(theta), 10) - 66861*pow(cos(theta), 8) + 40698*pow(cos(theta), 6) - 9690*pow(cos(theta), 4) + 765*pow(cos(theta), 2) - 9)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* Y_{13}^{-4} */
static double complex Ylm_l13m_4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(51051)*(10925*pow(cos(theta), 8) - 15732*pow(cos(theta), 6) + 7182*pow(cos(theta), 4) - 1140*pow(cos(theta), 2) + 45)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-5} */
static double complex Ylm_l13m_5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(102102)*(10925*pow(cos(theta), 8) - 12236*pow(cos(theta), 6) + 3990*pow(cos(theta), 4) - 380*pow(cos(theta), 2) + 5)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* Y_{13}^{-6} */
static double complex Ylm_l13m_6(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(969969)*(575*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-7} */
static double complex Ylm_l13m_7(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(692835)*(805*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 63*pow(cos(theta), 2) - 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* Y_{13}^{-8} */
static double complex Ylm_l13m_8(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(9699690)*(115*pow(cos(theta), 4) - 46*pow(cos(theta), 2) + 3)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-9} */
static double complex Ylm_l13m_9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(88179)*(575*pow(cos(theta), 4) - 138*pow(cos(theta), 2) + 3)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* Y_{13}^{-10} */
static double complex Ylm_l13m_10(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2028117)*(25*pow(cos(theta), 2) - 3)*cexp(-10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-11} */
static double complex Ylm_l13m_11(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(4056234)*(25*pow(cos(theta), 2) - 1)*cexp(-11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* Y_{13}^{-12} */
static double complex Ylm_l13m_12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(2028117)*cexp(-12*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* Y_{13}^{-13} */
static double complex Ylm_l13m_13(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(312018)*cexp(-13*I*phi)*pow(sin(theta), 13)/sqrt(M_PI);
}


/* Y_{14}^{-1} */
static double complex Ylm_l14m_1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*(334305*pow(cos(theta), 12) - 965770*pow(cos(theta), 10) + 1062347*pow(cos(theta), 8) - 554268*pow(cos(theta), 6) + 138567*pow(cos(theta), 4) - 14586*pow(cos(theta), 2) + 429)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-2} */
static double complex Ylm_l14m_2(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(79170)*(-334305*pow(cos(theta), 14) + 1151495*pow(cos(theta), 12) - 1552661*pow(cos(theta), 10) + 1033923*pow(cos(theta), 8) - 351747*pow(cos(theta), 6) + 56661*pow(cos(theta), 4) - 3399*pow(cos(theta), 2) + 33)*cexp(-2*I*phi)/sqrt(M_PI);
}


/* Y_{14}^{-3} */
static double complex Ylm_l14m_3(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(448630)*(58995*pow(cos(theta), 10) - 120175*pow(cos(theta), 8) + 86526*pow(cos(theta), 6) - 26334*pow(cos(theta), 4) + 3135*pow(cos(theta), 2) - 99)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-4} */
static double complex Ylm_l14m_4(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(2467465)*(6555*pow(cos(theta), 10) - 10925*pow(cos(theta), 8) + 6118*pow(cos(theta), 6) - 1330*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* Y_{14}^{-5} */
static double complex Ylm_l14m_5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(18752734)*(1725*pow(cos(theta), 8) - 2300*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 140*pow(cos(theta), 2) + 5)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-6} */
static double complex Ylm_l14m_6(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(93763670)*(3105*pow(cos(theta), 8) - 3220*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 84*pow(cos(theta), 2) + 1)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* Y_{14}^{-7} */
static double complex Ylm_l14m_7(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(20092215)*(1035*pow(cos(theta), 6) - 805*pow(cos(theta), 4) + 161*pow(cos(theta), 2) - 7)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-8} */
static double complex Ylm_l14m_8(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(25571910)*(1035*pow(cos(theta), 6) - 575*pow(cos(theta), 4) + 69*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* Y_{14}^{-9} */
static double complex Ylm_l14m_9(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(98025655)*(135*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 3)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-10} */
static double complex Ylm_l14m_10(const double theta, const double phi)
{
  return (1.0/16384.0)*sqrt(117630786)*(225*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 1)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* Y_{14}^{-11} */
static double complex Ylm_l14m_11(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(9*pow(cos(theta), 2) - 1)*cexp(-11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-12} */
static double complex Ylm_l14m_12(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(1508087)*(27*pow(cos(theta), 2) - 1)*cexp(-12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* Y_{14}^{-13} */
static double complex Ylm_l14m_13(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(9048522)*cexp(-13*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* Y_{14}^{-14} */
static double complex Ylm_l14m_14(const double theta, const double phi)
{
  return (15.0/16384.0)*sqrt(1292646)*cexp(-14*I*phi)*pow(sin(theta), 14)/sqrt(M_PI);
}


/* dY_dphi_{0}^{0} */
static double complex dYlm_dphi_l0m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{1}^{0} */
static double complex dYlm_dphi_l1m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{1}^{1} */
static double complex dYlm_dphi_l1m1(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(6)*I*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{2}^{0} */
static double complex dYlm_dphi_l2m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{2}^{1} */
static double complex dYlm_dphi_l2m1(const double theta, const double phi)
{
  return -1.0/8.0*sqrt(30)*I*cexp(I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* dY_dphi_{2}^{2} */
static double complex dYlm_dphi_l2m2(const double theta, const double phi)
{
  return (1.0/4.0)*sqrt(30)*I*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{3}^{0} */
static double complex dYlm_dphi_l3m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{3}^{1} */
static double complex dYlm_dphi_l3m1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*I*(1 - 5*pow(cos(theta), 2))*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{3}^{2} */
static double complex dYlm_dphi_l3m2(const double theta, const double phi)
{
  return (1.0/4.0)*sqrt(210)*I*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{3}^{3} */
static double complex dYlm_dphi_l3m3(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(35)*I*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{4}^{0} */
static double complex dYlm_dphi_l4m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{4}^{1} */
static double complex dYlm_dphi_l4m1(const double theta, const double phi)
{
  return -sqrt(5)*I*((3.0/32.0)*sin(2*theta) + (21.0/64.0)*sin(4*theta))*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dphi_{4}^{2} */
static double complex dYlm_dphi_l4m2(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(10)*I*(6 - 7*pow(sin(theta), 2))*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{4}^{3} */
static double complex dYlm_dphi_l4m3(const double theta, const double phi)
{
  return -9.0/8.0*sqrt(35)*I*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{4}^{4} */
static double complex dYlm_dphi_l4m4(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(70)*I*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{5}^{0} */
static double complex dYlm_dphi_l5m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{5}^{1} */
static double complex dYlm_dphi_l5m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*I*(-21*pow(cos(theta), 4) + 14*pow(cos(theta), 2) - 1)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{5}^{2} */
static double complex dYlm_dphi_l5m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(2310)*I*(2 - 3*pow(sin(theta), 2))*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{5}^{3} */
static double complex dYlm_dphi_l5m3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(385)*I*(1 - 9*pow(cos(theta), 2))*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{5}^{4} */
static double complex dYlm_dphi_l5m4(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(770)*I*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{5}^{5} */
static double complex dYlm_dphi_l5m5(const double theta, const double phi)
{
  return -15.0/32.0*sqrt(77)*I*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{6}^{0} */
static double complex dYlm_dphi_l6m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{6}^{1} */
static double complex dYlm_dphi_l6m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*I*(-33*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 5)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{6}^{2} */
static double complex dYlm_dphi_l6m2(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*I*(33*pow(sin(theta), 4) - 48*pow(sin(theta), 2) + 16)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{6}^{3} */
static double complex dYlm_dphi_l6m3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1365)*I*(3 - 11*pow(cos(theta), 2))*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{6}^{4} */
static double complex dYlm_dphi_l6m4(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(182)*I*(11*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{6}^{5} */
static double complex dYlm_dphi_l6m5(const double theta, const double phi)
{
  return -15.0/32.0*sqrt(1001)*I*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{6}^{6} */
static double complex dYlm_dphi_l6m6(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(3003)*I*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{7}^{0} */
static double complex dYlm_dphi_l7m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{7}^{1} */
static double complex dYlm_dphi_l7m1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*I*(-429*pow(cos(theta), 6) + 495*pow(cos(theta), 4) - 135*pow(cos(theta), 2) + 5)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{2} */
static double complex dYlm_dphi_l7m2(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(35)*I*(143*pow(sin(theta), 4) - 176*pow(sin(theta), 2) + 48)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{3} */
static double complex dYlm_dphi_l7m3(const double theta, const double phi)
{
  return (9.0/128.0)*sqrt(70)*I*(-143*pow(cos(theta), 4) + 66*pow(cos(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{7}^{4} */
static double complex dYlm_dphi_l7m4(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(770)*I*(13*pow(cos(theta), 2) - 3)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{5} */
static double complex dYlm_dphi_l7m5(const double theta, const double phi)
{
  return (15.0/128.0)*sqrt(770)*I*(1 - 13*pow(cos(theta), 2))*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{7}^{6} */
static double complex dYlm_dphi_l7m6(const double theta, const double phi)
{
  return (9.0/32.0)*sqrt(5005)*I*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{7} */
static double complex dYlm_dphi_l7m7(const double theta, const double phi)
{
  return -21.0/128.0*sqrt(1430)*I*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{8}^{0} */
static double complex dYlm_dphi_l8m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{8}^{1} */
static double complex dYlm_dphi_l8m1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*I*(-715*pow(cos(theta), 6) + 1001*pow(cos(theta), 4) - 385*pow(cos(theta), 2) + 35)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{2} */
static double complex dYlm_dphi_l8m2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(595)*I*(-143*pow(sin(theta), 4) + 253*pow(sin(theta), 2) + 143*pow(cos(theta), 6) - 111)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{8}^{3} */
static double complex dYlm_dphi_l8m3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(39270)*I*(-39*pow(cos(theta), 4) + 26*pow(cos(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{4} */
static double complex dYlm_dphi_l8m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(2618)*I*(65*pow(cos(theta), 4) - 26*pow(cos(theta), 2) + 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{8}^{5} */
static double complex dYlm_dphi_l8m5(const double theta, const double phi)
{
  return (15.0/128.0)*sqrt(34034)*I*(1 - 5*pow(cos(theta), 2))*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{6} */
static double complex dYlm_dphi_l8m6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(7293)*I*(15*pow(cos(theta), 2) - 1)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{8}^{7} */
static double complex dYlm_dphi_l8m7(const double theta, const double phi)
{
  return -21.0/128.0*sqrt(24310)*I*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{8} */
static double complex dYlm_dphi_l8m8(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(24310)*I*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{9}^{0} */
static double complex dYlm_dphi_l9m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{9}^{1} */
static double complex dYlm_dphi_l9m1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*I*(-2431*pow(cos(theta), 8) + 4004*pow(cos(theta), 6) - 2002*pow(cos(theta), 4) + 308*pow(cos(theta), 2) - 7)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{2} */
static double complex dYlm_dphi_l9m2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(1045)*I*(-273*pow(sin(theta), 4) + 455*pow(sin(theta), 2) + 221*pow(cos(theta), 6) - 189)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{3} */
static double complex dYlm_dphi_l9m3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(21945)*I*(-221*pow(cos(theta), 6) + 195*pow(cos(theta), 4) - 39*pow(cos(theta), 2) + 1)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{9}^{4} */
static double complex dYlm_dphi_l9m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(190190)*I*(17*pow(cos(theta), 4) - 10*pow(cos(theta), 2) + 1)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{5} */
static double complex dYlm_dphi_l9m5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(2717)*I*(-85*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 1)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{9}^{6} */
static double complex dYlm_dphi_l9m6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(40755)*I*(17*pow(cos(theta), 2) - 3)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{7} */
static double complex dYlm_dphi_l9m7(const double theta, const double phi)
{
  return (21.0/512.0)*sqrt(13585)*I*(1 - 17*pow(cos(theta), 2))*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{9}^{8} */
static double complex dYlm_dphi_l9m8(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(461890)*I*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{9} */
static double complex dYlm_dphi_l9m9(const double theta, const double phi)
{
  return -9.0/512.0*sqrt(230945)*I*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dphi_{10}^{0} */
static double complex dYlm_dphi_l10m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{10}^{1} */
static double complex dYlm_dphi_l10m1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*I*(-4199*pow(cos(theta), 8) + 7956*pow(cos(theta), 6) - 4914*pow(cos(theta), 4) + 1092*pow(cos(theta), 2) - 63)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{2} */
static double complex dYlm_dphi_l10m2(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(770)*I*(2730*pow(sin(theta), 4) - 5096*pow(sin(theta), 2) + 4199*pow(cos(theta), 8) - 6188*pow(cos(theta), 6) + 2373)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{10}^{3} */
static double complex dYlm_dphi_l10m3(const double theta, const double phi)
{
  return (9.0/256.0)*sqrt(5005)*I*(-323*pow(cos(theta), 6) + 357*pow(cos(theta), 4) - 105*pow(cos(theta), 2) + 7)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{4} */
static double complex dYlm_dphi_l10m4(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(10010)*I*(323*pow(cos(theta), 6) - 255*pow(cos(theta), 4) + 45*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{10}^{5} */
static double complex dYlm_dphi_l10m5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(1001)*I*(-323*pow(cos(theta), 4) + 170*pow(cos(theta), 2) - 15)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{6} */
static double complex dYlm_dphi_l10m6(const double theta, const double phi)
{
  return (9.0/512.0)*sqrt(5005)*I*(323*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 3)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{10}^{7} */
static double complex dYlm_dphi_l10m7(const double theta, const double phi)
{
  return (21.0/512.0)*sqrt(85085)*I*(3 - 19*pow(cos(theta), 2))*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{8} */
static double complex dYlm_dphi_l10m8(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(510510)*I*(19*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{10}^{9} */
static double complex dYlm_dphi_l10m9(const double theta, const double phi)
{
  return -9.0/512.0*sqrt(4849845)*I*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{10} */
static double complex dYlm_dphi_l10m10(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(969969)*I*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dphi_{11}^{0} */
static double complex dYlm_dphi_l11m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{11}^{1} */
static double complex dYlm_dphi_l11m1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*I*(-29393*pow(cos(theta), 10) + 62985*pow(cos(theta), 8) - 46410*pow(cos(theta), 6) + 13650*pow(cos(theta), 4) - 1365*pow(cos(theta), 2) + 21)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{2} */
static double complex dYlm_dphi_l11m2(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(98670)*I*(2142*pow(sin(theta), 4) - 3864*pow(sin(theta), 2) + 2261*pow(cos(theta), 8) - 3876*pow(cos(theta), 6) + 1743)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{3} */
static double complex dYlm_dphi_l11m3(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(345345)*I*(-969*pow(cos(theta), 8) + 1292*pow(cos(theta), 6) - 510*pow(cos(theta), 4) + 60*pow(cos(theta), 2) - 1)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{11}^{4} */
static double complex dYlm_dphi_l11m4(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(46046)*I*(323*pow(cos(theta), 6) - 323*pow(cos(theta), 4) + 85*pow(cos(theta), 2) - 5)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{5} */
static double complex dYlm_dphi_l11m5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(6578)*I*(-2261*pow(cos(theta), 6) + 1615*pow(cos(theta), 4) - 255*pow(cos(theta), 2) + 5)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{11}^{6} */
static double complex dYlm_dphi_l11m6(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(167739)*I*(399*pow(cos(theta), 4) - 190*pow(cos(theta), 2) + 15)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{7} */
static double complex dYlm_dphi_l11m7(const double theta, const double phi)
{
  return (7.0/2048.0)*sqrt(1677390)*I*(-133*pow(cos(theta), 4) + 38*pow(cos(theta), 2) - 1)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{11}^{8} */
static double complex dYlm_dphi_l11m8(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(31870410)*I*(7*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{9} */
static double complex dYlm_dphi_l11m9(const double theta, const double phi)
{
  return (9.0/2048.0)*sqrt(2124694)*I*(1 - 21*pow(cos(theta), 2))*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dphi_{11}^{10} */
static double complex dYlm_dphi_l11m10(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(22309287)*I*cexp(10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{11} */
static double complex dYlm_dphi_l11m11(const double theta, const double phi)
{
  return -11.0/2048.0*sqrt(4056234)*I*cexp(11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* dY_dphi_{12}^{0} */
static double complex dYlm_dphi_l12m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{12}^{1} */
static double complex dYlm_dphi_l12m1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*I*(-52003*pow(cos(theta), 10) + 124355*pow(cos(theta), 8) - 106590*pow(cos(theta), 6) + 39270*pow(cos(theta), 4) - 5775*pow(cos(theta), 2) + 231)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{2} */
static double complex dYlm_dphi_l12m2(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(6006)*I*(-2550*pow(sin(theta), 4) + 4875*pow(sin(theta), 2) + 7429*pow(cos(theta), 10) - 14535*pow(cos(theta), 8) + 9690*pow(cos(theta), 6) - 2328)*cexp(2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{12}^{3} */
static double complex dYlm_dphi_l12m3(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*I*(-7429*pow(cos(theta), 8) + 11628*pow(cos(theta), 6) - 5814*pow(cos(theta), 4) + 1020*pow(cos(theta), 2) - 45)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{4} */
static double complex dYlm_dphi_l12m4(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*I*(7429*pow(cos(theta), 8) - 9044*pow(cos(theta), 6) + 3230*pow(cos(theta), 4) - 340*pow(cos(theta), 2) + 5)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{12}^{5} */
static double complex dYlm_dphi_l12m5(const double theta, const double phi)
{
  return (75.0/2048.0)*sqrt(34034)*I*(-437*pow(cos(theta), 6) + 399*pow(cos(theta), 4) - 95*pow(cos(theta), 2) + 5)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{6} */
static double complex dYlm_dphi_l12m6(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(2431)*I*(3059*pow(cos(theta), 6) - 1995*pow(cos(theta), 4) + 285*pow(cos(theta), 2) - 5)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{12}^{7} */
static double complex dYlm_dphi_l12m7(const double theta, const double phi)
{
  return (35.0/2048.0)*sqrt(277134)*I*(-161*pow(cos(theta), 4) + 70*pow(cos(theta), 2) - 5)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{8} */
static double complex dYlm_dphi_l12m8(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(277134)*I*(161*pow(cos(theta), 4) - 42*pow(cos(theta), 2) + 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{12}^{9} */
static double complex dYlm_dphi_l12m9(const double theta, const double phi)
{
  return (45.0/2048.0)*sqrt(646646)*I*(3 - 23*pow(cos(theta), 2))*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{10} */
static double complex dYlm_dphi_l12m10(const double theta, const double phi)
{
  return (25.0/1024.0)*sqrt(88179)*I*(23*pow(cos(theta), 2) - 1)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dphi_{12}^{11} */
static double complex dYlm_dphi_l12m11(const double theta, const double phi)
{
  return -55.0/2048.0*sqrt(4056234)*I*cexp(11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{12} */
static double complex dYlm_dphi_l12m12(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(676039)*I*cexp(12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* dY_dphi_{13}^{0} */
static double complex dYlm_dphi_l13m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{13}^{1} */
static double complex dYlm_dphi_l13m1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*I*(-185725*pow(cos(theta), 12) + 490314*pow(cos(theta), 10) - 479655*pow(cos(theta), 8) + 213180*pow(cos(theta), 6) - 42075*pow(cos(theta), 4) + 2970*pow(cos(theta), 2) - 33)*cexp(I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{2} */
static double complex dYlm_dphi_l13m2(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(2730)*I*(-21318*pow(sin(theta), 4) + 39831*pow(sin(theta), 2) + 37145*pow(cos(theta), 10) - 81719*pow(cos(theta), 8) + 63954*pow(cos(theta), 6) - 18612)*cexp(2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{3} */
static double complex dYlm_dphi_l13m3(const double theta, const double phi)
{
  return (9.0/8192.0)*sqrt(30030)*I*(-37145*pow(cos(theta), 10) + 66861*pow(cos(theta), 8) - 40698*pow(cos(theta), 6) + 9690*pow(cos(theta), 4) - 765*pow(cos(theta), 2) + 9)*cexp(3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{13}^{4} */
static double complex dYlm_dphi_l13m4(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(51051)*I*(10925*pow(cos(theta), 8) - 15732*pow(cos(theta), 6) + 7182*pow(cos(theta), 4) - 1140*pow(cos(theta), 2) + 45)*cexp(4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{5} */
static double complex dYlm_dphi_l13m5(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(102102)*I*(-10925*pow(cos(theta), 8) + 12236*pow(cos(theta), 6) - 3990*pow(cos(theta), 4) + 380*pow(cos(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{13}^{6} */
static double complex dYlm_dphi_l13m6(const double theta, const double phi)
{
  return (9.0/1024.0)*sqrt(969969)*I*(575*pow(cos(theta), 6) - 483*pow(cos(theta), 4) + 105*pow(cos(theta), 2) - 5)*cexp(6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{7} */
static double complex dYlm_dphi_l13m7(const double theta, const double phi)
{
  return (21.0/4096.0)*sqrt(692835)*I*(-805*pow(cos(theta), 6) + 483*pow(cos(theta), 4) - 63*pow(cos(theta), 2) + 1)*cexp(7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{13}^{8} */
static double complex dYlm_dphi_l13m8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(9699690)*I*(115*pow(cos(theta), 4) - 46*pow(cos(theta), 2) + 3)*cexp(8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{9} */
static double complex dYlm_dphi_l13m9(const double theta, const double phi)
{
  return (27.0/4096.0)*sqrt(88179)*I*(-575*pow(cos(theta), 4) + 138*pow(cos(theta), 2) - 3)*cexp(9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dphi_{13}^{10} */
static double complex dYlm_dphi_l13m10(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(2028117)*I*(25*pow(cos(theta), 2) - 3)*cexp(10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{11} */
static double complex dYlm_dphi_l13m11(const double theta, const double phi)
{
  return (33.0/8192.0)*sqrt(4056234)*I*(1 - 25*pow(cos(theta), 2))*cexp(11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* dY_dphi_{13}^{12} */
static double complex dYlm_dphi_l13m12(const double theta, const double phi)
{
  return (45.0/1024.0)*sqrt(2028117)*I*cexp(12*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{13} */
static double complex dYlm_dphi_l13m13(const double theta, const double phi)
{
  return -195.0/8192.0*sqrt(312018)*I*cexp(13*I*phi)*pow(sin(theta), 13)/sqrt(M_PI);
}


/* dY_dphi_{14}^{0} */
static double complex dYlm_dphi_l14m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dphi_{14}^{1} */
static double complex dYlm_dphi_l14m1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*I*(-334305*pow(cos(theta), 12) + 965770*pow(cos(theta), 10) - 1062347*pow(cos(theta), 8) + 554268*pow(cos(theta), 6) - 138567*pow(cos(theta), 4) + 14586*pow(cos(theta), 2) - 429)*cexp(I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{2} */
static double complex dYlm_dphi_l14m2(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(79170)*I*(-334305*pow(cos(theta), 14) + 1151495*pow(cos(theta), 12) - 1552661*pow(cos(theta), 10) + 1033923*pow(cos(theta), 8) - 351747*pow(cos(theta), 6) + 56661*pow(cos(theta), 4) - 3399*pow(cos(theta), 2) + 33)*cexp(2*I*phi)/sqrt(M_PI);
}


/* dY_dphi_{14}^{3} */
static double complex dYlm_dphi_l14m3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(448630)*I*(-58995*pow(cos(theta), 10) + 120175*pow(cos(theta), 8) - 86526*pow(cos(theta), 6) + 26334*pow(cos(theta), 4) - 3135*pow(cos(theta), 2) + 99)*cexp(3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{4} */
static double complex dYlm_dphi_l14m4(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2467465)*I*(6555*pow(cos(theta), 10) - 10925*pow(cos(theta), 8) + 6118*pow(cos(theta), 6) - 1330*pow(cos(theta), 4) + 95*pow(cos(theta), 2) - 1)*cexp(4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{14}^{5} */
static double complex dYlm_dphi_l14m5(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(18752734)*I*(-1725*pow(cos(theta), 8) + 2300*pow(cos(theta), 6) - 966*pow(cos(theta), 4) + 140*pow(cos(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{6} */
static double complex dYlm_dphi_l14m6(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(93763670)*I*(3105*pow(cos(theta), 8) - 3220*pow(cos(theta), 6) + 966*pow(cos(theta), 4) - 84*pow(cos(theta), 2) + 1)*cexp(6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{14}^{7} */
static double complex dYlm_dphi_l14m7(const double theta, const double phi)
{
  return (7.0/4096.0)*sqrt(20092215)*I*(-1035*pow(cos(theta), 6) + 805*pow(cos(theta), 4) - 161*pow(cos(theta), 2) + 7)*cexp(7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{8} */
static double complex dYlm_dphi_l14m8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(25571910)*I*(1035*pow(cos(theta), 6) - 575*pow(cos(theta), 4) + 69*pow(cos(theta), 2) - 1)*cexp(8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{14}^{9} */
static double complex dYlm_dphi_l14m9(const double theta, const double phi)
{
  return (9.0/4096.0)*sqrt(98025655)*I*(-135*pow(cos(theta), 4) + 50*pow(cos(theta), 2) - 3)*cexp(9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{10} */
static double complex dYlm_dphi_l14m10(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*I*(225*pow(cos(theta), 4) - 50*pow(cos(theta), 2) + 1)*cexp(10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dphi_{14}^{11} */
static double complex dYlm_dphi_l14m11(const double theta, const double phi)
{
  return (55.0/8192.0)*sqrt(117630786)*I*(1 - 9*pow(cos(theta), 2))*cexp(11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{12} */
static double complex dYlm_dphi_l14m12(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(1508087)*I*(27*pow(cos(theta), 2) - 1)*cexp(12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* dY_dphi_{14}^{13} */
static double complex dYlm_dphi_l14m13(const double theta, const double phi)
{
  return -195.0/8192.0*sqrt(9048522)*I*cexp(13*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{14} */
static double complex dYlm_dphi_l14m14(const double theta, const double phi)
{
  return (105.0/8192.0)*sqrt(1292646)*I*cexp(14*I*phi)*pow(sin(theta), 14)/sqrt(M_PI);
}


/* dY_dphi_{1}^{-1} */
static double complex dYlm_dphi_l1m_1(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(6)*I*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{2}^{-1} */
static double complex dYlm_dphi_l2m_1(const double theta, const double phi)
{
  return -1.0/8.0*sqrt(30)*I*cexp(-I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* dY_dphi_{2}^{-2} */
static double complex dYlm_dphi_l2m_2(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(30)*I*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{3}^{-1} */
static double complex dYlm_dphi_l3m_1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*I*(1 - 5*pow(cos(theta), 2))*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{3}^{-2} */
static double complex dYlm_dphi_l3m_2(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(210)*I*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{3}^{-3} */
static double complex dYlm_dphi_l3m_3(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(35)*I*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{4}^{-1} */
static double complex dYlm_dphi_l4m_1(const double theta, const double phi)
{
  return -sqrt(5)*I*((3.0/32.0)*sin(2*theta) + (21.0/64.0)*sin(4*theta))*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dphi_{4}^{-2} */
static double complex dYlm_dphi_l4m_2(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(10)*I*(7*pow(sin(theta), 2) - 6)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{4}^{-3} */
static double complex dYlm_dphi_l4m_3(const double theta, const double phi)
{
  return -9.0/8.0*sqrt(35)*I*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{4}^{-4} */
static double complex dYlm_dphi_l4m_4(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(70)*I*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{5}^{-1} */
static double complex dYlm_dphi_l5m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*I*(-21*pow(cos(theta), 4) + 14*pow(cos(theta), 2) - 1)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{5}^{-2} */
static double complex dYlm_dphi_l5m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(2310)*I*(3*pow(sin(theta), 2) - 2)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{5}^{-3} */
static double complex dYlm_dphi_l5m_3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(385)*I*(1 - 9*pow(cos(theta), 2))*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{5}^{-4} */
static double complex dYlm_dphi_l5m_4(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(770)*I*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{5}^{-5} */
static double complex dYlm_dphi_l5m_5(const double theta, const double phi)
{
  return -15.0/32.0*sqrt(77)*I*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{6}^{-1} */
static double complex dYlm_dphi_l6m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*I*(-33*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 5)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{6}^{-2} */
static double complex dYlm_dphi_l6m_2(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*I*(-33*pow(sin(theta), 4) + 48*pow(sin(theta), 2) - 16)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{6}^{-3} */
static double complex dYlm_dphi_l6m_3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1365)*I*(3 - 11*pow(cos(theta), 2))*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{6}^{-4} */
static double complex dYlm_dphi_l6m_4(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(182)*I*(1 - 11*pow(cos(theta), 2))*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{6}^{-5} */
static double complex dYlm_dphi_l6m_5(const double theta, const double phi)
{
  return -15.0/32.0*sqrt(1001)*I*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{6}^{-6} */
static double complex dYlm_dphi_l6m_6(const double theta, const double phi)
{
  return -3.0/32.0*sqrt(3003)*I*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-1} */
static double complex dYlm_dphi_l7m_1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*I*(-429*pow(cos(theta), 6) + 495*pow(cos(theta), 4) - 135*pow(cos(theta), 2) + 5)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-2} */
static double complex dYlm_dphi_l7m_2(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(35)*I*(-143*pow(sin(theta), 4) + 176*pow(sin(theta), 2) - 48)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-3} */
static double complex dYlm_dphi_l7m_3(const double theta, const double phi)
{
  return (9.0/128.0)*sqrt(70)*I*(-143*pow(cos(theta), 4) + 66*pow(cos(theta), 2) - 3)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-4} */
static double complex dYlm_dphi_l7m_4(const double theta, const double phi)
{
  return (3.0/16.0)*sqrt(770)*I*(3 - 13*pow(cos(theta), 2))*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-5} */
static double complex dYlm_dphi_l7m_5(const double theta, const double phi)
{
  return (15.0/128.0)*sqrt(770)*I*(1 - 13*pow(cos(theta), 2))*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-6} */
static double complex dYlm_dphi_l7m_6(const double theta, const double phi)
{
  return -9.0/32.0*sqrt(5005)*I*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{7}^{-7} */
static double complex dYlm_dphi_l7m_7(const double theta, const double phi)
{
  return -21.0/128.0*sqrt(1430)*I*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-1} */
static double complex dYlm_dphi_l8m_1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*I*(-715*pow(cos(theta), 6) + 1001*pow(cos(theta), 4) - 385*pow(cos(theta), 2) + 35)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-2} */
static double complex dYlm_dphi_l8m_2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(595)*I*(143*pow(sin(theta), 4) - 253*pow(sin(theta), 2) - 143*pow(cos(theta), 6) + 111)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-3} */
static double complex dYlm_dphi_l8m_3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(39270)*I*(-39*pow(cos(theta), 4) + 26*pow(cos(theta), 2) - 3)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-4} */
static double complex dYlm_dphi_l8m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(2618)*I*(-65*pow(cos(theta), 4) + 26*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-5} */
static double complex dYlm_dphi_l8m_5(const double theta, const double phi)
{
  return (15.0/128.0)*sqrt(34034)*I*(1 - 5*pow(cos(theta), 2))*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-6} */
static double complex dYlm_dphi_l8m_6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(7293)*I*(1 - 15*pow(cos(theta), 2))*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-7} */
static double complex dYlm_dphi_l8m_7(const double theta, const double phi)
{
  return -21.0/128.0*sqrt(24310)*I*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{8}^{-8} */
static double complex dYlm_dphi_l8m_8(const double theta, const double phi)
{
  return -3.0/64.0*sqrt(24310)*I*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-1} */
static double complex dYlm_dphi_l9m_1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*I*(-2431*pow(cos(theta), 8) + 4004*pow(cos(theta), 6) - 2002*pow(cos(theta), 4) + 308*pow(cos(theta), 2) - 7)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-2} */
static double complex dYlm_dphi_l9m_2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(1045)*I*(273*pow(sin(theta), 4) - 455*pow(sin(theta), 2) - 221*pow(cos(theta), 6) + 189)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-3} */
static double complex dYlm_dphi_l9m_3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(21945)*I*(-221*pow(cos(theta), 6) + 195*pow(cos(theta), 4) - 39*pow(cos(theta), 2) + 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-4} */
static double complex dYlm_dphi_l9m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(190190)*I*(-17*pow(cos(theta), 4) + 10*pow(cos(theta), 2) - 1)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-5} */
static double complex dYlm_dphi_l9m_5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(2717)*I*(-85*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 1)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-6} */
static double complex dYlm_dphi_l9m_6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(40755)*I*(3 - 17*pow(cos(theta), 2))*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-7} */
static double complex dYlm_dphi_l9m_7(const double theta, const double phi)
{
  return (21.0/512.0)*sqrt(13585)*I*(1 - 17*pow(cos(theta), 2))*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-8} */
static double complex dYlm_dphi_l9m_8(const double theta, const double phi)
{
  return -3.0/64.0*sqrt(461890)*I*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{9}^{-9} */
static double complex dYlm_dphi_l9m_9(const double theta, const double phi)
{
  return -9.0/512.0*sqrt(230945)*I*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-1} */
static double complex dYlm_dphi_l10m_1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*I*(-4199*pow(cos(theta), 8) + 7956*pow(cos(theta), 6) - 4914*pow(cos(theta), 4) + 1092*pow(cos(theta), 2) - 63)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-2} */
static double complex dYlm_dphi_l10m_2(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(770)*I*(-2730*pow(sin(theta), 4) + 5096*pow(sin(theta), 2) - 4199*pow(cos(theta), 8) + 6188*pow(cos(theta), 6) - 2373)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-3} */
static double complex dYlm_dphi_l10m_3(const double theta, const double phi)
{
  return (9.0/256.0)*sqrt(5005)*I*(-323*pow(cos(theta), 6) + 357*pow(cos(theta), 4) - 105*pow(cos(theta), 2) + 7)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-4} */
static double complex dYlm_dphi_l10m_4(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(10010)*I*(-323*pow(cos(theta), 6) + 255*pow(cos(theta), 4) - 45*pow(cos(theta), 2) + 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-5} */
static double complex dYlm_dphi_l10m_5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(1001)*I*(-323*pow(cos(theta), 4) + 170*pow(cos(theta), 2) - 15)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-6} */
static double complex dYlm_dphi_l10m_6(const double theta, const double phi)
{
  return (9.0/512.0)*sqrt(5005)*I*(-323*pow(cos(theta), 4) + 102*pow(cos(theta), 2) - 3)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-7} */
static double complex dYlm_dphi_l10m_7(const double theta, const double phi)
{
  return (21.0/512.0)*sqrt(85085)*I*(3 - 19*pow(cos(theta), 2))*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-8} */
static double complex dYlm_dphi_l10m_8(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(510510)*I*(1 - 19*pow(cos(theta), 2))*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-9} */
static double complex dYlm_dphi_l10m_9(const double theta, const double phi)
{
  return -9.0/512.0*sqrt(4849845)*I*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{10}^{-10} */
static double complex dYlm_dphi_l10m_10(const double theta, const double phi)
{
  return -5.0/512.0*sqrt(969969)*I*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-1} */
static double complex dYlm_dphi_l11m_1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*I*(-29393*pow(cos(theta), 10) + 62985*pow(cos(theta), 8) - 46410*pow(cos(theta), 6) + 13650*pow(cos(theta), 4) - 1365*pow(cos(theta), 2) + 21)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-2} */
static double complex dYlm_dphi_l11m_2(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(98670)*I*(-2142*pow(sin(theta), 4) + 3864*pow(sin(theta), 2) - 2261*pow(cos(theta), 8) + 3876*pow(cos(theta), 6) - 1743)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-3} */
static double complex dYlm_dphi_l11m_3(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(345345)*I*(-969*pow(cos(theta), 8) + 1292*pow(cos(theta), 6) - 510*pow(cos(theta), 4) + 60*pow(cos(theta), 2) - 1)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-4} */
static double complex dYlm_dphi_l11m_4(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(46046)*I*(-323*pow(cos(theta), 6) + 323*pow(cos(theta), 4) - 85*pow(cos(theta), 2) + 5)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-5} */
static double complex dYlm_dphi_l11m_5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(6578)*I*(-2261*pow(cos(theta), 6) + 1615*pow(cos(theta), 4) - 255*pow(cos(theta), 2) + 5)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-6} */
static double complex dYlm_dphi_l11m_6(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(167739)*I*(-399*pow(cos(theta), 4) + 190*pow(cos(theta), 2) - 15)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-7} */
static double complex dYlm_dphi_l11m_7(const double theta, const double phi)
{
  return (7.0/2048.0)*sqrt(1677390)*I*(-133*pow(cos(theta), 4) + 38*pow(cos(theta), 2) - 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-8} */
static double complex dYlm_dphi_l11m_8(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(31870410)*I*(1 - 7*pow(cos(theta), 2))*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-9} */
static double complex dYlm_dphi_l11m_9(const double theta, const double phi)
{
  return (9.0/2048.0)*sqrt(2124694)*I*(1 - 21*pow(cos(theta), 2))*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-10} */
static double complex dYlm_dphi_l11m_10(const double theta, const double phi)
{
  return -5.0/512.0*sqrt(22309287)*I*cexp(-10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{11}^{-11} */
static double complex dYlm_dphi_l11m_11(const double theta, const double phi)
{
  return -11.0/2048.0*sqrt(4056234)*I*cexp(-11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-1} */
static double complex dYlm_dphi_l12m_1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*I*(-52003*pow(cos(theta), 10) + 124355*pow(cos(theta), 8) - 106590*pow(cos(theta), 6) + 39270*pow(cos(theta), 4) - 5775*pow(cos(theta), 2) + 231)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-2} */
static double complex dYlm_dphi_l12m_2(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(6006)*I*(2550*pow(sin(theta), 4) - 4875*pow(sin(theta), 2) - 7429*pow(cos(theta), 10) + 14535*pow(cos(theta), 8) - 9690*pow(cos(theta), 6) + 2328)*cexp(-2*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-3} */
static double complex dYlm_dphi_l12m_3(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*I*(-7429*pow(cos(theta), 8) + 11628*pow(cos(theta), 6) - 5814*pow(cos(theta), 4) + 1020*pow(cos(theta), 2) - 45)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-4} */
static double complex dYlm_dphi_l12m_4(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*I*(-7429*pow(cos(theta), 8) + 9044*pow(cos(theta), 6) - 3230*pow(cos(theta), 4) + 340*pow(cos(theta), 2) - 5)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-5} */
static double complex dYlm_dphi_l12m_5(const double theta, const double phi)
{
  return (75.0/2048.0)*sqrt(34034)*I*(-437*pow(cos(theta), 6) + 399*pow(cos(theta), 4) - 95*pow(cos(theta), 2) + 5)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-6} */
static double complex dYlm_dphi_l12m_6(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(2431)*I*(-3059*pow(cos(theta), 6) + 1995*pow(cos(theta), 4) - 285*pow(cos(theta), 2) + 5)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-7} */
static double complex dYlm_dphi_l12m_7(const double theta, const double phi)
{
  return (35.0/2048.0)*sqrt(277134)*I*(-161*pow(cos(theta), 4) + 70*pow(cos(theta), 2) - 5)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-8} */
static double complex dYlm_dphi_l12m_8(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(277134)*I*(-161*pow(cos(theta), 4) + 42*pow(cos(theta), 2) - 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-9} */
static double complex dYlm_dphi_l12m_9(const double theta, const double phi)
{
  return (45.0/2048.0)*sqrt(646646)*I*(3 - 23*pow(cos(theta), 2))*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-10} */
static double complex dYlm_dphi_l12m_10(const double theta, const double phi)
{
  return (25.0/1024.0)*sqrt(88179)*I*(1 - 23*pow(cos(theta), 2))*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-11} */
static double complex dYlm_dphi_l12m_11(const double theta, const double phi)
{
  return -55.0/2048.0*sqrt(4056234)*I*cexp(-11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{12}^{-12} */
static double complex dYlm_dphi_l12m_12(const double theta, const double phi)
{
  return -15.0/1024.0*sqrt(676039)*I*cexp(-12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-1} */
static double complex dYlm_dphi_l13m_1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*I*(-185725*pow(cos(theta), 12) + 490314*pow(cos(theta), 10) - 479655*pow(cos(theta), 8) + 213180*pow(cos(theta), 6) - 42075*pow(cos(theta), 4) + 2970*pow(cos(theta), 2) - 33)*cexp(-I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-2} */
static double complex dYlm_dphi_l13m_2(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(2730)*I*(21318*pow(sin(theta), 4) - 39831*pow(sin(theta), 2) - 37145*pow(cos(theta), 10) + 81719*pow(cos(theta), 8) - 63954*pow(cos(theta), 6) + 18612)*cexp(-2*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-3} */
static double complex dYlm_dphi_l13m_3(const double theta, const double phi)
{
  return (9.0/8192.0)*sqrt(30030)*I*(-37145*pow(cos(theta), 10) + 66861*pow(cos(theta), 8) - 40698*pow(cos(theta), 6) + 9690*pow(cos(theta), 4) - 765*pow(cos(theta), 2) + 9)*cexp(-3*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-4} */
static double complex dYlm_dphi_l13m_4(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(51051)*I*(-10925*pow(cos(theta), 8) + 15732*pow(cos(theta), 6) - 7182*pow(cos(theta), 4) + 1140*pow(cos(theta), 2) - 45)*cexp(-4*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-5} */
static double complex dYlm_dphi_l13m_5(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(102102)*I*(-10925*pow(cos(theta), 8) + 12236*pow(cos(theta), 6) - 3990*pow(cos(theta), 4) + 380*pow(cos(theta), 2) - 5)*cexp(-5*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-6} */
static double complex dYlm_dphi_l13m_6(const double theta, const double phi)
{
  return (9.0/1024.0)*sqrt(969969)*I*(-575*pow(cos(theta), 6) + 483*pow(cos(theta), 4) - 105*pow(cos(theta), 2) + 5)*cexp(-6*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-7} */
static double complex dYlm_dphi_l13m_7(const double theta, const double phi)
{
  return (21.0/4096.0)*sqrt(692835)*I*(-805*pow(cos(theta), 6) + 483*pow(cos(theta), 4) - 63*pow(cos(theta), 2) + 1)*cexp(-7*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-8} */
static double complex dYlm_dphi_l13m_8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(9699690)*I*(-115*pow(cos(theta), 4) + 46*pow(cos(theta), 2) - 3)*cexp(-8*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-9} */
static double complex dYlm_dphi_l13m_9(const double theta, const double phi)
{
  return (27.0/4096.0)*sqrt(88179)*I*(-575*pow(cos(theta), 4) + 138*pow(cos(theta), 2) - 3)*cexp(-9*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-10} */
static double complex dYlm_dphi_l13m_10(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(2028117)*I*(3 - 25*pow(cos(theta), 2))*cexp(-10*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-11} */
static double complex dYlm_dphi_l13m_11(const double theta, const double phi)
{
  return (33.0/8192.0)*sqrt(4056234)*I*(1 - 25*pow(cos(theta), 2))*cexp(-11*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-12} */
static double complex dYlm_dphi_l13m_12(const double theta, const double phi)
{
  return -45.0/1024.0*sqrt(2028117)*I*cexp(-12*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{13}^{-13} */
static double complex dYlm_dphi_l13m_13(const double theta, const double phi)
{
  return -195.0/8192.0*sqrt(312018)*I*cexp(-13*I*phi)*pow(sin(theta), 13)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-1} */
static double complex dYlm_dphi_l14m_1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*I*(-334305*pow(cos(theta), 12) + 965770*pow(cos(theta), 10) - 1062347*pow(cos(theta), 8) + 554268*pow(cos(theta), 6) - 138567*pow(cos(theta), 4) + 14586*pow(cos(theta), 2) - 429)*cexp(-I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-2} */
static double complex dYlm_dphi_l14m_2(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(79170)*I*(334305*pow(cos(theta), 14) - 1151495*pow(cos(theta), 12) + 1552661*pow(cos(theta), 10) - 1033923*pow(cos(theta), 8) + 351747*pow(cos(theta), 6) - 56661*pow(cos(theta), 4) + 3399*pow(cos(theta), 2) - 33)*cexp(-2*I*phi)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-3} */
static double complex dYlm_dphi_l14m_3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(448630)*I*(-58995*pow(cos(theta), 10) + 120175*pow(cos(theta), 8) - 86526*pow(cos(theta), 6) + 26334*pow(cos(theta), 4) - 3135*pow(cos(theta), 2) + 99)*cexp(-3*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-4} */
static double complex dYlm_dphi_l14m_4(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2467465)*I*(-6555*pow(cos(theta), 10) + 10925*pow(cos(theta), 8) - 6118*pow(cos(theta), 6) + 1330*pow(cos(theta), 4) - 95*pow(cos(theta), 2) + 1)*cexp(-4*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-5} */
static double complex dYlm_dphi_l14m_5(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(18752734)*I*(-1725*pow(cos(theta), 8) + 2300*pow(cos(theta), 6) - 966*pow(cos(theta), 4) + 140*pow(cos(theta), 2) - 5)*cexp(-5*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-6} */
static double complex dYlm_dphi_l14m_6(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(93763670)*I*(-3105*pow(cos(theta), 8) + 3220*pow(cos(theta), 6) - 966*pow(cos(theta), 4) + 84*pow(cos(theta), 2) - 1)*cexp(-6*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-7} */
static double complex dYlm_dphi_l14m_7(const double theta, const double phi)
{
  return (7.0/4096.0)*sqrt(20092215)*I*(-1035*pow(cos(theta), 6) + 805*pow(cos(theta), 4) - 161*pow(cos(theta), 2) + 7)*cexp(-7*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-8} */
static double complex dYlm_dphi_l14m_8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(25571910)*I*(-1035*pow(cos(theta), 6) + 575*pow(cos(theta), 4) - 69*pow(cos(theta), 2) + 1)*cexp(-8*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-9} */
static double complex dYlm_dphi_l14m_9(const double theta, const double phi)
{
  return (9.0/4096.0)*sqrt(98025655)*I*(-135*pow(cos(theta), 4) + 50*pow(cos(theta), 2) - 3)*cexp(-9*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-10} */
static double complex dYlm_dphi_l14m_10(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*I*(-225*pow(cos(theta), 4) + 50*pow(cos(theta), 2) - 1)*cexp(-10*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-11} */
static double complex dYlm_dphi_l14m_11(const double theta, const double phi)
{
  return (55.0/8192.0)*sqrt(117630786)*I*(1 - 9*pow(cos(theta), 2))*cexp(-11*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-12} */
static double complex dYlm_dphi_l14m_12(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(1508087)*I*(1 - 27*pow(cos(theta), 2))*cexp(-12*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-13} */
static double complex dYlm_dphi_l14m_13(const double theta, const double phi)
{
  return -195.0/8192.0*sqrt(9048522)*I*cexp(-13*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* dY_dphi_{14}^{-14} */
static double complex dYlm_dphi_l14m_14(const double theta, const double phi)
{
  return -105.0/8192.0*sqrt(1292646)*I*cexp(-14*I*phi)*pow(sin(theta), 14)/sqrt(M_PI);
}


/* dY_dtheta_{0}^{0} */
static double complex dYlm_dtheta_l0m0(const double theta, const double phi)
{
  UNUSED(phi);
  UNUSED(theta);
  return 0;
}


/* dY_dtheta_{1}^{0} */
static double complex dYlm_dtheta_l1m0(const double theta, const double phi)
{
  UNUSED(phi);
  return -1.0/2.0*sqrt(3)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{1}^{1} */
static double complex dYlm_dtheta_l1m1(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(6)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{2}^{0} */
static double complex dYlm_dtheta_l2m0(const double theta, const double phi)
{
  UNUSED(phi);
  return -3.0/4.0*sqrt(5)*sin(2*theta)/sqrt(M_PI);
}


/* dY_dtheta_{2}^{1} */
static double complex dYlm_dtheta_l2m1(const double theta, const double phi)
{
  return -1.0/4.0*sqrt(30)*cexp(I*phi)*cos(2*theta)/sqrt(M_PI);
}


/* dY_dtheta_{2}^{2} */
static double complex dYlm_dtheta_l2m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(2*I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{0} */
static double complex dYlm_dtheta_l3m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (3.0/4.0)*sqrt(7)*(1 - 5*pow(cos(theta), 2))*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{1} */
static double complex dYlm_dtheta_l3m1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*(15*pow(sin(theta), 2) - 4)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{2} */
static double complex dYlm_dtheta_l3m2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(210)*(3*pow(cos(theta), 2) - 1)*cexp(2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{3} */
static double complex dYlm_dtheta_l3m3(const double theta, const double phi)
{
  return -3.0/8.0*sqrt(35)*cexp(3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{0} */
static double complex dYlm_dtheta_l4m0(const double theta, const double phi)
{
  UNUSED(phi);
  return -1.0/32.0*(30*sin(2*theta) + 105*sin(4*theta))/sqrt(M_PI);
}


/* dY_dtheta_{4}^{1} */
static double complex dYlm_dtheta_l4m1(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(5)*(-28*pow(sin(theta), 4) + 29*pow(sin(theta), 2) - 4)*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{2} */
static double complex dYlm_dtheta_l4m2(const double theta, const double phi)
{
  return sqrt(10)*(-3.0/16.0*sin(2*theta) + (21.0/32.0)*sin(4*theta))*cexp(2*I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{3} */
static double complex dYlm_dtheta_l4m3(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(35)*(4*pow(sin(theta), 2) - 3)*cexp(3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{4} */
static double complex dYlm_dtheta_l4m4(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(70)*cexp(4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{0} */
static double complex dYlm_dtheta_l5m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (15.0/16.0)*sqrt(11)*(-21*pow(cos(theta), 4) + 14*pow(cos(theta), 2) - 1)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{1} */
static double complex dYlm_dtheta_l5m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*(-105*pow(sin(theta), 4) + 84*pow(sin(theta), 2) - 8)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{2} */
static double complex dYlm_dtheta_l5m2(const double theta, const double phi)
{
  return (1.0/16.0)*sqrt(2310)*(15*pow(cos(theta), 4) - 12*pow(cos(theta), 2) + 1)*cexp(2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{3} */
static double complex dYlm_dtheta_l5m3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(385)*(15*pow(sin(theta), 2) - 8)*cexp(3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{4} */
static double complex dYlm_dtheta_l5m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(5*cos(2*theta) + 3)*cexp(4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{5} */
static double complex dYlm_dtheta_l5m5(const double theta, const double phi)
{
  return -15.0/32.0*sqrt(77)*cexp(5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{0} */
static double complex dYlm_dtheta_l6m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (21.0/16.0)*sqrt(13)*(-33*pow(cos(theta), 4) + 30*pow(cos(theta), 2) - 5)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{1} */
static double complex dYlm_dtheta_l6m1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*(165*pow(sin(theta), 6) - 210*pow(sin(theta), 4) + 25*pow(sin(theta), 2) - 33*pow(cos(theta), 6) + 25)*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{2} */
static double complex dYlm_dtheta_l6m2(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*(99*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 19)*cexp(2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{3} */
static double complex dYlm_dtheta_l6m3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1365)*(-22*pow(sin(theta), 4) + 29*pow(sin(theta), 2) - 8)*cexp(3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{4} */
static double complex dYlm_dtheta_l6m4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(182)*(20 - 33*pow(sin(theta), 2))*cexp(4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{5} */
static double complex dYlm_dtheta_l6m5(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1001)*(6*pow(sin(theta), 2) - 5)*cexp(5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{6} */
static double complex dYlm_dtheta_l6m6(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(3003)*cexp(6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{0} */
static double complex dYlm_dtheta_l7m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (7.0/32.0)*sqrt(15)*(-429*pow(cos(theta), 6) + 495*pow(cos(theta), 4) - 135*pow(cos(theta), 2) + 5)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{1} */
static double complex dYlm_dtheta_l7m1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*(2574*pow(sin(theta), 6) - 2673*pow(sin(theta), 4) + 9*pow(sin(theta), 2) - 429*pow(cos(theta), 6) + 365)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{2} */
static double complex dYlm_dtheta_l7m2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(35)*(1001*pow(cos(theta), 6) - 1265*pow(cos(theta), 4) + 375*pow(cos(theta), 2) - 15)*cexp(2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{3} */
static double complex dYlm_dtheta_l7m3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(70)*(-1001*pow(sin(theta), 4) + 1100*pow(sin(theta), 2) - 240)*cexp(3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{4} */
static double complex dYlm_dtheta_l7m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(91*pow(sin(theta), 4) - 128*pow(sin(theta), 2) + 40)*cexp(4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{5} */
static double complex dYlm_dtheta_l7m5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(770)*(91*pow(sin(theta), 2) - 60)*cexp(5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{6} */
static double complex dYlm_dtheta_l7m6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(5005)*(6 - 7*pow(sin(theta), 2))*cexp(6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{7} */
static double complex dYlm_dtheta_l7m7(const double theta, const double phi)
{
  return -21.0/128.0*sqrt(1430)*cexp(7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{0} */
static double complex dYlm_dtheta_l8m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (9.0/32.0)*sqrt(17)*(-715*pow(cos(theta), 6) + 1001*pow(cos(theta), 4) - 385*pow(cos(theta), 2) + 35)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{1} */
static double complex dYlm_dtheta_l8m1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*(-5005*pow(sin(theta), 6) + 8470*pow(sin(theta), 4) - 3150*pow(sin(theta), 2) - 5720*pow(cos(theta), 8) + 6006*pow(cos(theta), 6) - 350)*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{2} */
static double complex dYlm_dtheta_l8m2(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(595)*(286*pow(cos(theta), 6) - 429*pow(cos(theta), 4) + 176*pow(cos(theta), 2) - 17)*cexp(2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{3} */
static double complex dYlm_dtheta_l8m3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(39270)*(65*pow(sin(theta), 6) - 78*pow(sin(theta), 4) - 9*pow(sin(theta), 2) - 39*pow(cos(theta), 6) + 23)*cexp(3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{4} */
static double complex dYlm_dtheta_l8m4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(2618)*(130*pow(sin(theta), 4) - 156*pow(sin(theta), 2) + 40)*cexp(4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{5} */
static double complex dYlm_dtheta_l8m5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34034)*(-40*pow(sin(theta), 4) + 59*pow(sin(theta), 2) - 20)*cexp(5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{6} */
static double complex dYlm_dtheta_l8m6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(7293)*(14 - 20*pow(sin(theta), 2))*cexp(6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{7} */
static double complex dYlm_dtheta_l8m7(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(24310)*(8*pow(sin(theta), 2) - 7)*cexp(7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{8} */
static double complex dYlm_dtheta_l8m8(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(24310)*cexp(8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{0} */
static double complex dYlm_dtheta_l9m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (45.0/256.0)*sqrt(19)*(-2431*pow(cos(theta), 8) + 4004*pow(cos(theta), 6) - 2002*pow(cos(theta), 4) + 308*pow(cos(theta), 2) - 7)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{1} */
static double complex dYlm_dtheta_l9m1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*(-24024*pow(sin(theta), 6) + 38038*pow(sin(theta), 4) - 12936*pow(sin(theta), 2) - 21879*pow(cos(theta), 8) + 23452*pow(cos(theta), 6) - 1701)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{2} */
static double complex dYlm_dtheta_l9m2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1045)*(1989*pow(cos(theta), 8) - 3458*pow(cos(theta), 6) + 1820*pow(cos(theta), 4) - 294*pow(cos(theta), 2) + 7)*cexp(2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{3} */
static double complex dYlm_dtheta_l9m3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(21945)*(442*pow(sin(theta), 6) - 429*pow(sin(theta), 4) - 143*pow(sin(theta), 2) - 221*pow(cos(theta), 6) + 157)*cexp(3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{4} */
static double complex dYlm_dtheta_l9m4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(190190)*(-85*pow(sin(theta), 6) + 100*pow(sin(theta), 4) + 20*pow(sin(theta), 2) + 68*pow(cos(theta), 6) - 36)*cexp(4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{5} */
static double complex dYlm_dtheta_l9m5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(2717)*(-153*pow(sin(theta), 4) + 196*pow(sin(theta), 2) - 56)*cexp(5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{6} */
static double complex dYlm_dtheta_l9m6(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(40755)*(51*pow(sin(theta), 4) - 78*pow(sin(theta), 2) + 28)*cexp(6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{7} */
static double complex dYlm_dtheta_l9m7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(13585)*(153*pow(sin(theta), 2) - 112)*cexp(7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{8} */
static double complex dYlm_dtheta_l9m8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(461890)*(8 - 9*pow(sin(theta), 2))*cexp(8*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{9} */
static double complex dYlm_dtheta_l9m9(const double theta, const double phi)
{
  return -9.0/512.0*sqrt(230945)*cexp(9*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{0} */
static double complex dYlm_dtheta_l10m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (55.0/256.0)*sqrt(21)*(-4199*pow(cos(theta), 8) + 7956*pow(cos(theta), 6) - 4914*pow(cos(theta), 4) + 1092*pow(cos(theta), 2) - 63)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{1} */
static double complex dYlm_dtheta_l10m1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*(24570*pow(sin(theta), 6) - 44772*pow(sin(theta), 4) + 19236*pow(sin(theta), 2) - 41990*pow(cos(theta), 10) + 101439*pow(cos(theta), 8) - 60606*pow(cos(theta), 6) + 1029)*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{2} */
static double complex dYlm_dtheta_l10m2(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(770)*(20995*pow(cos(theta), 8) - 41548*pow(cos(theta), 6) + 26754*pow(cos(theta), 4) - 6188*pow(cos(theta), 2) + 371)*cexp(2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{3} */
static double complex dYlm_dtheta_l10m3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(5005)*(-1785*pow(sin(theta), 6) + 2940*pow(sin(theta), 4) - 868*pow(sin(theta), 2) - 3230*pow(cos(theta), 8) + 3332*pow(cos(theta), 6) - 294)*cexp(3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{4} */
static double complex dYlm_dtheta_l10m4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(10010)*(-969*pow(sin(theta), 6) + 918*pow(sin(theta), 4) + 426*pow(sin(theta), 2) + 646*pow(cos(theta), 6) - 422)*cexp(4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{5} */
static double complex dYlm_dtheta_l10m5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(1001)*(323*pow(sin(theta), 6) - 374*pow(sin(theta), 4) - 101*pow(sin(theta), 2) - 323*pow(cos(theta), 6) + 155)*cexp(5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{6} */
static double complex dYlm_dtheta_l10m6(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(5005)*(1615*pow(sin(theta), 4) - 2176*pow(sin(theta), 2) + 672)*cexp(6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{7} */
static double complex dYlm_dtheta_l10m7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(85085)*(-190*pow(sin(theta), 4) + 299*pow(sin(theta), 2) - 112)*cexp(7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{8} */
static double complex dYlm_dtheta_l10m8(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(510510)*(72 - 95*pow(sin(theta), 2))*cexp(8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{9} */
static double complex dYlm_dtheta_l10m9(const double theta, const double phi)
{
  return -1.0/512.0*sqrt(4849845)*(5*cos(2*theta) + 4)*cexp(9*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{10} */
static double complex dYlm_dtheta_l10m10(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(969969)*cexp(10*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{0} */
static double complex dYlm_dtheta_l11m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (33.0/512.0)*sqrt(23)*(-29393*pow(cos(theta), 10) + 62985*pow(cos(theta), 8) - 46410*pow(cos(theta), 6) + 13650*pow(cos(theta), 4) - 1365*pow(cos(theta), 2) + 21)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{1} */
static double complex dYlm_dtheta_l11m1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*(278460*pow(sin(theta), 6) - 488670*pow(sin(theta), 4) + 200655*pow(sin(theta), 2) - 323323*pow(cos(theta), 10) + 860795*pow(cos(theta), 8) - 550290*pow(cos(theta), 6) + 12306)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{2} */
static double complex dYlm_dtheta_l11m2(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(98670)*(24871*pow(cos(theta), 10) - 55233*pow(cos(theta), 8) + 42126*pow(cos(theta), 6) - 12810*pow(cos(theta), 4) + 1323*pow(cos(theta), 2) - 21)*cexp(2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{3} */
static double complex dYlm_dtheta_l11m3(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(345345)*(-2584*pow(sin(theta), 6) + 3978*pow(sin(theta), 4) - 984*pow(sin(theta), 2) - 3553*pow(cos(theta), 8) + 3876*pow(cos(theta), 6) - 451)*cexp(3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{4} */
static double complex dYlm_dtheta_l11m4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(46046)*(1615*pow(sin(theta), 6) - 2635*pow(sin(theta), 4) + 705*pow(sin(theta), 2) + 3553*pow(cos(theta), 8) - 3553*pow(cos(theta), 6) + 320)*cexp(4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{5} */
static double complex dYlm_dtheta_l11m5(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(6578)*(13566*pow(sin(theta), 6) - 12597*pow(sin(theta), 4) - 7259*pow(sin(theta), 2) - 11305*pow(cos(theta), 6) + 6825)*cexp(5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{6} */
static double complex dYlm_dtheta_l11m6(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(167739)*(-665*pow(sin(theta), 6) + 760*pow(sin(theta), 4) + 250*pow(sin(theta), 2) + 798*pow(cos(theta), 6) - 350)*cexp(6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{7} */
static double complex dYlm_dtheta_l11m7(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(1677390)*(-1463*pow(sin(theta), 4) + 2052*pow(sin(theta), 2) - 672)*cexp(7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{8} */
static double complex dYlm_dtheta_l11m8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(31870410)*(77*pow(sin(theta), 4) - 124*pow(sin(theta), 2) + 48)*cexp(8*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{9} */
static double complex dYlm_dtheta_l11m9(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2124694)*(77*pow(sin(theta), 2) - 60)*cexp(9*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{10} */
static double complex dYlm_dtheta_l11m10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(22309287)*(10 - 11*pow(sin(theta), 2))*cexp(10*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{11} */
static double complex dYlm_dtheta_l11m11(const double theta, const double phi)
{
  return -11.0/2048.0*sqrt(4056234)*cexp(11*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{0} */
static double complex dYlm_dtheta_l12m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (195.0/512.0)*(-52003*pow(cos(theta), 10) + 124355*pow(cos(theta), 8) - 106590*pow(cos(theta), 6) + 39270*pow(cos(theta), 4) - 5775*pow(cos(theta), 2) + 231)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{1} */
static double complex dYlm_dtheta_l12m1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*(-196350*pow(sin(theta), 6) + 369600*pow(sin(theta), 4) - 167937*pow(sin(theta), 2) - 624036*pow(cos(theta), 12) + 1815583*pow(cos(theta), 10) - 1971915*pow(cos(theta), 8) + 785400*pow(cos(theta), 6) - 5544)*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{2} */
static double complex dYlm_dtheta_l12m2(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(6006)*(22287*pow(cos(theta), 10) - 54910*pow(cos(theta), 8) + 48450*pow(cos(theta), 6) - 18360*pow(cos(theta), 4) + 2775*pow(cos(theta), 2) - 114)*cexp(2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{3} */
static double complex dYlm_dtheta_l12m3(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*(9690*pow(sin(theta), 6) - 17340*pow(sin(theta), 4) + 6690*pow(sin(theta), 2) - 29716*pow(cos(theta), 10) + 61047*pow(cos(theta), 8) - 32946*pow(cos(theta), 6) + 975)*cexp(3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{4} */
static double complex dYlm_dtheta_l12m4(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*(13566*pow(sin(theta), 6) - 20672*pow(sin(theta), 4) + 4386*pow(sin(theta), 2) + 22287*pow(cos(theta), 8) - 23902*pow(cos(theta), 6) + 2895)*cexp(4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{5} */
static double complex dYlm_dtheta_l12m5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(34034)*(-1995*pow(sin(theta), 6) + 3230*pow(sin(theta), 4) - 790*pow(sin(theta), 2) - 5244*pow(cos(theta), 8) + 5054*pow(cos(theta), 6) - 450)*cexp(5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{6} */
static double complex dYlm_dtheta_l12m6(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(2431)*(-3059*pow(sin(theta), 6) + 2793*pow(sin(theta), 4) + 1881*pow(sin(theta), 2) + 3059*pow(cos(theta), 6) - 1715)*cexp(6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{7} */
static double complex dYlm_dtheta_l12m7(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(277134)*(805*pow(sin(theta), 6) - 910*pow(sin(theta), 4) - 345*pow(sin(theta), 2) - 1127*pow(cos(theta), 6) + 455)*cexp(7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{8} */
static double complex dYlm_dtheta_l12m8(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(277134)*(483*pow(sin(theta), 4) - 700*pow(sin(theta), 2) + 240)*cexp(8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{9} */
static double complex dYlm_dtheta_l12m9(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(646646)*(-92*pow(sin(theta), 4) + 151*pow(sin(theta), 2) - 60)*cexp(9*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{10} */
static double complex dYlm_dtheta_l12m10(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(88179)*(110 - 138*pow(sin(theta), 2))*cexp(10*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{11} */
static double complex dYlm_dtheta_l12m11(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(4056234)*(12*pow(sin(theta), 2) - 11)*cexp(11*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{12} */
static double complex dYlm_dtheta_l12m12(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(676039)*cexp(12*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{0} */
static double complex dYlm_dtheta_l13m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (273.0/2048.0)*sqrt(3)*(-185725*pow(cos(theta), 12) + 490314*pow(cos(theta), 10) - 479655*pow(cos(theta), 8) + 213180*pow(cos(theta), 6) - 42075*pow(cos(theta), 4) + 2970*pow(cos(theta), 2) - 33)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{1} */
static double complex dYlm_dtheta_l13m1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*(-1279080*pow(sin(theta), 6) + 2347785*pow(sin(theta), 4) - 1035540*pow(sin(theta), 2) - 2414425*pow(cos(theta), 12) + 7622154*pow(cos(theta), 10) - 9220035*pow(cos(theta), 8) + 4050420*pow(cos(theta), 6) - 39138)*cexp(I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{2} */
static double complex dYlm_dtheta_l13m2(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2730)*(482885*pow(cos(theta), 12) - 1307504*pow(cos(theta), 10) + 1311057*pow(cos(theta), 8) - 596904*pow(cos(theta), 6) + 120615*pow(cos(theta), 4) - 8712*pow(cos(theta), 2) + 99)*cexp(2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{3} */
static double complex dYlm_dtheta_l13m3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(30030)*(244188*pow(sin(theta), 6) - 420546*pow(sin(theta), 4) + 151113*pow(sin(theta), 2) - 482885*pow(cos(theta), 10) + 1106921*pow(cos(theta), 8) - 656982*pow(cos(theta), 6) + 26802)*cexp(3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{4} */
static double complex dYlm_dtheta_l13m4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(51051)*(-35910*pow(sin(theta), 6) + 63840*pow(sin(theta), 4) - 23595*pow(sin(theta), 2) + 142025*pow(cos(theta), 10) - 271377*pow(cos(theta), 8) + 138852*pow(cos(theta), 6) - 4380)*cexp(4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{5} */
static double complex dYlm_dtheta_l13m5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(102102)*(-73416*pow(sin(theta), 6) + 110922*pow(sin(theta), 4) - 20216*pow(sin(theta), 2) - 142025*pow(cos(theta), 8) + 148580*pow(cos(theta), 6) - 18075)*cexp(5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{6} */
static double complex dYlm_dtheta_l13m6(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(969969)*(2415*pow(sin(theta), 6) - 3885*pow(sin(theta), 4) + 875*pow(sin(theta), 2) + 7475*pow(cos(theta), 8) - 6923*pow(cos(theta), 6) + 600)*cexp(6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{7} */
static double complex dYlm_dtheta_l13m7(const double theta, const double phi)
{
  return (21.0/4096.0)*sqrt(692835)*(690*pow(sin(theta), 6) - 621*pow(sin(theta), 4) - 471*pow(sin(theta), 2) - 805*pow(cos(theta), 6) + 421)*cexp(7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{8} */
static double complex dYlm_dtheta_l13m8(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(9699690)*(-575*pow(sin(theta), 6) + 644*pow(sin(theta), 4) + 272*pow(sin(theta), 2) + 920*pow(cos(theta), 6) - 344)*cexp(8*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{9} */
static double complex dYlm_dtheta_l13m9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(88179)*(-7475*pow(sin(theta), 4) + 11132*pow(sin(theta), 2) - 3960)*cexp(9*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{10} */
static double complex dYlm_dtheta_l13m10(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2028117)*(325*pow(sin(theta), 4) - 542*pow(sin(theta), 2) + 220)*cexp(10*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{11} */
static double complex dYlm_dtheta_l13m11(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(4056234)*(325*pow(sin(theta), 2) - 264)*cexp(11*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{12} */
static double complex dYlm_dtheta_l13m12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(2028117)*(12 - 13*pow(sin(theta), 2))*cexp(12*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{13} */
static double complex dYlm_dtheta_l13m13(const double theta, const double phi)
{
  return -195.0/8192.0*sqrt(312018)*cexp(13*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{0} */
static double complex dYlm_dtheta_l14m0(const double theta, const double phi)
{
  UNUSED(phi);
  return (105.0/2048.0)*sqrt(29)*(-334305*pow(cos(theta), 12) + 965770*pow(cos(theta), 10) - 1062347*pow(cos(theta), 8) + 554268*pow(cos(theta), 6) - 138567*pow(cos(theta), 4) + 14586*pow(cos(theta), 2) - 429)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{1} */
static double complex dYlm_dtheta_l14m1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*(692835*pow(sin(theta), 6) - 1327326*pow(sin(theta), 4) + 620763*pow(sin(theta), 2) - 4680270*pow(cos(theta), 14) + 15935205*pow(cos(theta), 12) - 21246940*pow(cos(theta), 10) + 13995267*pow(cos(theta), 8) - 4018443*pow(cos(theta), 6) + 14157)*cexp(I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{2} */
static double complex dYlm_dtheta_l14m2(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(79170)*(2340135*pow(cos(theta), 12) - 6908970*pow(cos(theta), 10) + 7763305*pow(cos(theta), 8) - 4135692*pow(cos(theta), 6) + 1055241*pow(cos(theta), 4) - 113322*pow(cos(theta), 2) + 3399)*cexp(2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{3} */
static double complex dYlm_dtheta_l14m3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(448630)*(-43890*pow(sin(theta), 6) + 81510*pow(sin(theta), 4) - 34617*pow(sin(theta), 2) - 275310*pow(cos(theta), 12) + 697015*pow(cos(theta), 10) - 648945*pow(cos(theta), 8) + 228228*pow(cos(theta), 6) - 3036)*cexp(3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{4} */
static double complex dYlm_dtheta_l14m4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(2467465)*(-18354*pow(sin(theta), 6) + 31388*pow(sin(theta), 4) - 10659*pow(sin(theta), 2) + 45885*pow(cos(theta), 10) - 98325*pow(cos(theta), 8) + 55936*pow(cos(theta), 6) - 2472)*cexp(4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{5} */
static double complex dYlm_dtheta_l14m5(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(18752734)*(966*pow(sin(theta), 6) - 1708*pow(sin(theta), 4) + 608*pow(sin(theta), 2) - 4830*pow(cos(theta), 10) + 8625*pow(cos(theta), 8) - 4186*pow(cos(theta), 6) + 135)*cexp(5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{6} */
static double complex dYlm_dtheta_l14m6(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(93763670)*(3220*pow(sin(theta), 6) - 4830*pow(sin(theta), 4) + 756*pow(sin(theta), 2) + 7245*pow(cos(theta), 8) - 7360*pow(cos(theta), 6) + 883)*cexp(6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{7} */
static double complex dYlm_dtheta_l14m7(const double theta, const double phi)
{
  return (7.0/4096.0)*sqrt(20092215)*(-575*pow(sin(theta), 6) + 920*pow(sin(theta), 4) - 192*pow(sin(theta), 2) - 2070*pow(cos(theta), 8) + 1840*pow(cos(theta), 6) - 154)*cexp(7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{8} */
static double complex dYlm_dtheta_l14m8(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(25571910)*(-3105*pow(sin(theta), 6) + 2760*pow(sin(theta), 4) + 2300*pow(sin(theta), 2) + 4140*pow(cos(theta), 6) - 2028)*cexp(8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{9} */
static double complex dYlm_dtheta_l14m9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(98025655)*(225*pow(sin(theta), 6) - 250*pow(sin(theta), 4) - 115*pow(sin(theta), 2) - 405*pow(cos(theta), 6) + 141)*cexp(9*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{10} */
static double complex dYlm_dtheta_l14m10(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(315*pow(sin(theta), 4) - 480*pow(sin(theta), 2) + 176)*cexp(10*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{11} */
static double complex dYlm_dtheta_l14m11(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(-126*pow(sin(theta), 4) + 213*pow(sin(theta), 2) - 88)*cexp(11*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{12} */
static double complex dYlm_dtheta_l14m12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(1508087)*(52 - 63*pow(sin(theta), 2))*cexp(12*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{13} */
static double complex dYlm_dtheta_l14m13(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(9048522)*(14*pow(sin(theta), 2) - 13)*cexp(13*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{14} */
static double complex dYlm_dtheta_l14m14(const double theta, const double phi)
{
  return (105.0/8192.0)*sqrt(1292646)*cexp(14*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{1}^{-1} */
static double complex dYlm_dtheta_l1m_1(const double theta, const double phi)
{
  return (1.0/4.0)*sqrt(6)*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{2}^{-1} */
static double complex dYlm_dtheta_l2m_1(const double theta, const double phi)
{
  return (1.0/4.0)*sqrt(30)*cexp(-I*phi)*cos(2*theta)/sqrt(M_PI);
}


/* dY_dtheta_{2}^{-2} */
static double complex dYlm_dtheta_l2m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(30)*cexp(-2*I*phi)*sin(2*theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{-1} */
static double complex dYlm_dtheta_l3m_1(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(21)*(4 - 15*pow(sin(theta), 2))*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{-2} */
static double complex dYlm_dtheta_l3m_2(const double theta, const double phi)
{
  return (1.0/8.0)*sqrt(210)*(3*pow(cos(theta), 2) - 1)*cexp(-2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{3}^{-3} */
static double complex dYlm_dtheta_l3m_3(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(35)*cexp(-3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{-1} */
static double complex dYlm_dtheta_l4m_1(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(5)*(28*pow(sin(theta), 4) - 29*pow(sin(theta), 2) + 4)*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{-2} */
static double complex dYlm_dtheta_l4m_2(const double theta, const double phi)
{
  return sqrt(10)*(-3.0/16.0*sin(2*theta) + (21.0/32.0)*sin(4*theta))*cexp(-2*I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{-3} */
static double complex dYlm_dtheta_l4m_3(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(35)*(3 - 4*pow(sin(theta), 2))*cexp(-3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{4}^{-4} */
static double complex dYlm_dtheta_l4m_4(const double theta, const double phi)
{
  return (3.0/8.0)*sqrt(70)*cexp(-4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{-1} */
static double complex dYlm_dtheta_l5m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(330)*(105*pow(sin(theta), 4) - 84*pow(sin(theta), 2) + 8)*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{-2} */
static double complex dYlm_dtheta_l5m_2(const double theta, const double phi)
{
  return (1.0/16.0)*sqrt(2310)*(15*pow(cos(theta), 4) - 12*pow(cos(theta), 2) + 1)*cexp(-2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{-3} */
static double complex dYlm_dtheta_l5m_3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(385)*(8 - 15*pow(sin(theta), 2))*cexp(-3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{-4} */
static double complex dYlm_dtheta_l5m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(5*cos(2*theta) + 3)*cexp(-4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{5}^{-5} */
static double complex dYlm_dtheta_l5m_5(const double theta, const double phi)
{
  return (15.0/32.0)*sqrt(77)*cexp(-5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{-1} */
static double complex dYlm_dtheta_l6m_1(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(546)*(-165*pow(sin(theta), 6) + 210*pow(sin(theta), 4) - 25*pow(sin(theta), 2) + 33*pow(cos(theta), 6) - 25)*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{-2} */
static double complex dYlm_dtheta_l6m_2(const double theta, const double phi)
{
  return (1.0/32.0)*sqrt(1365)*(99*pow(cos(theta), 4) - 102*pow(cos(theta), 2) + 19)*cexp(-2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{-3} */
static double complex dYlm_dtheta_l6m_3(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1365)*(22*pow(sin(theta), 4) - 29*pow(sin(theta), 2) + 8)*cexp(-3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{-4} */
static double complex dYlm_dtheta_l6m_4(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(182)*(20 - 33*pow(sin(theta), 2))*cexp(-4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{-5} */
static double complex dYlm_dtheta_l6m_5(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(1001)*(5 - 6*pow(sin(theta), 2))*cexp(-5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{6}^{-6} */
static double complex dYlm_dtheta_l6m_6(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(3003)*cexp(-6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-1} */
static double complex dYlm_dtheta_l7m_1(const double theta, const double phi)
{
  return (1.0/128.0)*sqrt(210)*(-2574*pow(sin(theta), 6) + 2673*pow(sin(theta), 4) - 9*pow(sin(theta), 2) + 429*pow(cos(theta), 6) - 365)*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-2} */
static double complex dYlm_dtheta_l7m_2(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(35)*(1001*pow(cos(theta), 6) - 1265*pow(cos(theta), 4) + 375*pow(cos(theta), 2) - 15)*cexp(-2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-3} */
static double complex dYlm_dtheta_l7m_3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(70)*(1001*pow(sin(theta), 4) - 1100*pow(sin(theta), 2) + 240)*cexp(-3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-4} */
static double complex dYlm_dtheta_l7m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(770)*(91*pow(sin(theta), 4) - 128*pow(sin(theta), 2) + 40)*cexp(-4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-5} */
static double complex dYlm_dtheta_l7m_5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(770)*(60 - 91*pow(sin(theta), 2))*cexp(-5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-6} */
static double complex dYlm_dtheta_l7m_6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(5005)*(6 - 7*pow(sin(theta), 2))*cexp(-6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{7}^{-7} */
static double complex dYlm_dtheta_l7m_7(const double theta, const double phi)
{
  return (21.0/128.0)*sqrt(1430)*cexp(-7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-1} */
static double complex dYlm_dtheta_l8m_1(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34)*(5005*pow(sin(theta), 6) - 8470*pow(sin(theta), 4) + 3150*pow(sin(theta), 2) + 5720*pow(cos(theta), 8) - 6006*pow(cos(theta), 6) + 350)*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-2} */
static double complex dYlm_dtheta_l8m_2(const double theta, const double phi)
{
  return (3.0/32.0)*sqrt(595)*(286*pow(cos(theta), 6) - 429*pow(cos(theta), 4) + 176*pow(cos(theta), 2) - 17)*cexp(-2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-3} */
static double complex dYlm_dtheta_l8m_3(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(39270)*(-65*pow(sin(theta), 6) + 78*pow(sin(theta), 4) + 9*pow(sin(theta), 2) + 39*pow(cos(theta), 6) - 23)*cexp(-3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-4} */
static double complex dYlm_dtheta_l8m_4(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(2618)*(130*pow(sin(theta), 4) - 156*pow(sin(theta), 2) + 40)*cexp(-4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-5} */
static double complex dYlm_dtheta_l8m_5(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(34034)*(40*pow(sin(theta), 4) - 59*pow(sin(theta), 2) + 20)*cexp(-5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-6} */
static double complex dYlm_dtheta_l8m_6(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(7293)*(14 - 20*pow(sin(theta), 2))*cexp(-6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-7} */
static double complex dYlm_dtheta_l8m_7(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(24310)*(7 - 8*pow(sin(theta), 2))*cexp(-7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{8}^{-8} */
static double complex dYlm_dtheta_l8m_8(const double theta, const double phi)
{
  return (3.0/64.0)*sqrt(24310)*cexp(-8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-1} */
static double complex dYlm_dtheta_l9m_1(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(190)*(24024*pow(sin(theta), 6) - 38038*pow(sin(theta), 4) + 12936*pow(sin(theta), 2) + 21879*pow(cos(theta), 8) - 23452*pow(cos(theta), 6) + 1701)*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-2} */
static double complex dYlm_dtheta_l9m_2(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(1045)*(1989*pow(cos(theta), 8) - 3458*pow(cos(theta), 6) + 1820*pow(cos(theta), 4) - 294*pow(cos(theta), 2) + 7)*cexp(-2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-3} */
static double complex dYlm_dtheta_l9m_3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(21945)*(-442*pow(sin(theta), 6) + 429*pow(sin(theta), 4) + 143*pow(sin(theta), 2) + 221*pow(cos(theta), 6) - 157)*cexp(-3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-4} */
static double complex dYlm_dtheta_l9m_4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(190190)*(-85*pow(sin(theta), 6) + 100*pow(sin(theta), 4) + 20*pow(sin(theta), 2) + 68*pow(cos(theta), 6) - 36)*cexp(-4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-5} */
static double complex dYlm_dtheta_l9m_5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(2717)*(153*pow(sin(theta), 4) - 196*pow(sin(theta), 2) + 56)*cexp(-5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-6} */
static double complex dYlm_dtheta_l9m_6(const double theta, const double phi)
{
  return (3.0/128.0)*sqrt(40755)*(51*pow(sin(theta), 4) - 78*pow(sin(theta), 2) + 28)*cexp(-6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-7} */
static double complex dYlm_dtheta_l9m_7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(13585)*(112 - 153*pow(sin(theta), 2))*cexp(-7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-8} */
static double complex dYlm_dtheta_l9m_8(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(461890)*(8 - 9*pow(sin(theta), 2))*cexp(-8*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dtheta_{9}^{-9} */
static double complex dYlm_dtheta_l9m_9(const double theta, const double phi)
{
  return (9.0/512.0)*sqrt(230945)*cexp(-9*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-1} */
static double complex dYlm_dtheta_l10m_1(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(2310)*(-24570*pow(sin(theta), 6) + 44772*pow(sin(theta), 4) - 19236*pow(sin(theta), 2) + 41990*pow(cos(theta), 10) - 101439*pow(cos(theta), 8) + 60606*pow(cos(theta), 6) - 1029)*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-2} */
static double complex dYlm_dtheta_l10m_2(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(770)*(20995*pow(cos(theta), 8) - 41548*pow(cos(theta), 6) + 26754*pow(cos(theta), 4) - 6188*pow(cos(theta), 2) + 371)*cexp(-2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-3} */
static double complex dYlm_dtheta_l10m_3(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(5005)*(1785*pow(sin(theta), 6) - 2940*pow(sin(theta), 4) + 868*pow(sin(theta), 2) + 3230*pow(cos(theta), 8) - 3332*pow(cos(theta), 6) + 294)*cexp(-3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-4} */
static double complex dYlm_dtheta_l10m_4(const double theta, const double phi)
{
  return (3.0/256.0)*sqrt(10010)*(-969*pow(sin(theta), 6) + 918*pow(sin(theta), 4) + 426*pow(sin(theta), 2) + 646*pow(cos(theta), 6) - 422)*cexp(-4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-5} */
static double complex dYlm_dtheta_l10m_5(const double theta, const double phi)
{
  return (15.0/256.0)*sqrt(1001)*(-323*pow(sin(theta), 6) + 374*pow(sin(theta), 4) + 101*pow(sin(theta), 2) + 323*pow(cos(theta), 6) - 155)*cexp(-5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-6} */
static double complex dYlm_dtheta_l10m_6(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(5005)*(1615*pow(sin(theta), 4) - 2176*pow(sin(theta), 2) + 672)*cexp(-6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-7} */
static double complex dYlm_dtheta_l10m_7(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(85085)*(190*pow(sin(theta), 4) - 299*pow(sin(theta), 2) + 112)*cexp(-7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-8} */
static double complex dYlm_dtheta_l10m_8(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(510510)*(72 - 95*pow(sin(theta), 2))*cexp(-8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-9} */
static double complex dYlm_dtheta_l10m_9(const double theta, const double phi)
{
  return (1.0/512.0)*sqrt(4849845)*(5*cos(2*theta) + 4)*cexp(-9*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dtheta_{10}^{-10} */
static double complex dYlm_dtheta_l10m_10(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(969969)*cexp(-10*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-1} */
static double complex dYlm_dtheta_l11m_1(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(759)*(-278460*pow(sin(theta), 6) + 488670*pow(sin(theta), 4) - 200655*pow(sin(theta), 2) + 323323*pow(cos(theta), 10) - 860795*pow(cos(theta), 8) + 550290*pow(cos(theta), 6) - 12306)*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-2} */
static double complex dYlm_dtheta_l11m_2(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(98670)*(24871*pow(cos(theta), 10) - 55233*pow(cos(theta), 8) + 42126*pow(cos(theta), 6) - 12810*pow(cos(theta), 4) + 1323*pow(cos(theta), 2) - 21)*cexp(-2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-3} */
static double complex dYlm_dtheta_l11m_3(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(345345)*(2584*pow(sin(theta), 6) - 3978*pow(sin(theta), 4) + 984*pow(sin(theta), 2) + 3553*pow(cos(theta), 8) - 3876*pow(cos(theta), 6) + 451)*cexp(-3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-4} */
static double complex dYlm_dtheta_l11m_4(const double theta, const double phi)
{
  return (3.0/512.0)*sqrt(46046)*(1615*pow(sin(theta), 6) - 2635*pow(sin(theta), 4) + 705*pow(sin(theta), 2) + 3553*pow(cos(theta), 8) - 3553*pow(cos(theta), 6) + 320)*cexp(-4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-5} */
static double complex dYlm_dtheta_l11m_5(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(6578)*(-13566*pow(sin(theta), 6) + 12597*pow(sin(theta), 4) + 7259*pow(sin(theta), 2) + 11305*pow(cos(theta), 6) - 6825)*cexp(-5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-6} */
static double complex dYlm_dtheta_l11m_6(const double theta, const double phi)
{
  return (3.0/1024.0)*sqrt(167739)*(-665*pow(sin(theta), 6) + 760*pow(sin(theta), 4) + 250*pow(sin(theta), 2) + 798*pow(cos(theta), 6) - 350)*cexp(-6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-7} */
static double complex dYlm_dtheta_l11m_7(const double theta, const double phi)
{
  return (1.0/2048.0)*sqrt(1677390)*(1463*pow(sin(theta), 4) - 2052*pow(sin(theta), 2) + 672)*cexp(-7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-8} */
static double complex dYlm_dtheta_l11m_8(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(31870410)*(77*pow(sin(theta), 4) - 124*pow(sin(theta), 2) + 48)*cexp(-8*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-9} */
static double complex dYlm_dtheta_l11m_9(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2124694)*(60 - 77*pow(sin(theta), 2))*cexp(-9*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-10} */
static double complex dYlm_dtheta_l11m_10(const double theta, const double phi)
{
  return (1.0/1024.0)*sqrt(22309287)*(10 - 11*pow(sin(theta), 2))*cexp(-10*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dtheta_{11}^{-11} */
static double complex dYlm_dtheta_l11m_11(const double theta, const double phi)
{
  return (11.0/2048.0)*sqrt(4056234)*cexp(-11*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-1} */
static double complex dYlm_dtheta_l12m_1(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(39)*(196350*pow(sin(theta), 6) - 369600*pow(sin(theta), 4) + 167937*pow(sin(theta), 2) + 624036*pow(cos(theta), 12) - 1815583*pow(cos(theta), 10) + 1971915*pow(cos(theta), 8) - 785400*pow(cos(theta), 6) + 5544)*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-2} */
static double complex dYlm_dtheta_l12m_2(const double theta, const double phi)
{
  return (5.0/512.0)*sqrt(6006)*(22287*pow(cos(theta), 10) - 54910*pow(cos(theta), 8) + 48450*pow(cos(theta), 6) - 18360*pow(cos(theta), 4) + 2775*pow(cos(theta), 2) - 114)*cexp(-2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-3} */
static double complex dYlm_dtheta_l12m_3(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*(-9690*pow(sin(theta), 6) + 17340*pow(sin(theta), 4) - 6690*pow(sin(theta), 2) + 29716*pow(cos(theta), 10) - 61047*pow(cos(theta), 8) + 32946*pow(cos(theta), 6) - 975)*cexp(-3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-4} */
static double complex dYlm_dtheta_l12m_4(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(1001)*(13566*pow(sin(theta), 6) - 20672*pow(sin(theta), 4) + 4386*pow(sin(theta), 2) + 22287*pow(cos(theta), 8) - 23902*pow(cos(theta), 6) + 2895)*cexp(-4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-5} */
static double complex dYlm_dtheta_l12m_5(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(34034)*(1995*pow(sin(theta), 6) - 3230*pow(sin(theta), 4) + 790*pow(sin(theta), 2) + 5244*pow(cos(theta), 8) - 5054*pow(cos(theta), 6) + 450)*cexp(-5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-6} */
static double complex dYlm_dtheta_l12m_6(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(2431)*(-3059*pow(sin(theta), 6) + 2793*pow(sin(theta), 4) + 1881*pow(sin(theta), 2) + 3059*pow(cos(theta), 6) - 1715)*cexp(-6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-7} */
static double complex dYlm_dtheta_l12m_7(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(277134)*(-805*pow(sin(theta), 6) + 910*pow(sin(theta), 4) + 345*pow(sin(theta), 2) + 1127*pow(cos(theta), 6) - 455)*cexp(-7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-8} */
static double complex dYlm_dtheta_l12m_8(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(277134)*(483*pow(sin(theta), 4) - 700*pow(sin(theta), 2) + 240)*cexp(-8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-9} */
static double complex dYlm_dtheta_l12m_9(const double theta, const double phi)
{
  return (15.0/2048.0)*sqrt(646646)*(92*pow(sin(theta), 4) - 151*pow(sin(theta), 2) + 60)*cexp(-9*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-10} */
static double complex dYlm_dtheta_l12m_10(const double theta, const double phi)
{
  return (5.0/1024.0)*sqrt(88179)*(110 - 138*pow(sin(theta), 2))*cexp(-10*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-11} */
static double complex dYlm_dtheta_l12m_11(const double theta, const double phi)
{
  return (5.0/2048.0)*sqrt(4056234)*(11 - 12*pow(sin(theta), 2))*cexp(-11*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dtheta_{12}^{-12} */
static double complex dYlm_dtheta_l12m_12(const double theta, const double phi)
{
  return (15.0/1024.0)*sqrt(676039)*cexp(-12*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-1} */
static double complex dYlm_dtheta_l13m_1(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(546)*(1279080*pow(sin(theta), 6) - 2347785*pow(sin(theta), 4) + 1035540*pow(sin(theta), 2) + 2414425*pow(cos(theta), 12) - 7622154*pow(cos(theta), 10) + 9220035*pow(cos(theta), 8) - 4050420*pow(cos(theta), 6) + 39138)*cexp(-I*phi)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-2} */
static double complex dYlm_dtheta_l13m_2(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2730)*(482885*pow(cos(theta), 12) - 1307504*pow(cos(theta), 10) + 1311057*pow(cos(theta), 8) - 596904*pow(cos(theta), 6) + 120615*pow(cos(theta), 4) - 8712*pow(cos(theta), 2) + 99)*cexp(-2*I*phi)*sin(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-3} */
static double complex dYlm_dtheta_l13m_3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(30030)*(-244188*pow(sin(theta), 6) + 420546*pow(sin(theta), 4) - 151113*pow(sin(theta), 2) + 482885*pow(cos(theta), 10) - 1106921*pow(cos(theta), 8) + 656982*pow(cos(theta), 6) - 26802)*cexp(-3*I*phi)*pow(sin(theta), 2)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-4} */
static double complex dYlm_dtheta_l13m_4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(51051)*(-35910*pow(sin(theta), 6) + 63840*pow(sin(theta), 4) - 23595*pow(sin(theta), 2) + 142025*pow(cos(theta), 10) - 271377*pow(cos(theta), 8) + 138852*pow(cos(theta), 6) - 4380)*cexp(-4*I*phi)*pow(sin(theta), 3)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-5} */
static double complex dYlm_dtheta_l13m_5(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(102102)*(73416*pow(sin(theta), 6) - 110922*pow(sin(theta), 4) + 20216*pow(sin(theta), 2) + 142025*pow(cos(theta), 8) - 148580*pow(cos(theta), 6) + 18075)*cexp(-5*I*phi)*pow(sin(theta), 4)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-6} */
static double complex dYlm_dtheta_l13m_6(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(969969)*(2415*pow(sin(theta), 6) - 3885*pow(sin(theta), 4) + 875*pow(sin(theta), 2) + 7475*pow(cos(theta), 8) - 6923*pow(cos(theta), 6) + 600)*cexp(-6*I*phi)*pow(sin(theta), 5)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-7} */
static double complex dYlm_dtheta_l13m_7(const double theta, const double phi)
{
  return (21.0/4096.0)*sqrt(692835)*(-690*pow(sin(theta), 6) + 621*pow(sin(theta), 4) + 471*pow(sin(theta), 2) + 805*pow(cos(theta), 6) - 421)*cexp(-7*I*phi)*pow(sin(theta), 6)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-8} */
static double complex dYlm_dtheta_l13m_8(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(9699690)*(-575*pow(sin(theta), 6) + 644*pow(sin(theta), 4) + 272*pow(sin(theta), 2) + 920*pow(cos(theta), 6) - 344)*cexp(-8*I*phi)*pow(sin(theta), 7)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-9} */
static double complex dYlm_dtheta_l13m_9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(88179)*(7475*pow(sin(theta), 4) - 11132*pow(sin(theta), 2) + 3960)*cexp(-9*I*phi)*pow(sin(theta), 8)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-10} */
static double complex dYlm_dtheta_l13m_10(const double theta, const double phi)
{
  return (3.0/2048.0)*sqrt(2028117)*(325*pow(sin(theta), 4) - 542*pow(sin(theta), 2) + 220)*cexp(-10*I*phi)*pow(sin(theta), 9)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-11} */
static double complex dYlm_dtheta_l13m_11(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(4056234)*(264 - 325*pow(sin(theta), 2))*cexp(-11*I*phi)*pow(sin(theta), 10)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-12} */
static double complex dYlm_dtheta_l13m_12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(2028117)*(12 - 13*pow(sin(theta), 2))*cexp(-12*I*phi)*pow(sin(theta), 11)/sqrt(M_PI);
}


/* dY_dtheta_{13}^{-13} */
static double complex dYlm_dtheta_l13m_13(const double theta, const double phi)
{
  return (195.0/8192.0)*sqrt(312018)*cexp(-13*I*phi)*pow(sin(theta), 12)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-1} */
static double complex dYlm_dtheta_l14m_1(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(6090)*(-692835*pow(sin(theta), 6) + 1327326*pow(sin(theta), 4) - 620763*pow(sin(theta), 2) + 4680270*pow(cos(theta), 14) - 15935205*pow(cos(theta), 12) + 21246940*pow(cos(theta), 10) - 13995267*pow(cos(theta), 8) + 4018443*pow(cos(theta), 6) - 14157)*cexp(-I*phi)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-2} */
static double complex dYlm_dtheta_l14m_2(const double theta, const double phi)
{
  return (1.0/8192.0)*sqrt(79170)*(2340135*pow(cos(theta), 12) - 6908970*pow(cos(theta), 10) + 7763305*pow(cos(theta), 8) - 4135692*pow(cos(theta), 6) + 1055241*pow(cos(theta), 4) - 113322*pow(cos(theta), 2) + 3399)*cexp(-2*I*phi)*sin(theta)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-3} */
static double complex dYlm_dtheta_l14m_3(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(448630)*(43890*pow(sin(theta), 6) - 81510*pow(sin(theta), 4) + 34617*pow(sin(theta), 2) + 275310*pow(cos(theta), 12) - 697015*pow(cos(theta), 10) + 648945*pow(cos(theta), 8) - 228228*pow(cos(theta), 6) + 3036)*cexp(-3*I*phi)*pow(sin(theta), 2)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-4} */
static double complex dYlm_dtheta_l14m_4(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(2467465)*(-18354*pow(sin(theta), 6) + 31388*pow(sin(theta), 4) - 10659*pow(sin(theta), 2) + 45885*pow(cos(theta), 10) - 98325*pow(cos(theta), 8) + 55936*pow(cos(theta), 6) - 2472)*cexp(-4*I*phi)*pow(sin(theta), 3)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-5} */
static double complex dYlm_dtheta_l14m_5(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(18752734)*(-966*pow(sin(theta), 6) + 1708*pow(sin(theta), 4) - 608*pow(sin(theta), 2) + 4830*pow(cos(theta), 10) - 8625*pow(cos(theta), 8) + 4186*pow(cos(theta), 6) - 135)*cexp(-5*I*phi)*pow(sin(theta), 4)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-6} */
static double complex dYlm_dtheta_l14m_6(const double theta, const double phi)
{
  return (3.0/8192.0)*sqrt(93763670)*(3220*pow(sin(theta), 6) - 4830*pow(sin(theta), 4) + 756*pow(sin(theta), 2) + 7245*pow(cos(theta), 8) - 7360*pow(cos(theta), 6) + 883)*cexp(-6*I*phi)*pow(sin(theta), 5)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-7} */
static double complex dYlm_dtheta_l14m_7(const double theta, const double phi)
{
  return (7.0/4096.0)*sqrt(20092215)*(575*pow(sin(theta), 6) - 920*pow(sin(theta), 4) + 192*pow(sin(theta), 2) + 2070*pow(cos(theta), 8) - 1840*pow(cos(theta), 6) + 154)*cexp(-7*I*phi)*pow(sin(theta), 6)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-8} */
static double complex dYlm_dtheta_l14m_8(const double theta, const double phi)
{
  return (1.0/4096.0)*sqrt(25571910)*(-3105*pow(sin(theta), 6) + 2760*pow(sin(theta), 4) + 2300*pow(sin(theta), 2) + 4140*pow(cos(theta), 6) - 2028)*cexp(-8*I*phi)*pow(sin(theta), 7)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-9} */
static double complex dYlm_dtheta_l14m_9(const double theta, const double phi)
{
  return (3.0/4096.0)*sqrt(98025655)*(-225*pow(sin(theta), 6) + 250*pow(sin(theta), 4) + 115*pow(sin(theta), 2) + 405*pow(cos(theta), 6) - 141)*cexp(-9*I*phi)*pow(sin(theta), 8)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-10} */
static double complex dYlm_dtheta_l14m_10(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(315*pow(sin(theta), 4) - 480*pow(sin(theta), 2) + 176)*cexp(-10*I*phi)*pow(sin(theta), 9)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-11} */
static double complex dYlm_dtheta_l14m_11(const double theta, const double phi)
{
  return (5.0/8192.0)*sqrt(117630786)*(126*pow(sin(theta), 4) - 213*pow(sin(theta), 2) + 88)*cexp(-11*I*phi)*pow(sin(theta), 10)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-12} */
static double complex dYlm_dtheta_l14m_12(const double theta, const double phi)
{
  return (15.0/4096.0)*sqrt(1508087)*(52 - 63*pow(sin(theta), 2))*cexp(-12*I*phi)*pow(sin(theta), 11)*cos(theta)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-13} */
static double complex dYlm_dtheta_l14m_13(const double theta, const double phi)
{
  return (15.0/8192.0)*sqrt(9048522)*(13 - 14*pow(sin(theta), 2))*cexp(-13*I*phi)*pow(sin(theta), 12)/sqrt(M_PI);
}


/* dY_dtheta_{14}^{-14} */
static double complex dYlm_dtheta_l14m_14(const double theta, const double phi)
{
  return (105.0/8192.0)*sqrt(1292646)*cexp(-14*I*phi)*pow(sin(theta), 13)*cos(theta)/sqrt(M_PI);
}

