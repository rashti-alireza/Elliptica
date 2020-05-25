#include "bbn_headers.h"

void bbn_populate_free_data(Grid_T *const grid);
void bbn_free_data_gammas(Grid_T *const grid);
void bbn_free_conformal_metric_derivatives(Patch_T *const patch);
void bbn_preparing_conformal_metric_derivatives(Patch_T *const patch);
void bbn_free_data_Gamma(Grid_T *const grid);
void bbn_free_data_Ricci(Grid_T *const grid);
void bbn_free_data_KS_trKij(Patch_T *const patch);
double bbn_KerrSchild_H(const double M_BH,const double rbar,const double a,const double z);
double bbn_KerrShcild_r(const double x,const double y,const double z,const double a);
void bbn_free_data_dGamma(Grid_T *const grid);
void bbn_free_data_tr_KSKij(Grid_T *const grid);
static void populating_KSGamma(Patch_T *const patch);
static void populate_KSgammas_KSalpha_KSBeta(Patch_T *const patch);
static void partial_derivative_KSBeta(Patch_T *const patch);
static void free_KSfields(Patch_T *const patch);
void bbn_boost_rotation_transform(Transformation_T *const tB,Transformation_T *const tR,const int IsInverse,const double in[4],double out[4]);
void bbn_populate_tB_tR(Transformation_T *const tB,Transformation_T *const tR);






