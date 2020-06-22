#include "bbn_headers.h"

void bbn_populate_free_data(Grid_T *const grid);
void bbn_free_data_gammas(Grid_T *const grid);
void bbn_free_conformal_metric_derivatives(Patch_T *const patch);
void bbn_preparing_conformal_metric_derivatives(Patch_T *const patch);
void bbn_free_data_Gamma(Grid_T *const grid);
void bbn_free_data_Ricci(Grid_T *const grid);
void bbn_free_data_KS_trKij(Patch_T *const patch);
static double KerrSchild_H(const double M_BH,const double _r,const double a,const double _z);
static double KerrSchild_r(const double _x,const double _y,const double _z,const double _a);
void bbn_free_data_dGamma(Grid_T *const grid);
void bbn_free_data_tr_KSKij(Grid_T *const grid);
static void populating_KSGamma(Patch_T *const patch);
static void populate_KSgammas_KSalpha_KSBeta(Patch_T *const patch);
static void partial_derivative_KSBeta(Patch_T *const patch);
static void free_KSfields(Patch_T *const patch);
static void execute_boost_and_rotation(Transformation_T *const tB,Transformation_T *const tR,const int IsInverse,const double *const in,double *const out);
void bbn_transform_populate_boost_rotation(Transformation_T *const tB,Transformation_T *const tR);
void bbn_transform_get_k_and_H_KerrSchild(const double x,const double y,const double z,const double a,const double m,Transformation_T *const tB,Transformation_T *const tR,double *const kt,double *const k0,double *const k1,double *const k2,double *const H);



