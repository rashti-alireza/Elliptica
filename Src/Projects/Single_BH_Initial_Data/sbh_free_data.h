#include "sbh_headers.h"

void sbh_populate_free_data(Grid_T *const grid);
void sbh_free_data_gammas(Grid_T *const grid);
void sbh_free_conformal_metric_derivatives(Patch_T *const patch);
void sbh_preparing_conformal_metric_derivatives(Patch_T *const patch);
void sbh_free_data_Gamma(Grid_T *const grid);
void sbh_free_data_Ricci(Grid_T *const grid);
void sbh_free_data_KS_trKij(Patch_T *const patch);
static void sbh_free_data_dGamma(Grid_T *const grid);
static void populating_KSGamma(Patch_T *const patch);
static void populate_KSgammas_KSalpha_KSBeta(Patch_T *const patch);
static void partial_derivative_KSBeta(Patch_T *const patch);
static void free_KSfields(Patch_T *const patch);
static void sbh_free_data_tr_KSKij(Grid_T *const grid);





