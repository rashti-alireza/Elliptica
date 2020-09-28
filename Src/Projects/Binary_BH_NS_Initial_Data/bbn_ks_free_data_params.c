/*
// Alireza Rashti
// September 2020
*/

/* setting parameters e.g. global variables for 
// analytic Kerr-Schild free data  NOTE: NOT thread safe. */

#include "bbn_headers.h"

/* for external variables. NOTE: MUST be the same as 
// the one in 'bbn_ks_free_data_analytic.h' DON'T CHANGE. */
#define KS_glob_var(x) bbn_ks_glob##x

#define M_BH    KS_glob_var(M_BH) /* BH mass */
#define a_BH    KS_glob_var(a_BH) /* BH spin */
#define phiy    KS_glob_var(phiy) /* rotation angel */
#define phiz    KS_glob_var(phiz) /* rotation angel */
#define Bx      KS_glob_var(Bx) /* B=v/c (boost) */
#define By      KS_glob_var(By) /* B=v/c (boost) */
#define Bz      KS_glob_var(Bz) /* B=v/c (boost) */
#define B2      KS_glob_var(B2) /* B^i B_i */
#define r0      KS_glob_var(r0) /* roll off radius */
#define Lambda  KS_glob_var(Lambda) /* flat data => 0, kerr-schild => 1 */

/* global variables */
double M_BH,a_BH;/* mass and spin of BH */
double phiy,phiz;/* rotation angels */
double Bx,By,Bz;/* B=v/c (boost) */
double B2;/* B^i B_i */
double r0;/* roll off radius */
double Lambda;/* flat data => 0, kerr-schild => 1 */


void bbn_ks_free_data_set_params(Grid_T *const grid);
void bbn_ks_free_data_set_params(Grid_T *const grid)
{
  const double BH_center_x = Pgetd("BH_center_x");
  const double BH_center_y = Pgetd("BH_center_y");
  const double chi_U0      = Pgetd("BH_chi_U0");
  const double chi_U1      = Pgetd("BH_chi_U1");
  const double chi_U2      = Pgetd("BH_chi_U2");
  const double y_CM        = Pgetd("y_CM");
  const double x_CM        = Pgetd("x_CM");
  const double Omega_BHNS  = Pgetd("BH_NS_angular_velocity");
  const double chi         = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  
  r0   = Pgetd("BH_KerrSchild_RollOff");
  M_BH = Pgetd("BH_irreducible_mass");
  a_BH = Pgetd("BH_net_spin");
  
  if (Pcmps("BH_NS_free_data_metric","conformally_flat_metric"))
    Lambda = 0;
  else if (Pcmps("BH_NS_free_data_metric","Boosted_KerrSchild_metric"))
    Lambda = 1;
  else
    Error0(NO_OPTION);

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = -Omega_BHNS*(BH_center_y-y_CM);
  By =  Omega_BHNS*(BH_center_x-x_CM);
  Bz = Pgetd("BH_Vz");
  B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);

  assert(!EQL(B2,0));
  /* rotation */
  if (!EQL(chi,0))/* otherwise R is 0 */
  {
    phiz = arctan(chi_U1,chi_U0);
    phiy = acos(chi_U2/chi);
    assert(isfinite(phiy));
  }
  else
  {
    phiz = 0.;
    phiy = 0.;
  }

  UNUSED(grid);
}

/* ->: pointer to field to write.
// if does not exist add field. 
// note: this version is simple and params is just the patch. */
double *bbn_ks_read_analytic(const char *const name, void *params);
double *bbn_ks_read_analytic(const char *const name, void *params)
{
  Patch_T *const patch = params;
  Field_T *u   = 0;
  const int fn = _Ind(name);
  
  if (fn < 0)/* if there is not such a field */
  {
    u = add_field(name,0,patch,YES);
  }
  else
  {
    u = patch->pool[fn];
    empty_field(u);
    u->v = alloc_double(patch->nn);
  }
  
  return u->v;
}



