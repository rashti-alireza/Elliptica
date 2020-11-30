/*
// Alireza Rashti
// September 2020
*/

/* setting parameters e.g. global variables for 
// analytic Kerr-Schild free data.
// NOTE: before each call and for each BH one must updata the paramters.
// NOTE: NOT thread safe. */

#include "frda_header.h"

/* for external variables. NOTE: MUST be the same as 
// the one in 'frda_ks_free_data_header.h' DON'T CHANGE. */
#define KS_glob_var(x) frda_ks_glob##x

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


void frda_ks_free_data_set_params(Physics_T *const phys);
void frda_ks_free_data_set_params(Physics_T *const phys)
{
  const double BH_center_x = Getd("center_x");
  const double BH_center_y = Getd("center_y");
  const double chi_U0      = Getd("chi_U0");
  const double chi_U1      = Getd("chi_U1");
  const double chi_U2      = Getd("chi_U2");
  const double y_CM        = Getd("y_CM");
  const double x_CM        = Getd("x_CM");
  const double Omega       = sysGetd("angular_velocity");
  const double chi         = sqrt(Pow2(chi_U0)+Pow2(chi_U1)+Pow2(chi_U2));
  
  r0   = Getd("KerrSchild_RollOff");
  M_BH = Getd("irreducible_mass");
  a_BH = Getd("net_spin");
  Lambda = 1;

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = -Omega*(BH_center_y-y_CM);
  By =  Omega*(BH_center_x-x_CM);
  Bz = Getd("Vz");
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

}

