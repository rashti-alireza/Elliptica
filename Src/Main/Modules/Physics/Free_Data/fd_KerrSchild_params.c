/*
// Alireza Rashti
// September 2020
*/

/* setting parameters e.g. global variables for 
// analytic Kerr-Schild free data.
// NOTE: before each call and for each BH one must updata the paramters.
// NOTE: NOT thread safe. */

#include "fd_header.h"

/* for external variables. NOTE: MUST be the same as 
// the one in 'fd_KerrSchild_header.h' DON'T CHANGE. */
#define KS_glob_var(x) fd_ks_glob##x

#define M_BH    KS_glob_var(M_BH) /* BH mass */
#define a_BH    KS_glob_var(a_BH) /* BH spin */
#define phiy    KS_glob_var(phiy) /* rotation angel */
#define phiz    KS_glob_var(phiz) /* rotation angel */
#define Bx      KS_glob_var(Bx) /* B=v/c (boost) */
#define By      KS_glob_var(By) /* B=v/c (boost) */
#define Bz      KS_glob_var(Bz) /* B=v/c (boost) */
#define B2      KS_glob_var(B2) /* B^i B_i */
#define Lambda  KS_glob_var(Lambda) /* flat data => 0, kerr-schild => 1 */

/* global variables */
double M_BH,a_BH;/* mass and spin of BH */
double phiy,phiz;/* rotation angels */
double Bx,By,Bz;/* B=v/c (boost) */
double B2;/* B^i B_i */
double Lambda;/* flat data => 0, kerr-schild => 1 */


void fd_KerrSchild_set_params(Physics_T *const phys);
void fd_KerrSchild_set_params(Physics_T *const phys)
{
  const double Epsilon = 1E-14;
  const double chi_x   = Getd("chi_x");
  const double chi_y   = Getd("chi_y");
  const double chi_z   = Getd("chi_z");
  const double chi     = sqrt(Pow2(chi_x)+Pow2(chi_y)+Pow2(chi_z));
  
  M_BH = Getd("Christodoulou_mass");
  a_BH = Getd("spin_a");
  Lambda = 1.;

  assert(LSSEQL(chi,1));

  /* boost */
  Bx = Getd("boost_Vx");
  By = Getd("boost_Vy");
  Bz = Getd("boost_Vz");
  B2 = Pow2(Bx)+Pow2(By)+Pow2(Bz);

 /* avoid division by zero */
  if(EQL(B2,0.)) 
  {
    Bx = 0.0;
    By = 0.0;
    Bz = 0.0;
    B2 = Epsilon;
  }
  
  /* rotation */
  if (!EQL(chi,0.))/* otherwise R is 0 */
  {
    phiz = arctan(chi_y,chi_x);
    phiy = acos(chi_z/chi);
    assert(isfinite(phiy));
  }
  else
  {
    phiz = 0.;
    phiy = 0.;
  }

}

