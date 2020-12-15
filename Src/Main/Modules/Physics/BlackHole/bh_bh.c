/*
// Alireza Rashti
// November 2020
*/


/* collection of black hole affairs */

#include "bh_bh.h"


/* use KerrSchild assuming pefect S2 to start off black hole
// paramters and domain shape etc.
// NOTE: this is perfect sphere, so chi = 0. and spin = 0. */
void bh_start_off_KerrSchild_perfect_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  const double bh_chi_x    = 0.;
  const double bh_chi_y    = 0.;
  const double bh_chi_z    = 0.;
  const double bh_irr_mass = Getd("irreducible_mass");
  const double bh_R        = 2.0*bh_irr_mass;/* approximate initial radius */
  const double bh_a        = 0.;
  
  /* set initial params */
  Setd("perfect_S2_radius",bh_R);
  Setd("min_radius",bh_R);
  Setd("max_radius",bh_R);
  Setd("Christodoulou_mass",bh_irr_mass);
  Setd("spin_a",bh_a);
  
  printf(Pretty0"%s properties:\n",phys->stype);
  printf(Pretty0"%s radius (Kerr-Schild Coords.) = %+e\n",phys->stype,bh_R);
  printf(Pretty0"%s irreducible mass             = %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,bh_chi_x);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,bh_chi_y);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,bh_chi_z);
  printf(Pretty0"%s spin/M(= a)                  = %+e\n",phys->stype,bh_a);
  
  FUNC_TOC
}

/* use a general KerrSchild to start off black hole
// print properties.
// NOTE: this is NOT assuming perfect sphere. */
void bh_start_off_KerrSchild_general_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  const double bh_chi_x    = Getd("chi_x");
  const double bh_chi_y    = Getd("chi_y");
  const double bh_chi_z    = Getd("chi_z");
  const double bh_irr_mass = Getd("irreducible_mass");
  const double bh_chi = sqrt(Pow2(bh_chi_x)+Pow2(bh_chi_y)+Pow2(bh_chi_z));
  const double bh_chr_mass = bh_irr_mass*sqrt(2./(1.+sqrt(1-Pow2(bh_chi))));
  const double bh_R = bh_chr_mass*(1+sqrt(1-Pow2(bh_chi)));/* approximate initial radius */
  const double bh_a = bh_chi*bh_chr_mass;
  
  /* check size of bh_chi */
  if (GRT(bh_chi,1.))
    Errors("%s spin is too large!\n",phys->stype);
  
  /* set initial params */
  Setd("Christodoulou_mass",bh_irr_mass);
  Setd("spin_a",bh_a);
  
  printf(Pretty0"%s properties:\n",phys->stype);
  printf(Pretty0"%s irreducible mass             = %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s Christodoulou_mass           = %+e\n",phys->stype,bh_chr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,bh_chi_x);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,bh_chi_y);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,bh_chi_z);
  printf(Pretty0"%s net dimensionless spin       = %+e\n",phys->stype,bh_chi);
  printf(Pretty0"%s Spin/M (= a)                 = %+e\n",phys->stype,bh_a);
  printf(Pretty0"%s approximate radius           ~ %+e\n",phys->stype,bh_R);
  
  FUNC_TOC
}

/* use Schwarzchild in isotropic coords (pefect S2) to start off 
// black hole paramters and domain shape etc.
// NOTE: chi = 0. and spin = 0. */
void bh_start_off_IsoSchild_perfect_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  const double bh_irr_mass = Getd("irreducible_mass");
  const double bh_R        = 0.5*bh_irr_mass;/* approximate initial radius */
  
  /* set initial params */
  Setd("perfect_S2_radius",bh_R);
  Setd("min_radius",bh_R);
  Setd("max_radius",bh_R);
  Setd("Christodoulou_mass",bh_irr_mass);
  Setd("spin_a",0.);

  printf(Pretty0"%s properties:\n",phys->stype);
  printf(Pretty0"%s radius (isotropic Coords.)   = %+e\n",phys->stype,bh_R);
  printf(Pretty0"%s irreducible mass             = %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,0.);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,0.);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,0.);
  
  FUNC_TOC
}

/* use Schwarzchild in Painleve-Gullstrand coords (pefect S2) to start off 
// black hole paramters and domain shape etc.
// NOTE: this is perfect sphere, and chi = 0. and spin = 0. */
void bh_start_off_PGSchild_perfect_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  const double bh_irr_mass = Getd("irreducible_mass");
  const double bh_R        = 2*bh_irr_mass;/* approximate initial radius */
  
  /* set initial params */
  Setd("perfect_S2_radius",bh_R);
  Setd("min_radius",bh_R);
  Setd("max_radius",bh_R);
  Setd("Christodoulou_mass",bh_irr_mass);
  Setd("spin_a",0.);

  printf(Pretty0"%s properties:\n",phys->stype);
  printf(Pretty0"%s radius (Painleve-Gullstrand) = %+e\n",phys->stype,bh_R);
  printf(Pretty0"%s irreducible mass             = %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,0.);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,0.);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,0.);
  
  FUNC_TOC
}


/* find BH surface and then set grid characteristic */
void bh_find_bh_surface_perfect_s2(Physics_T *const phys)
{
  Grid_Char_T *grid_char = phys->grid_char;
  const Uint lmax   = (Uint)Geti("surface_Ylm_max_l");
  const Uint Ntheta = Ntheta_Ylm(lmax);
  const Uint Nphi   = Nphi_Ylm(lmax);
  const Uint Ntot   = Ntotal_Ylm(lmax);
  const double R_BH     = Getd("perfect_S2_radius");
  double *rbh = alloc_double(Ntot);/* surface function r = r(th,ph). */
  double *reClm_rbh = alloc_ClmYlm(lmax),
         *imClm_rbh = alloc_ClmYlm(lmax);
  Uint ij;
  
  init_Legendre_root_function();
  for (ij = 0; ij < Ntot; ++ij)
  {
    rbh[ij] = R_BH;
  }
  /* calculating coeffs */
  get_Ylm_coeffs(reClm_rbh,imClm_rbh,rbh,Ntheta,Nphi,lmax);
  
  assert(!grid_char->params[phys->igc]->occupied);
  grid_char->params[phys->igc]->obj    = phys->stype;
  grid_char->params[phys->igc]->dir    = phys->spos;
  grid_char->params[phys->igc]->relClm = reClm_rbh;
  grid_char->params[phys->igc]->imgClm = imClm_rbh;
  grid_char->params[phys->igc]->r_min  = Getd("min_radius");
  grid_char->params[phys->igc]->r_max  = Getd("max_radius");
  grid_char->params[phys->igc]->lmax   = lmax;
  grid_char->params[phys->igc]->occupied = 1;
  
  Free(rbh);
}

/* set BH surface for a general Kerr-Schild BH with arbitrary 
// spin and boost and then set grid characteristic */
void bh_find_bh_surface_KerrSchild_s2(Physics_T *const phys)
{
  Grid_Char_T *grid_char = phys->grid_char;
  const Uint lmax   = (Uint)Geti("surface_Ylm_max_l");
  const Uint Ntheta = Ntheta_Ylm(lmax);
  const Uint Nphi   = Nphi_Ylm(lmax);
  const Uint Ntot   = Ntotal_Ylm(lmax);
  double *rbh       = fd_find_KerrSchild_BH_surface(phys);
  double rmin = DBL_MAX,rmax = 0.;
  double *reClm_rbh = alloc_ClmYlm(lmax),
         *imClm_rbh = alloc_ClmYlm(lmax);

  for (Uint ij = 0; ij < Ntot; ++ij)
  {
    rmin = (rbh[ij] < rmin ? rbh[ij]: rmin);
    rmax = (rbh[ij] > rmax ? rbh[ij]: rmax);
  }
  
  Setd("min_radius",rmin);
  Setd("max_radius",rmax);
  
  /* calculating coeffs */
  get_Ylm_coeffs(reClm_rbh,imClm_rbh,rbh,Ntheta,Nphi,lmax);
  
  assert(!grid_char->params[phys->igc]->occupied);
  grid_char->params[phys->igc]->obj    = phys->stype;
  grid_char->params[phys->igc]->dir    = phys->spos;
  grid_char->params[phys->igc]->relClm = reClm_rbh;
  grid_char->params[phys->igc]->imgClm = imClm_rbh;
  grid_char->params[phys->igc]->r_min  = rmin;
  grid_char->params[phys->igc]->r_max  = rmax;
  grid_char->params[phys->igc]->lmax   = lmax;
  grid_char->params[phys->igc]->occupied = 1;
  
  Free(rbh);
}


/* adjust the BH radius to acquire the desired BH irreducible mass
// when the bh is a perfect sphere in x coords. */
void bh_tune_BH_radius_irreducible_mass_perfect_s2(Physics_T *const phys)
{
  const double target_bh_mass = Getd("irreducible_mass");
  const double current_r_bh   = Getd("perfect_S2_radius");
  const double W              = Getd("radius_update_weight");
  const double dM_tolerance   = Getd("mass_tolerance");
  double kommar_mass, adm_mass,irr_mass;
  double dr,r_bh,dM;
  
  observe(phys,"Irreducible(M)",&irr_mass);
  observe(phys,"ADM(M)",&adm_mass);
  observe(phys,"Kommar(M)",&kommar_mass);
  
  printf(Pretty0"current BH irreducible mass = %e\n",irr_mass);
  printf(Pretty0"current BH ADM mass         = %e\n",adm_mass);
  printf(Pretty0"current BH Kommar mass      = %e\n",kommar_mass);
  
  Setd("irreducible_mass_current",irr_mass);
  Setd("ADM_mass",adm_mass);
  Setd("Kommar_mass",kommar_mass);
  
  dM = fabs(irr_mass/target_bh_mass-1);
  dr = -current_r_bh*(irr_mass/target_bh_mass-1);
  if (EQL(W,0))
  {
    dr = 0;
    printf(Pretty0"updating weight factor is zero.\n");
  }
  if (LSSEQL(dM,dM_tolerance)) 
  {
    dr = 0;
    printf(Pretty0"|dM/M| = %g < Tol. = %g\n",dM,dM_tolerance);
  }

  r_bh = current_r_bh + W*dr;
  
  Setd("min_radius",r_bh);
  Setd("max_radius",r_bh);
  
  if (EQL(dr,0))/* => no change in AH surface */
  {
    Seti("did_AH_surface_change?",0);
  }
  else          /* => change in AH surface */
  {
    Seti("did_AH_surface_change?",1);
  }
  
}

