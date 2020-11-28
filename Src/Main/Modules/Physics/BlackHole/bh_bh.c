/*
// Alireza Rashti
// November 2020
*/


/* collection of black hole affairs */

#include "bh_bh.h"


/* adjust AH radius to meet a criteria for instant the mass is fixed */
int bh_tune_black_hole_radius(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("tune_BH_radius_criteria","fix_irreducible_mass")
  {
    IF_sval("surface_type","perfect_s2")
      tune_BH_radius_irreducible_mass_perfect_s2(phys);
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* find bh surface, mainly for grid setup */
int bh_find_black_hole_surface(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("surface_type","perfect_s2")
  {
    find_bh_surface_perfect_s2(phys);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* set initial feature of blach holes, mostly used 
// for the very first time that we need to make grid. */
int bh_start_off(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("start_off","KerrSchild")
  {
    IF_sval("surface_type","perfect_s2")
      start_off_KerrSchild_perfect_s2(phys);
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* use KerrSchild assuming pefect S2 to start off black hole
// paramters and domain shape etc. */
static void start_off_KerrSchild_perfect_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  const double bh_chi_x    = Getd("chi_U0");
  const double bh_chi_y    = Getd("chi_U1");
  const double bh_chi_z    = Getd("chi_U2");
  const double bh_irr_mass = Getd("irreducible_mass");
  const double bh_R        = 1.5*bh_irr_mass;/* approximate initial radius */
  const double bh_chi      = sqrt(Pow2(bh_chi_x)+Pow2(bh_chi_y)+Pow2(bh_chi_z));
  const double bh_a        = bh_chi*bh_irr_mass;
  
  /* check size of bh_chi */
  if (GRT(bh_chi,1))
    Error0("BH spin is too large!\n");
  
  /* set initial grid parameters */
  Setd("perfect_S2_radius",bh_R);
  Setd("min_radius",bh_R);
  Setd("max_radius",bh_R);
  
  printf("%s properties:\n",phys->stype);
  printf(Pretty0"%s radius (Kerr-Schild Coords.) ~ %+e\n",phys->stype,bh_R);
  printf(Pretty0"%s irreducible mass             ~ %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,bh_chi_x);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,bh_chi_y);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,bh_chi_z);
  printf(Pretty0"%s approximate net spin         ~ %+e\n",phys->stype,bh_a);
  
  FUNC_TOC
}


/* find BH surface and then set grid characteristic */
static void find_bh_surface_perfect_s2(Physics_T *const phys)
{
  Grid_Char_T *grid_char = phys->grid_char;
  const unsigned lmax   = (unsigned)Geti("surface_Ylm_expansion_max_l");
  const unsigned Ntheta = Ntheta_Ylm(lmax);
  const unsigned Nphi   = Nphi_Ylm(lmax);
  const unsigned Ntot   = Ntotal_Ylm(lmax);
  const double R_BH     = Getd("perfect_S2_radius");
  double *rbh = alloc_double(Ntot);/* surface function r = r(th,ph). */
  double *reClm_rbh = alloc_ClmYlm(lmax),
         *imClm_rbh = alloc_ClmYlm(lmax);
  unsigned ij;
  
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
    
}


/* adjust the BH radius to acquire the desired BH irreducible mass
// when the bh is a perfect sphere in x coords. */
static void tune_BH_radius_irreducible_mass_perfect_s2(Physics_T *const phys)
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

/* adding default parameters. */
int bh_add_params(Physics_T *const phys)
{
  FUNC_TIC
  
  UNUSED(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* adding fields. */
int bh_add_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  UNUSED(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}


