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

/* setting grid characteristic */
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
  
  Psetd("min_radius",r_bh);
  Psetd("max_radius",r_bh);
  
  if (EQL(dr,0))/* => no change in AH surface */
  {
    Seti("did_AH_surface_change?",0);
  }
  else          /* => change in AH surface */
  {
    Seti("did_AH_surface_change?",1);
  }
  
}