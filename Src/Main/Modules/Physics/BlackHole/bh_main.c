/*
// Alireza Rashti
// November 2020
*/

/* black hole general affairs. one can add new different function readily
// by adding new parameter and the name of the function as shown. */

#include "bh_main.h"

/* main function to issue commands */
int bh_main(Physics_T *const phys)
{
  int ret = EXIT_SUCCESS;
  
  switch (phys->cmd)
  {
    case BH_TUNE_RADIUS:
      AssureType(phys->ctype == BH);
      ret = tune_black_hole_radius(phys);
    break;
    
    case BH_FIND_SURFACE:
      AssureType(phys->ctype == BH);
      ret = find_black_hole_surface(phys);
    break;
    
    case BH_FILL:
      AssureType(phys->ctype == BH);
      ret = bh_fill_inside_black_hole(phys);
    break;
    
    case BH_START:
      AssureType(phys->ctype == BH);
      ret = start_off_black_hole(phys);
    break;
    
    case BH_SET_PARAMS:
      ret = set_black_hole_params(phys);
    break;
    
    case BH_ADD_FIELDS:
      ret = add_black_hole_fields(phys);
    break;
    
    case BH_UPDATE_sConf:
      AssureType(phys->ctype == BH);
      ret = update_conformal_normal_vector_on_AH(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* set default parameters. */
static int set_black_hole_params(Physics_T *const phys)
{
  FUNC_TIC

  /* these parameters have a prefix BH.?_ and 
  // they are supposed to be added in parameter file:
  //
  // params:
  // =======
  //
  // surface_type:
  // 	o. perfect_s2: perfect sphere
  //    o. KerrSchild_s2: Kerr-Schild with arbitrary spin and boost
  //
  //
  // tune_BH_radius_criteria:
  // 	o. fix_irreducible_mass: keep irreducible mass constant
  //
  //
  // start_off:
  //	o. KerrSchild: a kerr-schild black hole
  //    o. IsoSchild : a Schwarzchild in isotropic coords.
  //    o. PGSchild  : a Schwarzchild in Painleve-Gullstrand coords. */
  
  
  UNUSED(phys);
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* adding fields. */
static int add_black_hole_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  bh_add_fields(phys->grid);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}


/* adjust AH radius to meet a criteria for instant the mass is fixed */
static int tune_black_hole_radius(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("tune_BH_radius_criteria","fix_irreducible_mass")
  {
    IF_sval("surface_type","perfect_s2")
      bh_tune_BH_radius_irreducible_mass_perfect_s2(phys);
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* find bh surface, mainly for grid setup */
static int find_black_hole_surface(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("surface_type","perfect_s2")
  {
    bh_find_bh_surface_perfect_s2(phys);
  }
  else IF_sval("surface_type","KerrSchild_s2")
  {
    bh_find_bh_surface_KerrSchild_s2(phys);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* set initial feature of blach holes, mostly used 
// for the very first time that we need to make grid. */
static int start_off_black_hole(Physics_T *const phys)
{
  FUNC_TIC
  
  IF_sval("start_off","KerrSchild")
  {
    IF_sval("surface_type","perfect_s2")
      bh_start_off_KerrSchild_perfect_s2(phys);
    else IF_sval("surface_type","KerrSchild_s2")
      bh_start_off_KerrSchild_general_s2(phys);
    else
      Error0(NO_OPTION);
  }
  else IF_sval("start_off","IsoSchild")
  {
    IF_sval("surface_type","perfect_s2")
      bh_start_off_IsoSchild_perfect_s2(phys);
    else
      Error0(NO_OPTION);
  }
  else IF_sval("start_off","PGSchild")
  {
    IF_sval("surface_type","perfect_s2")
      bh_start_off_PGSchild_perfect_s2(phys);
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* updating conformal normal and its derivatives on apparent horizon 
// note: we need gConf ready for this function. */
static int update_conformal_normal_vector_on_AH(Physics_T *const phys)
{
  FUNC_TIC
  
  bh_update_sConf_dsConf(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}
