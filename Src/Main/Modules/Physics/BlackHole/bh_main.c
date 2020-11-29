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
  
  assert(phys->ctype == BH);
  
  switch (phys->cmd)
  {
    case BH_TUNE_RADIUS:
      ret = tune_black_hole_radius(phys);
    break;
    
    case BH_FIND_SURFACE:
      ret = find_black_hole_surface(phys);
    break;
    
    case BH_FILL:
      ret = bh_fill_inside_black_hole(phys);
    break;
    
    case BH_START:
      ret = start_off_black_hole(phys);
    break;
    
    case BH_ADD_PARAMS:
      ret = add_black_hole_params(phys);
    break;
    
    case BH_ADD_FIELDS:
      ret = add_black_hole_fields(phys);
    break;
    
    default:
      Error0(NO_OPTION);
  }
  
  return ret;
}

/* adding default parameters. */
static int add_black_hole_params(Physics_T *const phys)
{
  FUNC_TIC
  
  UNUSED(phys);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

/* adding fields. */
static int add_black_hole_fields(Physics_T *const phys)
{
  FUNC_TIC
  
  assert(phys->grid);
  
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
    else
      Error0(NO_OPTION);
  }
  else
    Error0(NO_OPTION);
  
  FUNC_TOC
  return EXIT_SUCCESS;
}

