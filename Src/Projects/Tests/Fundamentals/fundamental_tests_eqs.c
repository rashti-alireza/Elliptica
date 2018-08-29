/*
// Alireza Rashti
// August 2018
*/

#include "fundamental_tests_eqs.h"

/* adding fundamental_tests equations to data base */
void fundamental_tests_fill_db_eqs(sEquation_T ***const field_eq,sEquation_T ***const bc_eq,sEquation_T ***const jacobian_eq)
{
  /* adding field and boundary condition equations */
  *field_eq      = init_eq();
  *bc_eq         = init_eq();
  *jacobian_eq   = init_eq();

  add_eq(field_eq,eq_alpha,"eq_alpha");
  add_eq(bc_eq,bc_alpha,"bc_alpha");
  add_eq(jacobian_eq,jacobian_alpha,"jacobian_alpha");
}

/* the following functions are just for test purposes */
static void *eq_alpha(void *vp1,void *vp2)
{
  UNUSED(vp1);
  UNUSED(vp2);
  return 0;
}

static void *bc_alpha(void *vp1,void *vp2)
{
  UNUSED(vp1);
  UNUSED(vp2);
  return 0;
}

static void *jacobian_alpha(void *vp1,void *vp2)
{
  UNUSED(vp1);
  UNUSED(vp2);
  return 0;
}
