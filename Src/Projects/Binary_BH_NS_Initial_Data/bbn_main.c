/*
// Alireza Rashti
// June 2019
*/

#include "bbn_main.h"

/* constructing initial data for system of binary black hole neutron star */
int Binary_BH_NS_Initial_Data(void)
{
  /* if this is a convergent test */
  if (strcmp_i(PgetsEZ("elliptic_convergence_test"),"yes"))
    bbn_elliptic_eqs_convergence_test();
    
  /* if this is a BAM call */
  else if (strcmp_i(PgetsEZ("bbn_bam_export_id"),"yes"))
    bbn_bam_export_id();
    
  /* construct id */
  else
    bbn_construct_id();
 
  return EXIT_SUCCESS;
}

