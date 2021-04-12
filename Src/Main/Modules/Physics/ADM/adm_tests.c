/*
// Alireza Rashti
// December 2020
*/


/* various doc test to make sure things works correctly */

#include "adm_tests.h"


/* test Aconf^{ij} and adm_Kij
// NOTE: it overwrites many fields thus fields are not supposed to be 
// used again and one must exit from project or reset fields.
// NOTE: It assumes AConf^{ij} are already populated. */
void adm_doctest_AConfIJ(Physics_T *const phys)
{
  FUNC_TIC
  
  if (phys->sys                       == SBH         &&
      Pcmps(P_"doctest_AConfIJ_compare","KerrSchild"))
  {
    
    /* calculate analytic adm_K_{ij} using beta name test_adm_Kij */
    if (1)
    {
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      Grid_T *const grid  = mygrid(phys,".*");

      add_3x3_symmetric_field(grid,"test_adm_Kij","down");
      
      fd_populate_gConf_igConf_dgConf_KerrSchild(bh,".*","gConf",
                                                  "igConf","dgConf");
                                                  
      fd_1st_derivative_Christoffel_symbol(bh,".*","dChrisConf");
      
      fd_conformal_Ricci(bh,".*","igConf","ChrisConf","dChrisConf",
                           "RicciConf","trRicciConf");
      
      fd_extrinsic_curvature_KerrSchild(bh,".*","igConf","ChrisConf",
                                          "test_adm_Kij","trK","dtrK");
      
      /* compute adm_K_{ij} using AConf^{ij} */
      physics(phys,ADM_UPDATE_Kij);
      
      /* compare */
      diff_3x3_symmetric_fields
        (grid,"test_adm_Kij","adm_Kij","down",1);
      
      /* remove test_adm_Kij */
      FOR_ALL_p(grid->np)
      {
        Patch_T *patch = grid->patch[p];
        remove_field_regex(patch,"^test_adm_Kij_D.*");
      }
      free_physics(bh);
    }
    
    /* calculate analytic K_{ij} using gr-qc/9805023.pdf */
    if (1)
    {
      /* important to have dedicated BH physics to read correct parameters */
      Physics_T *const bh = init_physics(phys,BH);
      KerrSchild_analytic_adm_Kij_comparison(bh);
      free_physics(bh);
    }
  }
  else
    Error0(NO_OPTION);
    
  FUNC_TOC
}


/* compare analytic and numeric value for adm_K_{ij} for
// a KerrSchild BH with NO spin and boost using gr-qc/9805023.pdf */
static void KerrSchild_analytic_adm_Kij_comparison(Physics_T *const phys)
{
  const int pr     = 0;
  const double M   = Getd("irreducible_mass");
  const double BHx = Getd("center_x");
  const double BHy = Getd("center_y");
  const double BHz = Getd("center_z");
  const double Cx  = Getd("chi_x");
  const double Cy  = Getd("chi_y");
  const double Cz  = Getd("chi_z");
  const double Bx  = Getd("boost_Vx");
  const double By  = Getd("boost_Vy");
  const double Bz  = Getd("boost_Vz");
  Grid_T *const grid = mygrid(phys,".*");
  double max = 0.;
  
  if (GRT(Cx,0.) ||
      GRT(Cy,0.) ||
      GRT(Cz,0.) ||
      GRT(Bx,0.) ||
      GRT(By,0.) ||
      GRT(Bz,0.))
  {
    printf(Pretty0"Spin or boost > 0 => quit the test.\n");
    return;
  }
  
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    const double KD[2] = {0,1};
    READ_v(adm_Kij_D2D2)
    READ_v(adm_Kij_D0D1)
    READ_v(adm_Kij_D0D0)
    READ_v(adm_Kij_D0D2)
    READ_v(adm_Kij_D1D1)
    READ_v(adm_Kij_D1D2)
    
    FOR_ALL_ijk
    {
      double x,y,z,r,r2;
      double fac;
      double Kxx,Kxy,Kxz,Kyy,Kyz,Kzz;
      double dKxx,dKxy,dKxz,dKyy,dKyz,dKzz;
      
      x = patch->node[ijk]->x[0]-BHx;
      y = patch->node[ijk]->x[1]-BHy;
      z = patch->node[ijk]->x[2]-BHz;
      r2= (Pow2(x)+Pow2(y)+Pow2(z));
      r = sqrt(r2);
      
      fac = 2*M/Pow2(r2)/sqrt(1+2*M/r);
      Kxx = fac*(r2*KD[1]-(2+M/r)*x*x);
      Kxy = fac*(r2*KD[0]-(2+M/r)*x*y);
      Kxz = fac*(r2*KD[0]-(2+M/r)*x*z);
      Kyy = fac*(r2*KD[1]-(2+M/r)*y*y);
      Kyz = fac*(r2*KD[0]-(2+M/r)*y*z);
      Kzz = fac*(r2*KD[1]-(2+M/r)*z*z);
      
      dKxx = fabs(Kxx-adm_Kij_D0D0[ijk]);
      dKxy = fabs(Kxy-adm_Kij_D0D1[ijk]);
      dKxz = fabs(Kxz-adm_Kij_D0D2[ijk]);
      dKyy = fabs(Kyy-adm_Kij_D1D1[ijk]);
      dKyz = fabs(Kyz-adm_Kij_D1D2[ijk]);
      dKzz = fabs(Kzz-adm_Kij_D2D2[ijk]);
      
      Is_Different(dKxx,pr);
      Is_Different(dKxy,pr);
      Is_Different(dKxz,pr);
      Is_Different(dKyy,pr);
      Is_Different(dKyz,pr);
      Is_Different(dKzz,pr);
    }
  }
  printf(Pretty0"Max difference = %e\n",max);
}

