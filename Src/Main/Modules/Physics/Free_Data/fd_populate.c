/*
// Alireza Rashti
// November 2020
*/

/* various functions and ways to populate free data */

#include "fd_populate.h"

/* compute K_{ij}, trK = ig^{ij} K_{ij} and its partial derivatives dtrK 
// for KerrSchild. */
void fd_extrinsic_curvature_KerrSchild(Physics_T *const phys,
                                         const char *const region,
                                         const char *const ig,
                                         const char *const Chris,
                                         const char *const Kij,
                                         const char *const trK,
                                         const char *const dtrK)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  const double BHx   = Getd("center_x");
  const double BHy   = Getd("center_y");
  const double BHz   = Getd("center_z");
  Uint p;
  
  fd_KerrSchild_set_params(phys);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    fd_Kij_trK_KerrSchild(patch,BHx,BHy,BHz,ig,Chris,Kij,trK);
    
    dField_di_STEM(dtrK_D0,dtrK);
    dField_di_STEM(dtrK_D1,dtrK);
    dField_di_STEM(dtrK_D2,dtrK);
    
  }
  
  FUNC_TOC
}

/* compute K_{ij}, trK = ig^{ij} K_{ij} and its partial derivatives dtrK 
// for Schwarzchild in isotropic coords. */
void fd_extrinsic_curvature_IsoSchild(Physics_T *const phys,
                                      const char *const region,
                                      const char *const ig,
                                      const char *const Chris,
                                      const char *const Kij,
                                      const char *const trK,
                                      const char *const dtrK)
{
  FUNC_TIC
  
  fd_extrinsic_curvature_Minkowski(phys,region,Kij,trK,dtrK);
  
  UNUSED(ig);
  UNUSED(Chris);
  
  FUNC_TOC
}

/* compute K_{ij}, trK = ig^{ij} K_{ij} and its partial derivatives dtrK 
// for Minkowski, they are all zero. */
void fd_extrinsic_curvature_Minkowski(Physics_T *const phys,
                                      const char *const region,
                                      const char *const Kij,
                                      const char *const trK,
                                      const char *const dtrK)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    /* they all are zero */
    REALLOC_v_WRITE_v_STEM(K_D2D2,Kij)
    REALLOC_v_WRITE_v_STEM(K_D0D2,Kij)
    REALLOC_v_WRITE_v_STEM(K_D0D0,Kij)
    REALLOC_v_WRITE_v_STEM(K_D0D1,Kij)
    REALLOC_v_WRITE_v_STEM(K_D1D2,Kij)
    REALLOC_v_WRITE_v_STEM(K_D1D1,Kij)
    
    REALLOC_v_WRITE_v_STEM(tr,trK)
    REALLOC_v_WRITE_v_STEM(dtr_D0,dtrK)
    REALLOC_v_WRITE_v_STEM(dtr_D1,dtrK)
    REALLOC_v_WRITE_v_STEM(dtr_D2,dtrK)
  }
  
  FUNC_TOC
}

/* compute confromal Ricci_{ij} and its trace */
void fd_conformal_Ricci(Physics_T *const phys,
                          const char *const region,
                          const char *const ig,
                          const char *const Chris,
                          const char *const dChris,
                          const char *const RicciConf,
                          const char *const trRicciConf)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Ricci_3d(patch,ig,Chris,dChris,RicciConf,trRicciConf);
  }
  
  FUNC_TOC
}

/* compute Christoffel symbol compatible with given metric */
void fd_compatible_Christoffel_symbol(Physics_T *const phys,
                                        const char *const region,
                                        const char *const ig,
                                        const char *const dg,
                                        const char *const Chris
                                       )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Christoffel_symbol_3d(patch,ig,dg,Chris);
  }
  
  FUNC_TOC
}

/* compute 1st derivative Christoffel symbol */
void fd_1st_derivative_Christoffel_symbol(Physics_T *const phys,
                                            const char *const region,
                                            const char *const dChris)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    dChristoffel_symbol_3d(patch,dChris);
  }
  
  FUNC_TOC
}

/* populate confromal metric, inverse of confromal metric 
// and first order derivative of confromal metric for KerrSchild.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_dgConf_igConf_KerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  const double BHx   = Getd("center_x");
  const double BHy   = Getd("center_y");
  const double BHz   = Getd("center_z");
  Uint p;
  
  fd_KerrSchild_set_params(phys);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    fd_kerr_schild_g_analytic(patch,BHx,BHy,BHz,gConf);
    fd_kerr_schild_dg_analytic(patch,BHx,BHy,BHz,dgConf);
    
    READ_v_STEM(gConf_D2D2,gConf)
    READ_v_STEM(gConf_D0D2,gConf)
    READ_v_STEM(gConf_D0D0,gConf)
    READ_v_STEM(gConf_D0D1,gConf)
    READ_v_STEM(gConf_D1D2,gConf)
    READ_v_STEM(gConf_D1D1,gConf)
    
    REALLOC_v_WRITE_v_STEM(igConf_U2U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U0,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U1,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U1U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U1U1,igConf)
    
    FOR_ALL_ijk
    {
      Matrix_Inverse_3x3_Symmetric_Field(gConf,D,igConf,U,ijk);
    }
  }
  
  FUNC_TOC
}

/* populate confromal metric, inverse of confromal metric 
// and first order derivative of confromal metric for Schwarzchild in
// isotropic coords.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_dgConf_igConf_IsoSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    char regex[99] = {'\0'};
    
    REALLOC_v_WRITE_v_STEM(gConf_D2D2,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D0D2,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D0D0,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D0D1,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D1D2,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D1D1,gConf)
    
    FOR_ALL_ijk
    {
      gConf_D0D0[ijk] = 1.;
      gConf_D1D1[ijk] = 1.;
      gConf_D2D2[ijk] = 1.;
    }
    /* since gConf is constant dgConf is machine precision exact */
    sprintf(regex,"^%s_D.D.D.$",dgConf);
    partial_derivative_with_regex(patch,regex);
    
    /* calculate the inverse of gConf */
    REALLOC_v_WRITE_v_STEM(igConf_U2U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U0,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U1,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U1U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U1U1,igConf)
    
    FOR_ALL_ijk
    {
      Matrix_Inverse_3x3_Symmetric_Field(gConf,D,igConf,U,ijk);
    }
  }
  
  FUNC_TOC
}

/* populate psi, alpha and beta KerrSchild value. */
void 
fd_populate_psi_alphaPsi_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig/*(inverse metric) if ig is null, it makes them */
 )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  const double BHx   = Getd("center_x");
  const double BHy   = Getd("center_y");
  const double BHz   = Getd("center_z");
  Uint p;
  
  fd_KerrSchild_set_params(phys);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    /* if there is no ig, make one */
    if (!ig)
    {
     /* populate g */
     ADD_AND_ALLOC_FIELD(fd_ks__g_D2D2)
     ADD_AND_ALLOC_FIELD(fd_ks__g_D0D2)
     ADD_AND_ALLOC_FIELD(fd_ks__g_D0D0)
     ADD_AND_ALLOC_FIELD(fd_ks__g_D0D1)
     ADD_AND_ALLOC_FIELD(fd_ks__g_D1D2)
     ADD_AND_ALLOC_FIELD(fd_ks__g_D1D1)
     
     fd_kerr_schild_g_analytic(patch,BHx,BHy,BHz,"fd_ks__g");
     
     READ_v_STEM(g_D2D2,"fd_ks__g")
     READ_v_STEM(g_D0D2,"fd_ks__g")
     READ_v_STEM(g_D0D0,"fd_ks__g")
     READ_v_STEM(g_D0D1,"fd_ks__g")
     READ_v_STEM(g_D1D2,"fd_ks__g")
     READ_v_STEM(g_D1D1,"fd_ks__g")
     
     /* populate ig */
     ADD_FIELD(fd_ks__ig_U2U2)
     ADD_FIELD(fd_ks__ig_U0U2)
     ADD_FIELD(fd_ks__ig_U0U0)
     ADD_FIELD(fd_ks__ig_U0U1)
     ADD_FIELD(fd_ks__ig_U1U2)
     ADD_FIELD(fd_ks__ig_U1U1)
     
     REALLOC_v_WRITE_v_STEM(ig_U2U2,"fd_ks__ig")
     REALLOC_v_WRITE_v_STEM(ig_U0U2,"fd_ks__ig")
     REALLOC_v_WRITE_v_STEM(ig_U0U0,"fd_ks__ig")
     REALLOC_v_WRITE_v_STEM(ig_U0U1,"fd_ks__ig")
     REALLOC_v_WRITE_v_STEM(ig_U1U2,"fd_ks__ig")
     REALLOC_v_WRITE_v_STEM(ig_U1U1,"fd_ks__ig")
     
     FOR_ALL_ijk
     {
       Matrix_Inverse_3x3_Symmetric_Field(g,D,ig,U,ijk);
     }
     fd_psi_alphaPsi_beta_KerrSchild_patch(patch,BHx,BHy,BHz,
                                       "fd_ks__ig",Psi,AlphaPsi,Beta);
     /* remove redundant */
     remove_field_with_regex(patch,"^fd_ks__g_D.+");
     remove_field_with_regex(patch,"^fd_ks__ig_U.+");
    }
    else
     fd_psi_alphaPsi_beta_KerrSchild_patch(patch,BHx,BHy,BHz,
                                        ig,Psi,AlphaPsi,Beta);
  }
  
  FUNC_TOC
}

/* populate psi, alpha and beta Schwarzchild in isotropic coords value. */
void 
fd_populate_psi_alphaPsi_beta_IsoSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig/*(inverse metric) if ig is null, it makes them */
 )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  const double   M   = Getd("irreducible_mass");
  const double BHx   = Getd("center_x");
  const double BHy   = Getd("center_y");
  const double BHz   = Getd("center_z");
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    REALLOC_v_WRITE_v_STEM(psi,Psi)
    REALLOC_v_WRITE_v_STEM(alphaPsi,AlphaPsi)
    
    /* beta identically zero */
    REALLOC_v_WRITE_v_STEM(beta_U0,Beta)
    REALLOC_v_WRITE_v_STEM(beta_U1,Beta)
    REALLOC_v_WRITE_v_STEM(beta_U2,Beta)
    
    FOR_ALL_ijk
    {
      double x,y,z,r;
      
      x = patch->node[ijk]->x[0]-BHx;
      y = patch->node[ijk]->x[1]-BHy;
      z = patch->node[ijk]->x[2]-BHz;
      r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      
      double alpha = (1-M/(2*r))/(1+M/(2*r));
      psi[ijk]     = (1+M/(2*r));
      alphaPsi[ijk]= alpha*psi[ijk];
      
    }
    UNUSED(beta_U0);
    UNUSED(beta_U1);
    UNUSED(beta_U2);
  }
  
  UNUSED(ig);
  FUNC_TOC
  
}
