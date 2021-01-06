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


/* compute beta and dbeta for Schwarzchild in Painleve-Gullstrand coords. 
// note: if dbeta is null, it won't compute derivatives. */
void fd_beta_and_dbeta_PGSchild(Physics_T *const phys,
                                const char *const region,
                                const char *const beta,
                                const char *const dbeta)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  const double M     = Getd("irreducible_mass");
  const double BHx   = Getd("center_x");
  const double BHy   = Getd("center_y");
  const double BHz   = Getd("center_z");
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Transformation_T *t = init_transformation();
    
    t->spheNcart->s2c = 1;/* spherical to Cartesian */
    
    REALLOC_v_WRITE_v_STEM(b_U0,beta)
    REALLOC_v_WRITE_v_STEM(b_U1,beta)
    REALLOC_v_WRITE_v_STEM(b_U2,beta)
    
    FOR_ALL_ijk
    {
      double x,y,z,r;
      double b_cart[3] = {0};
      double b_sphe[3] = {0};
      
      x = patch->node[ijk]->x[0]-BHx;
      y = patch->node[ijk]->x[1]-BHy;
      z = patch->node[ijk]->x[2]-BHz;
      r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      
      t->spheNcart->r  = r;
      t->spheNcart->th = acos(z/r);
      t->spheNcart->ph = arctan(y,x);
      
      /* beta in spherical coords */
      b_sphe[0] = sqrt(2*M/r);
      b_sphe[1] = 0.;
      b_sphe[2] = 0.;
      
      spheNcart_transformation(t,b_sphe,b_cart);
      
      /* beta in Cartesian coords */
      b_U0[ijk] = b_cart[0];
      b_U1[ijk] = b_cart[1];
      b_U2[ijk] = b_cart[2];
    }
    free_transformation(t);
    
    if (dbeta)
    {
      dField_di_STEM(dbeta_U0D0,dbeta);
      dField_di_STEM(dbeta_U0D1,dbeta);
      dField_di_STEM(dbeta_U0D2,dbeta);
      
      dField_di_STEM(dbeta_U1D0,dbeta);
      dField_di_STEM(dbeta_U1D1,dbeta);
      dField_di_STEM(dbeta_U1D2,dbeta);
      
      dField_di_STEM(dbeta_U2D0,dbeta);
      dField_di_STEM(dbeta_U2D1,dbeta);
      dField_di_STEM(dbeta_U2D2,dbeta);
    }
  }

  FUNC_TOC
}                                

/* compute K_{ij}, trK = ig^{ij} K_{ij} and its partial derivatives dtrK 
// for Schwarzchild in Painleve-Gullstrand coords. */
void fd_extrinsic_curvature_PGSchild(Physics_T *const phys,
                                      const char *const region,
                                      const char *const ig,
                                      const char *const Chris,
                                      const char *const Kij,
                                      const char *const trK,
                                      const char *const dtrK)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  /* add fields */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    /* beta */
    add_field("pg__beta_U0",0,patch,NO);
    add_field("pg__beta_U1",0,patch,NO);
    add_field("pg__beta_U2",0,patch,NO);
    
    
    /* add dbeta/dx? */
    add_field("dpg__beta_U0D0",0,patch,NO);
    add_field("dpg__beta_U0D1",0,patch,NO);
    add_field("dpg__beta_U0D2",0,patch,NO);
    
    add_field("dpg__beta_U1D0",0,patch,NO);
    add_field("dpg__beta_U1D1",0,patch,NO);
    add_field("dpg__beta_U1D2",0,patch,NO);
    
    add_field("dpg__beta_U2D0",0,patch,NO);
    add_field("dpg__beta_U2D1",0,patch,NO);
    add_field("dpg__beta_U2D2",0,patch,NO);
  }
  
  /* populate beta and dbeta for PG */
  fd_beta_and_dbeta_PGSchild(phys,region,"pg__beta","dpg__beta");
  
  /* populate Kij, trK and dtrK and clean up */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    fd_Kij_trK_PGSchild(patch,"dpg__beta",Kij,trK);
    dField_di_STEM(dtrK_D0,dtrK);
    dField_di_STEM(dtrK_D1,dtrK);
    dField_di_STEM(dtrK_D2,dtrK);
    
    /* clean up */
    remove_field_with_regex(patch,"^pg__beta_U.$");
    remove_field_with_regex(patch,"^dpg__beta_U.D.$");
  }
    
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
    
    /* UNUSED */
    UNUSED(K_D2D2)
    UNUSED(K_D0D2)
    UNUSED(K_D0D0)
    UNUSED(K_D0D1)
    UNUSED(K_D1D2)
    UNUSED(K_D1D1)
    
    UNUSED(tr)
    UNUSED(dtr_D0)
    UNUSED(dtr_D1)
    UNUSED(dtr_D2)
  }
  
  FUNC_TOC
}

/* compute conformal Ricci_{ij} and its trace */
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

/* populate conformal metric, inverse of conformal metric 
// and first order derivative of conformal metric for KerrSchild.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_igConf_dgConf_KerrSchild
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
    
    /* if asked for the metric */
    if (gConf)
      fd_kerr_schild_g_analytic(patch,BHx,BHy,BHz,gConf);
    
    /* if asked for the derivatives */
    if(dgConf)
      fd_kerr_schild_dg_analytic(patch,BHx,BHy,BHz,dgConf);
    
    
    /* if asked for the inverse */
    if (gConf && igConf)
    {
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
  }
  
  FUNC_TOC
}

/* populate conformal metric, inverse of conformal metric 
// and first order derivative of conformal metric for parameter
// "flat + exp(-r^4)*(KerrSchild-flat)"  which is:
// gConf_{ij} = delta_{ij} + atten * (gKS_{ij} - delta_{ij}).
// this free data mainly is used for BHNS system.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_igConf_dgConf_flat_expm4KS
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 )
{
  FUNC_TIC
  
  Grid_T *const grid  = mygrid(phys,region);
  Physics_T *const bh = init_physics(phys,BH);
  const double BHx    = Pgetd("BH_center_x");
  const double BHy    = Pgetd("BH_center_y");
  const double BHz    = Pgetd("BH_center_z");
  const double R04    = pow(Pgetd(P_"RollOff_radius"),4.);
  
  /* add some auxiliary fields */
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    ADD_AND_ALLOC_FIELD(fd_ks__g_D2D2)
    ADD_AND_ALLOC_FIELD(fd_ks__g_D0D2)
    ADD_AND_ALLOC_FIELD(fd_ks__g_D0D0)
    ADD_AND_ALLOC_FIELD(fd_ks__g_D0D1)
    ADD_AND_ALLOC_FIELD(fd_ks__g_D1D2)
    ADD_AND_ALLOC_FIELD(fd_ks__g_D1D1)
  }
  
  /* populate fd_ks__g */
  fd_populate_gConf_igConf_dgConf_KerrSchild(bh,".*","fd_ks__g",0,0);
  
  /* superimpose */
  OpenMP_Patch_Pragma(omp parallel for)
  for (Uint p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    READ_v(fd_ks__g_D2D2)
    READ_v(fd_ks__g_D0D2)
    READ_v(fd_ks__g_D0D0)
    READ_v(fd_ks__g_D0D1)
    READ_v(fd_ks__g_D1D2)
    READ_v(fd_ks__g_D1D1)
    
    REALLOC_v_WRITE_v_STEM(gConf_D2D2,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D0D2,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D0D0,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D0D1,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D1D2,gConf)
    REALLOC_v_WRITE_v_STEM(gConf_D1D1,gConf)
    
    /* g = delta_{ij} + atten * (gKS_{ij} - delta_{ij}) */
    FOR_ALL_ijk
    {
      double x   = patch->node[ijk]->x[0] - BHx;
      double y   = patch->node[ijk]->x[1] - BHy;
      double z   = patch->node[ijk]->x[2] - BHz;
      double r2  = (Pow2(x)+Pow2(y)+Pow2(z));
      double att = exp(-pow(r2,2.)/R04);
      
      /* diagonal */
      gConf_D0D0[ijk] = 1.+att*(fd_ks__g_D0D0[ijk]-1.);
      gConf_D1D1[ijk] = 1.+att*(fd_ks__g_D1D1[ijk]-1.);
      gConf_D2D2[ijk] = 1.+att*(fd_ks__g_D2D2[ijk]-1.);
      
      /* off diagonal */
      gConf_D0D1[ijk] = att*(fd_ks__g_D0D1[ijk]);
      gConf_D0D2[ijk] = att*(fd_ks__g_D0D2[ijk]);
      gConf_D1D2[ijk] = att*(fd_ks__g_D1D2[ijk]);
    }
  }
  
  /* populate derivatives and inverse */
  OpenMP_Patch_Pragma(omp parallel for)
  for (Uint p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    char regex[STR_LEN];
    
    sprintf(regex,"^%s_D.D.D.$",dgConf);
    partial_derivative_with_regex(patch,regex);
    
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
  
  /* remove added fields */
  for (Uint p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    remove_field_with_regex(patch,"^fd_ks__g__D.D.$");
  }
  free_physics(bh);
  
  FUNC_TOC
}

/* populate conformal metric, inverse of conformal metric 
// and first order derivative of conformal metric for 
// conformal KerrSchild, in this case the det(conformal metric) = 1
// and conformal factor is no longer 1 as we have for KerrSchild itself.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_igConf_dgConf_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *g/* given metric stem (if is null it makes it)*/,
 const char *const gConf/* stem of conformal metric */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  if (g)/* if no KerrSchild metric given, make it */
  {
    const double BHx   = Getd("center_x");
    const double BHy   = Getd("center_y");
    const double BHz   = Getd("center_z");
    
    fd_KerrSchild_set_params(phys);
    
    OpenMP_Patch_Pragma(omp parallel for)
    for (p = 0; p < grid->np; ++p)
    {
      Patch_T *patch = grid->patch[p];
      fd_kerr_schild_g_analytic(patch,BHx,BHy,BHz,gConf);
    }
    g = gConf;
  }
  
  /* scale it to have det = 1. */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    char regex[99] = {'\0'};
  
    scale_to_BSSN_metric_3d(patch,g,gConf,igConf,0);
    
    /* derivatives */
    sprintf(regex,"^%s_D.D.D.$",dgConf);
    partial_derivative_with_regex(patch,regex);
  }
  
  FUNC_TOC
}

/* populate conformal metric, inverse of conformal metric 
// and first order derivative of conformal metric for Schwarzchild in
// isotropic coords.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_igConf_dgConf_IsoSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 )
{
  FUNC_TIC
  
  fd_populate_gConf_igConf_dgConf_flat
   (phys,region,gConf,igConf,dgConf);
  
  FUNC_TOC
}

/* populate conformal metric, inverse of conformal metric 
// and first order derivative of conformal metric for Schwarzchild in
// Painleve-Gullstrand coords.
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_igConf_dgConf_PGSchild
 (
 Physics_T *const phys,
 const char *const region/* where computations take place */,
 const char *const gConf/* metric stem */,
 const char *const igConf/* inverse of metric stem */,
 const char *const dgConf/* derivative of metric stem */
 )
{
  FUNC_TIC
  
  fd_populate_gConf_igConf_dgConf_flat
   (phys,region,gConf,igConf,dgConf);
  
  FUNC_TOC
}

/* populate conformal metric, inverse of conformal metric 
// and first order derivative of conformal metric for flat space .
// the nomenclature of fields determined by the passed stems */
void 
fd_populate_gConf_igConf_dgConf_flat
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
    
    REALLOC_v_WRITE_v_STEM(igConf_U2U2,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U2,igConf) 
    REALLOC_v_WRITE_v_STEM(igConf_U0U0,igConf)
    REALLOC_v_WRITE_v_STEM(igConf_U0U1,igConf) 
    REALLOC_v_WRITE_v_STEM(igConf_U1U2,igConf) 
    REALLOC_v_WRITE_v_STEM(igConf_U1U1,igConf)
    
    FOR_ALL_ijk
    {
      gConf_D0D0[ijk] = 1.;
      gConf_D1D1[ijk] = 1.;
      gConf_D2D2[ijk] = 1.;
      
      igConf_U0U0[ijk] = 1.;
      igConf_U1U1[ijk] = 1.;
      igConf_U2U2[ijk] = 1.;
      
    }
    UNUSED(igConf_U0U2);
    UNUSED(igConf_U0U1);
    UNUSED(igConf_U1U2);
    
    UNUSED(gConf_D0D2);
    UNUSED(gConf_D0D1);
    UNUSED(gConf_D1D2);
    
    /* since gConf is constant dgConf is machine precision exact */
    sprintf(regex,"^%s_D.D.D.$",dgConf);
    partial_derivative_with_regex(patch,regex);
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


/* populate beta ConfKerrSchild value. */
void 
fd_populate_beta_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Beta
 )
{
 FUNC_TIC
 
 fd_populate_beta_KerrSchild(phys,region,Beta);
 
 FUNC_TOC
}

/* populate beta KerrSchild value. */
void 
fd_populate_beta_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Beta
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
    fd_beta_KerrSchild_patch(patch,BHx,BHy,BHz,"fd_ks__ig",Beta);
    
    /* remove redundant */
    remove_field_with_regex(patch,"^fd_ks__g_D.+");
    remove_field_with_regex(patch,"^fd_ks__ig_U.+");
  }
  
  FUNC_TOC
}


/* populate psi, alpha and beta conformal KerrSchild value. */
void 
fd_populate_psi_alphaPsi_beta_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta
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
    
    /* make a KerrSchild metric */
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
    /* populate fields using KerrSchild */
    fd_psi_alphaPsi_beta_KerrSchild_patch
     (patch,BHx,BHy,BHz,"fd_ks__ig",Psi,AlphaPsi,Beta);
    
    /* rescale psi and alphaPsi */
    ADD_AND_ALLOC_FIELD(fd_ks__psi)
    scale_to_BSSN_metric_3d(patch,"fd_ks__g",0,0,"fd_ks__psi");
    
    WRITE_v_STEM(psi,Psi)
    WRITE_v_STEM(alphaPsi,AlphaPsi)
    READ_v_STEM(scaled_psi,"fd_ks__psi");
    
    FOR_ALL_ijk
    {
      psi[ijk]      *= scaled_psi[ijk];
      alphaPsi[ijk] *= scaled_psi[ijk];
    }
    
    /* remove redundant */
    remove_field_with_regex(patch,"^fd_ks__g_D.+");
    remove_field_with_regex(patch,"^fd_ks__ig_U.+");
    remove_field_with_regex(patch,"^fd_ks__psi$");
  }
  
  FUNC_TOC
}


/* populate alpha conformal KerrSchild value.
// it's the same as fd_populate_alpha_KerrSchild */
void 
fd_populate_alpha_ConfKerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Alpha
 )
{
  FUNC_TIC
  fd_populate_alpha_KerrSchild(phys,region,Alpha);
  FUNC_TOC
}

/* populate alpha of KerrSchild value. */
void 
fd_populate_alpha_KerrSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Alpha
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
    fd_alpha_KerrSchild_patch(patch,BHx,BHy,BHz,Alpha);
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
 const char *const ig
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

/* populate psi, alpha and beta Schwarzchild in Painleve-Gullstrand 
// coords value. */
void 
fd_populate_psi_alphaPsi_beta_PGSchild
 (
 Physics_T *const phys,
 const char *const region,
 const char *const Psi,
 const char *const AlphaPsi,
 const char *const Beta,
 const char *const ig
 )
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  /* beta */
  fd_beta_and_dbeta_PGSchild(phys,region,Beta,0);
  
  /* psi and alphaPsi */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    REALLOC_v_WRITE_v_STEM(psi,Psi)
    REALLOC_v_WRITE_v_STEM(alphaPsi,AlphaPsi)
    
    FOR_ALL_ijk
    {
      psi[ijk]     = 1.;
      alphaPsi[ijk]= 1.;
    }
  }
  
  UNUSED(ig);
  FUNC_TOC
  
}
