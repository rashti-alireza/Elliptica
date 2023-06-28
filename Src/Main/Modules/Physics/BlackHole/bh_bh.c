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

/* use close to a KerrSchild BH assuming pefect S2 to start off 
// some black hole paramters and domain shape etc. 
// Note: this mainly used to initiate a binary system. */
void bh_start_off_CloseKerrSchild_perfect_s2(Physics_T *const phys)
{
  FUNC_TIC
  
  const double Fac         = 1.5;/* this is an experimental factor */
  const double bh_chi_x    = Getd("chi_x");
  const double bh_chi_y    = Getd("chi_y");
  const double bh_chi_z    = Getd("chi_z");
  const double bh_irr_mass = Getd("irreducible_mass");
  const double bh_chi = sqrt(Pow2(bh_chi_x)+Pow2(bh_chi_y)+Pow2(bh_chi_z));
  const double bh_chr_mass = bh_irr_mass*sqrt(2./(1.+sqrt(1-Pow2(bh_chi))));
  const double bh_R        = Fac*bh_irr_mass;/* approximate initial radius */
  const double bh_a        = bh_chi*bh_chr_mass;
  
  /* check size of bh_chi */
  if (GRT(bh_chi,1.))
    Errors("%s spin is too large!\n",phys->stype);
  
  /* set initial params */
  Setd("perfect_S2_radius",bh_R);
  Setd("min_radius",bh_R);
  Setd("max_radius",bh_R);
  Setd("Christodoulou_mass",bh_irr_mass);
  Setd("spin_a",bh_a);
  
  printf(Pretty0"%s properties:\n",phys->stype);
  printf(Pretty0"%s radius (Kerr-Schild Coords.) = %+e\n",phys->stype,bh_R);
  printf(Pretty0"%s irreducible mass             = %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s Christodoulou_mass           ~ %+e\n",phys->stype,bh_chr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,bh_chi_x);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,bh_chi_y);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,bh_chi_z);
  printf(Pretty0"%s net dimensionless spin       = %+e\n",phys->stype,bh_chi);
  printf(Pretty0"%s spin/M (= a)                 = %+e\n",phys->stype,bh_a);
  
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
  Setd("Christodoulou_mass",bh_chr_mass);
  Setd("spin_a",bh_a);
  
  printf(Pretty0"%s properties:\n",phys->stype);
  printf(Pretty0"%s radius (Kerr-Schild Coords.) ~ %+e\n",phys->stype,bh_R);
  printf(Pretty0"%s irreducible mass             ~ %+e\n",phys->stype,bh_irr_mass);
  printf(Pretty0"%s Christodoulou_mass           ~ %+e\n",phys->stype,bh_chr_mass);
  printf(Pretty0"%s dimensionless spin (x comp.) = %+e\n",phys->stype,bh_chi_x);
  printf(Pretty0"%s dimensionless spin (y comp.) = %+e\n",phys->stype,bh_chi_y);
  printf(Pretty0"%s dimensionless spin (z comp.) = %+e\n",phys->stype,bh_chi_z);
  printf(Pretty0"%s net dimensionless spin       = %+e\n",phys->stype,bh_chi);
  printf(Pretty0"%s spin/M (= a)                 = %+e\n",phys->stype,bh_a);
  
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
  const double R_BH = Getd("perfect_S2_radius");
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
  double obs_mass[2] = {0.};
  double dr,r_bh,dM,irr_mass;
  
  observe(phys,"Irreducible(M)",Gets("Observe_Irreducible_M"),obs_mass);
  irr_mass = obs_mass[0];
  Setd("irreducible_mass_current",irr_mass);
  
  printf(Pretty0"current BH irreducible mass = %e\n",irr_mass);
  printf(Pretty0"update weight               = %e\n",W);
  printf(Pretty0"mass tolerance              = %e\n",dM_tolerance);
  
  dM = fabs(irr_mass/target_bh_mass-1);
  dr = -current_r_bh*(irr_mass/target_bh_mass-1);
  if (EQL(W,0))
  {
    dr = 0;
  }
  if (LSSEQL(dM,dM_tolerance)) 
  {
    dr = 0;
    printf(Pretty0"|dM/M| = %g < Tol. = %g\n",dM,dM_tolerance);
  }
  r_bh = current_r_bh + W*dr;
  
  Setd("min_radius",r_bh);
  Setd("max_radius",r_bh);
  Setd("perfect_S2_radius",r_bh);
  
  if (EQL(dr,0))/* => no change in AH surface */
  {
    Seti("did_BH_surface_change?",0);
  }
  else          /* => change in AH surface */
  {
    Seti("did_BH_surface_change?",1);
  }
  
}

/* adjust BH omega to meet BH's chi (dimensionless spin) target value */
void bh_tune_BH_chi_simple(Physics_T *const phys)
{
  const double W1        = Getd("spin_update_weight");
  const double dchi_tol  = Getd("spin_tolerance");
  const double chi_x     = Getd("chi_x");
  const double chi_y     = Getd("chi_y");
  const double chi_z     = Getd("chi_z");
  const double omega_x   = Getd("Omega_x");
  const double omega_y   = Getd("Omega_y");
  const double omega_z   = Getd("Omega_z");
  const double omega_s   = fabs(sysGetd("angular_velocity"));/* scale */
  double omega[3]        = {0.};
  double s[3]            = {0.};
  double chi_current[3]  = {0.};
  double dchi[3]         = {0.};
  
  observe(phys,"spin",Gets("Observe_spin"),s);

  /* calculate BH current Christodoulou mass.
  // NOTE: s[?] depend on BH spin on top */
  double irr_mass = Getd("irreducible_mass_current");
  double net_spin = sqrt(Pow2(s[0])+Pow2(s[1])+Pow2(s[2]));
  double m = sqrt(Pow2(irr_mass)+Pow2(net_spin)/(4*Pow2(irr_mass)));

  /* current chi */
  chi_current[0] = s[0]/Pow2(m);
  chi_current[1] = s[1]/Pow2(m);
  chi_current[2] = s[2]/Pow2(m);
  
  dchi[0] = chi_current[0] - chi_x;
  dchi[1] = chi_current[1] - chi_y;
  dchi[2] = chi_current[2] - chi_z;
  
  printf(Pretty0"current %s chi_x = %e\n",phys->stype,chi_current[0]);
  printf(Pretty0"current %s chi_y = %e\n",phys->stype,chi_current[1]);
  printf(Pretty0"current %s chi_z = %e\n",phys->stype,chi_current[2]);
  printf(Pretty0"update weight    = %e\n",W1);
  printf(Pretty0"chi tolerance    = %e\n",dchi_tol);
  
  if (EQL(W1,0.))
  {
    return;
  }
  if (GRT(fabs(dchi[0]),dchi_tol))
  {
    omega[0] = omega_x+W1*dchi[0]*omega_s;
    Setd("Omega_x",omega[0]);
  }
  if (GRT(fabs(dchi[1]),dchi_tol))
  {
    omega[1] = omega_y+W1*dchi[1]*omega_s;
    Setd("Omega_y",omega[1]);
  }
  if (GRT(fabs(dchi[2]),dchi_tol))
  {
    omega[2] = omega_z+W1*dchi[2]*omega_s;
    Setd("Omega_z",omega[2]);
  }
}

/* update inner BC values, for example the values of alpha on AH */
void bh_update_inner_BC(Physics_T *const phys)
{
  
  IF_sval("Eq_inner_BC_fields","none")
  {
    ;
  }
  else IF_sval("Eq_inner_BC_fields","XCTS")
  {
    /* alpha field */
    IF_sval("Eq_inner_BC_alpha","exact_ConfKerrSchild")
    {
      fd_populate_alpha_ConfKerrSchild(phys,"BH_around_IB","ibc_alpha");
    }
    else IF_sval("Eq_inner_BC_alpha","exact_KerrSchild")
    {
      fd_populate_alpha_KerrSchild(phys,"BH_around_IB","ibc_alpha");
    }
    else IF_sval("Eq_inner_BC_alpha","w*KerrSchild")
    {
      fd_populate_alpha_wKerrSchild(phys,"BH_around_IB","ibc_alpha");
    }
    else IF_sval("Eq_inner_BC_alpha","none")
    {
      ;
    }
    else
    {
      Error0(NO_OPTION);
    }
    
    /* beta fields */
    IF_sval("Eq_inner_BC_beta","exact_KerrSchild")
    {
      fd_populate_beta_KerrSchild(phys,"BH_around_IB","ibc_beta");
    }
    else IF_sval("Eq_inner_BC_beta","exact_ConfKerrSchild")
    {
      fd_populate_beta_ConfKerrSchild(phys,"BH_around_IB","ibc_beta");
    }
    else IF_sval("Eq_inner_BC_beta","alpha+Omega*r")
    {
      set_beta_inner_bc_alpha_omegaXr(phys,"BH_around_IB","ibc_beta");
    }
    else IF_sval("Eq_inner_BC_beta","none")
    {
      ;
    }
    else
    {
      Error0(NO_OPTION);
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
}


/* updating conformal normal vector and its derivatives on 
// apparent horizon. 
// NOTE: bh_sConf^i = psi^2 bh_s^i and bh_sConf_i = psi^-2 bh_s_i.
// NOTE: bh_sConf^i is independent of psi.
// NOTE: it's difficult to resolve derivatives of normal vector in 
// cubed spherical coords. , thus, one must to resolve this 
// analytically or using another coords. system. */
void bh_update_sConf_dsConf(Physics_T *const phys)
{
  /* it assumes perfect S2 to use analytic formulae */
  if (!strcmp_i(Gets("surface_type"),"perfect_s2"))
    Warning("this function only supports 'surface_type = perfect_s2'!");
  
  Grid_T *const grid = mygrid(phys,Ftype("BH_around_IB"));
  const double BH_center_x = Getd("center_x");
  const double BH_center_y = Getd("center_y");
  const double BH_center_z = Getd("center_z");
  const double kd[2]       = {0.,1.};/* Kronecker Delta */
  
  add_aux_fields(grid,"bh__n,dbh__n_D0,dbh__n_D1,dbh__n_D2");
  
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    double n[3] = {0};/* normal of perfect sphere */
    double dn[3][3] = {{0}};/* derivative of n, dn^i/dx^j */
    double *const n_U0 = &n[0];
    double *const n_U1 = &n[1];
    double *const n_U2 = &n[2];
      
    /* populate norm: */  
    
    /* conformal metric */
    READ_v(gConf_D2D2)
    READ_v(gConf_D0D2)
    READ_v(gConf_D0D0)
    READ_v(gConf_D0D1)
    READ_v(gConf_D1D2)
    READ_v(gConf_D1D1)
    
    /* norm */
    REALLOC_v_WRITE_v(bh__n);
    
    FOR_ALL_ijk
    {
      double x = patch->node[ijk]->x[0]-BH_center_x;
      double y = patch->node[ijk]->x[1]-BH_center_y;
      double z = patch->node[ijk]->x[2]-BH_center_z;
      double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      
      /* n = x^i/r */
      n[0] = x/r;
      n[1] = y/r;
      n[2] = z/r;
      
      /* norm^2: used cpi: gConf(-i,-j)*n(i)*n(j); */  
      double N2 = 
 gConf_D0D0[ijk]*pow(n_U0[0], 2) + 2.0*gConf_D0D1[ijk]*n_U0[0]*n_U1[0] +
 2.0*gConf_D0D2[ijk]*n_U0[0]*n_U2[0] + gConf_D1D1[ijk]*pow(n_U1[0], 2) +
 2.0*gConf_D1D2[ijk]*n_U1[0]*n_U2[0] + gConf_D2D2[ijk]*pow(n_U2[0], 2);

      bh__n[ijk] = sqrt(N2);
    }
    
    /* populate sConf and dsConf: */
    
    /* normal derivatives */
    dField_di(dbh__n_D0);
    dField_di(dbh__n_D1);
    dField_di(dbh__n_D2);
    READ_v(dbh__n_D0);
    READ_v(dbh__n_D1);
    READ_v(dbh__n_D2);
    
    /* normal vector on horizon */
    REALLOC_v_WRITE_v(bh_sConf_U0);
    REALLOC_v_WRITE_v(bh_sConf_U1);
    REALLOC_v_WRITE_v(bh_sConf_U2);
    
    /* derivatives of normal vector on horizon */
    REALLOC_v_WRITE_v(dbh_sConf_U0D0);
    REALLOC_v_WRITE_v(dbh_sConf_U0D1);
    REALLOC_v_WRITE_v(dbh_sConf_U0D2);
    REALLOC_v_WRITE_v(dbh_sConf_U1D0);
    REALLOC_v_WRITE_v(dbh_sConf_U1D1);
    REALLOC_v_WRITE_v(dbh_sConf_U1D2);
    REALLOC_v_WRITE_v(dbh_sConf_U2D0);
    REALLOC_v_WRITE_v(dbh_sConf_U2D1);
    REALLOC_v_WRITE_v(dbh_sConf_U2D2);
    
    FOR_ALL_ijk
    {
      double x = patch->node[ijk]->x[0]-BH_center_x;
      double y = patch->node[ijk]->x[1]-BH_center_y;
      double z = patch->node[ijk]->x[2]-BH_center_z;
      double r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));
      
      /* n = x^i/r */
      n[0] = x/r;
      n[1] = y/r;
      n[2] = z/r;
      
      /* normalize */
      bh_sConf_U0[ijk] = n[0]/bh__n[ijk];
      bh_sConf_U1[ijk] = n[1]/bh__n[ijk];
      bh_sConf_U2[ijk] = n[2]/bh__n[ijk];
      
      /* dn^i/dx^j = (kd^{ij}-n^i n^j)/r */
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          dn[i][j] = (kd[i==j]-n[i]*n[j])/r;
      
      dbh_sConf_U0D0[ijk] = (dn[0][0] - n[0]*dbh__n_D0[ijk]/bh__n[ijk])/bh__n[ijk];
      dbh_sConf_U0D1[ijk] = (dn[0][1] - n[0]*dbh__n_D1[ijk]/bh__n[ijk])/bh__n[ijk];
      dbh_sConf_U0D2[ijk] = (dn[0][2] - n[0]*dbh__n_D2[ijk]/bh__n[ijk])/bh__n[ijk];
      
      dbh_sConf_U1D0[ijk] = (dn[1][0] - n[1]*dbh__n_D0[ijk]/bh__n[ijk])/bh__n[ijk];
      dbh_sConf_U1D1[ijk] = (dn[1][1] - n[1]*dbh__n_D1[ijk]/bh__n[ijk])/bh__n[ijk];
      dbh_sConf_U1D2[ijk] = (dn[1][2] - n[1]*dbh__n_D2[ijk]/bh__n[ijk])/bh__n[ijk];
    
      dbh_sConf_U2D0[ijk] = (dn[2][0] - n[2]*dbh__n_D0[ijk]/bh__n[ijk])/bh__n[ijk];
      dbh_sConf_U2D1[ijk] = (dn[2][1] - n[2]*dbh__n_D1[ijk]/bh__n[ijk])/bh__n[ijk];
      dbh_sConf_U2D2[ijk] = (dn[2][2] - n[2]*dbh__n_D2[ijk]/bh__n[ijk])/bh__n[ijk];
    }
  }
  
  remove_aux_fields(grid,"bh__n,dbh__n_D0,dbh__n_D1,dbh__n_D2");
}

/* set inner BC for beta such that:
// beta^i = alpha*s^i + Omega x r. */
static void 
set_beta_inner_bc_alpha_omegaXr
 (
 Physics_T *const phys,
 const char *const region,
 const char *const ib_Beta
 )
{
 FUNC_TIC
 
 AssureType(phys->ctype == BH)
 
 Grid_T *const grid = mygrid(phys,region);
 const double BH_center_x = Getd("center_x");
 const double BH_center_y = Getd("center_y");
 const double BH_center_z = Getd("center_z");
 const double BH_Omega_U0 = Getd("Omega_x");
 const double BH_Omega_U1 = Getd("Omega_y");
 const double BH_Omega_U2 = Getd("Omega_z");
 
 FOR_ALL_p(grid->np)
 {
   Patch_T *patch = grid->patch[p];
   
   READ_v(alphaPsi)
   READ_v(psi)
   READ_v(bh_sConf_U0)
   READ_v(bh_sConf_U1)
   READ_v(bh_sConf_U2)
   REALLOC_v_WRITE_v_STEM(ibc_U0,ib_Beta)
   REALLOC_v_WRITE_v_STEM(ibc_U1,ib_Beta)
   REALLOC_v_WRITE_v_STEM(ibc_U2,ib_Beta)
   
   FOR_ALL_ijk
   {
     double x = patch->node[ijk]->x[0]-BH_center_x;
     double y = patch->node[ijk]->x[1]-BH_center_y;
     double z = patch->node[ijk]->x[2]-BH_center_z; 
     
     double alpha = alphaPsi[ijk]/psi[ijk];
     double S_U0  = bh_sConf_U0[ijk]/pow(psi[ijk], 2);
     double S_U1  = bh_sConf_U1[ijk]/pow(psi[ijk], 2);
     double S_U2  = bh_sConf_U2[ijk]/pow(psi[ijk], 2);
     
     double OmegaXr_U0 = BH_Omega_U1*z - BH_Omega_U2*y;
     double OmegaXr_U1 = -BH_Omega_U0*z + BH_Omega_U2*x;
     double OmegaXr_U2 = BH_Omega_U0*y - BH_Omega_U1*x;

     ibc_U0[ijk] =  OmegaXr_U0 + S_U0*alpha;
     ibc_U1[ijk] =  OmegaXr_U1 + S_U1*alpha;
     ibc_U2[ijk] =  OmegaXr_U2 + S_U2*alpha;
   }
 }
 
 FUNC_TOC
}
