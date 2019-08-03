/*
// Alireza Rashti
// June 2019
*/

#include "bbn_initialize.h"

/* initialize this system according to the previous grid. 
// ->return value: the next grid as a result of this initialization. */
Grid_T *bbn_initialize_next_grid(Grid_T *const grid_prev)
{
  Grid_T *grid_next = 0;
  
  if (!grid_prev)/* if grid is empty come up with an approximation */
  {
    /* if we use TOV and Kerr-Schil black hole approximation */
    if (strcmp_i(GetParameterS_E("BH_NS_initialization"),"TOV_KerrShild"))
      grid_next = TOV_KerrShild_approximation();
    else
      abortEr(NO_OPTION);
  }
  else/* use previous grid to make the next one */
  {
    printf("Not made yet!\n");
  }
  
  return grid_next;   
}

/* use TOV and Kerr-Schil black hole approximation.
// ->return value: resultant grid from this approximation */
static Grid_T *TOV_KerrShild_approximation(void)
{
  Grid_T *grid = 0;
  
  /* solve fields for a TOV star located at left side of y axis */
  TOV_T *tov = TOV_init();
  tov->bar_m = GetParameterD_E("NS_baryonic_mass");
  tov->description = "Estimating NS";
  tov = TOV_solution(tov);
  const double ns_R = tov->rbar[tov->N-1];
  
  /* basics of Kerr Shild black hole located at right side of y axis */
  pr_line_custom('=');
  printf("Acquiring Black Hole properties ...\n");
  const double bh_chi  = GetParameterD_E("BH_dimensionless_spin");
  const double bh_mass = GetParameterD_E("BH_mass");
  const double bh_R    = bh_mass*(1+sqrt(1-SQR(bh_chi)));
  printf("BH properties:\n");
  printf("--> BH radius (Kerr-Schild Coords.) = %e\n",bh_R);
  printf("--> BH dimensionless spin           = %e\n",bh_chi);
  printf("--> BH ADM mass                     = %e\n",bh_mass);
  printf("Acquiring Black Hole properties ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
 
  /* combining these two geometry to create the grid */
  grid = creat_grid_TOV_KerrShild(ns_R,bh_R,bh_chi*bh_mass/* a = chi*M */);
  
  /* creating all of the fields needed for construction of Initial Data */
  bbn_allocate_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid);
  
  /* initialize the fields using TOV and Kerr-Shild solution */
  init_field_TOV_plus_KerrSchild(grid,tov,bh_chi*bh_mass,bh_mass);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid);
  
  /* updating u0, _J^i, _E and _S */
  Tij_IF_build_psi6Sources(grid);
  
  TOV_free(tov);
  
  return grid;
}

/* initialize the fields using TOV and Kerr-Shild solution.
// the idea is to superimpose two fields of each solutions. */
static void init_field_TOV_plus_KerrSchild(Grid_T *const grid,const TOV_T *const tov, const double a_BH, const double M_BH)
{
  pr_line_custom('=');
  printf("Initializing the fields using TOV and Kerr-Schild solution ...\n");

  const double M_NS = tov->ADM_m;/* NS adm mass */
  const double D = GetParameterD_E("BH_NS_separation");
  const double C_BH = 0.5*GetParameterD_E("BH_NS_separation");/* center of BH it's on +y axis */
  const double C_NS = -C_BH;/* center of NS it's on -y axis*/
  const double R_Schwar = tov->r[tov->N-1];/* NS's Schwarzchild radius */
  const double a2_BH = SQR(a_BH);/* spin vector of BH */
  const double y_CM = (M_NS*C_NS+M_BH*C_BH)/(M_NS+M_BH);/* center of rotation, approx. Center of Mass */
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double Omega_NS_x = GetParameterD_E("NS_spin_vector_x");
  const double Omega_NS_y = GetParameterD_E("NS_spin_vector_y");
  const double Omega_NS_z = GetParameterD_E("NS_spin_vector_z");
  const double Vr = GetParameterD_E("BH_NS_infall_velocity");
  unsigned p;
  
  add_parameter_double("NS_Center",C_NS);
  add_parameter_double("y_CM",y_CM);
  
  /* black hole parts */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    GET_FIELD(Beta_U0)
    GET_FIELD(Beta_U1)
    GET_FIELD(Beta_U2)
    GET_FIELD(_gammaI_U0U2)
    GET_FIELD(_gammaI_U0U0)
    GET_FIELD(_gammaI_U0U1)
    GET_FIELD(_gammaI_U1U2)
    GET_FIELD(_gammaI_U1U1)
    GET_FIELD(_gammaI_U2U2)

    ADD_FIELD(KSbeta_D0)
    ADD_FIELD(KSbeta_D1)
    ADD_FIELD(KSbeta_D2)
    GET_FIELD(KSbeta_D0)
    GET_FIELD(KSbeta_D1)
    GET_FIELD(KSbeta_D2)
    
    ADD_FIELD(KSalpha)
    GET_FIELD(KSalpha)
    
    /* beta and alpha needed */
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x   = patch->node[ijk]->x[0];
      double y   = patch->node[ijk]->x[1]-C_BH;
      double z   = patch->node[ijk]->x[2];
      double r2 = SQR(x)+SQR(y)+SQR(z);
      double rbar2  = 0.5*(r2-a2_BH+sqrt(SQR(r2-a2_BH)+4*a2_BH*SQR(z)));
      double rbar   = sqrt(rbar2);
      double k0 = (rbar*x+a_BH*y)/(rbar2+a2_BH);
      double k1 = (rbar*y-a_BH*x)/(rbar2+a2_BH);
      double k2 = z/rbar;
      double H  = M_BH*rbar/(rbar2+a2_BH*SQR(k2));
      double C = 2.*H;
      
      KSalpha[ijk] = 1/sqrt(1+C);
      KSbeta_D0[ijk]  = C*k0;
      KSbeta_D1[ijk]  = C*k1;
      KSbeta_D2[ijk]  = C*k2;
      
      /* note the followings are multiplied by _gammaI, 
      // they need also multiplication by (psi)^-4 to make gammaI 
      // which will be done after psi is made */
      double shift_U0 = 
KSbeta_D0[ijk]*_gammaI_U0U0[ijk] + KSbeta_D1[ijk]*_gammaI_U0U1[ijk] + 
KSbeta_D2[ijk]*_gammaI_U0U2[ijk];

      double shift_U1 = 
KSbeta_D0[ijk]*_gammaI_U0U1[ijk] + KSbeta_D1[ijk]*_gammaI_U1U1[ijk] + 
KSbeta_D2[ijk]*_gammaI_U1U2[ijk];

      double shift_U2 = 
KSbeta_D0[ijk]*_gammaI_U0U2[ijk] + KSbeta_D1[ijk]*_gammaI_U1U2[ijk] + 
KSbeta_D2[ijk]*_gammaI_U2U2[ijk];


      /* populating: */
      Beta_U1[ijk] = shift_U1;
      Beta_U0[ijk] = shift_U0;
      Beta_U2[ijk] = shift_U2;
      
    }
    
  }/* end of black hole part */
  
  /* initialization psi, eta and matter fields: */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    GET_FIELD(psi)
    GET_FIELD(eta)
    GET_FIELD(KSalpha)
    
    if (IsItNSPatch(patch))
    {
      Interpolation_T *interp_psi = init_interpolation();
      interp_psi->method          = "Natural_Cubic_Spline_1D";
      interp_psi->N_cubic_spline_1d->f = tov->psi;
      interp_psi->N_cubic_spline_1d->x = tov->rbar;
      interp_psi->N_cubic_spline_1d->N = tov->N;
      plan_interpolation(interp_psi);

      Interpolation_T *interp_h = init_interpolation();
      interp_h->method          = "Natural_Cubic_Spline_1D";
      interp_h->N_cubic_spline_1d->f = tov->h;
      interp_h->N_cubic_spline_1d->x = tov->rbar;
      interp_h->N_cubic_spline_1d->N = tov->N;
      plan_interpolation(interp_h);

      EoS_T *eos = initialize_EoS();
      
      GET_FIELD(enthalpy)
      GET_FIELD(rho0)
      GET_FIELD(phi)
      GET_FIELD(W_U0)
      GET_FIELD(W_U1)
      GET_FIELD(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        /* note that we naturally using isotropic coords. 
        // for our coordiante, so bar in rbar is dropped */
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1]-C_NS;
        double z = patch->node[ijk]->x[2];
        double r = sqrt(SQR(x)+SQR(y)+SQR(z));
        double alpha;
        double enthalpy_h;
        
        interp_psi->N_cubic_spline_1d->h = r;
        interp_h->N_cubic_spline_1d->h = r;
        
        /* psi */
        psi[ijk] = execute_interpolation(interp_psi);
        /* + 1 for KS but we won't add and so we won't subtrac 1 either */
        
        /* eta */
        enthalpy_h = execute_interpolation(interp_h);
        alpha = sqrt(1-0.5*M_NS/R_Schwar)/enthalpy_h/* NS part */ + 
                KSalpha[ijk]/* BH part */ - 1./* supper-position */;
        eta[ijk] = psi[ijk]*alpha;
        
        /* enthalpy */
        enthalpy[ijk] = enthalpy_h;
        
        /* rho0 */
        eos->h = enthalpy_h;
        rho0[ijk] = eos->rest_mass_density(eos);
        
        /* phi corrotating approximation */
        phi[ijk] = 0;
        
        /* spin part */
        W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
        W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
        W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
      }
      free_interpolation(interp_psi);
      free_interpolation(interp_h);
      free_EoS(eos);
    }
    else/* outside NS */
    {
      for (ijk = 0; ijk < nn; ++ijk)
      {
        /* note that we naturally using isotropic for our coordiante, 
        // so bar in rbar is dropped */
        double x    = patch->node[ijk]->x[0];
        double y    = patch->node[ijk]->x[1]-C_NS;
        double z    = patch->node[ijk]->x[2];
        double r = sqrt(SQR(x)+SQR(y)+SQR(z));
        double alpha;
        
        /* psi */
        psi[ijk] = 1+0.5*M_NS/r;
        /* + 1 for KS but we won't add and so we won't subtrac 1 either */
        
        /* eta */
        alpha = (1-0.5*M_NS/r)/(1+0.5*M_NS/r)+KSalpha[ijk]-1;
        eta[ijk] = psi[ijk]*alpha;
        
      }
    }
    
  }/* end of initialization psi, eta and matter fields */
  
  /* initializing Beta and B's */
  FOR_ALL_PATCHES(p,grid)
  {
     Patch_T *patch = grid->patch[p];
     unsigned nn = patch->nn;
     unsigned ijk;
     double psim4;/* psi^-4 */
      
     GET_FIELD(psi)
     GET_FIELD(B0_U0)
     GET_FIELD(B0_U1)
     GET_FIELD(B0_U2)
     GET_FIELD(B1_U0)
     GET_FIELD(B1_U1)
     GET_FIELD(B1_U2)
     GET_FIELD(Beta_U0)
     GET_FIELD(Beta_U1)
     GET_FIELD(Beta_U2)
    
     for (ijk = 0; ijk < nn; ++ijk)
     {
       double x   = patch->node[ijk]->x[0];
       double y   = patch->node[ijk]->x[1];
       
       psim4 = pow(psi[ijk],-4);
       
       /* Beta */
       Beta_U0[ijk] *= psim4;
       Beta_U1[ijk] *= psim4;
       Beta_U2[ijk] *= psim4;
       
       /* B1 */
       B1_U0[ijk] = Omega_BHNS*(-y+y_CM)+Vr*x/D;
       B1_U1[ijk] = Omega_BHNS*x+Vr*(y-y_CM)/D;
       B1_U2[ijk] = 0;
       
       /* B0 */
       B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
       B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
       B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
   }
      
  }/* end of * initializing Beta and B */
  
  /* freeing */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    DECLARE_FIELD(KSbeta_D0)
    DECLARE_FIELD(KSbeta_D1)
    DECLARE_FIELD(KSbeta_D2)
    DECLARE_FIELD(KSalpha)

    REMOVE_FIELD(KSbeta_D0)
    REMOVE_FIELD(KSbeta_D1)
    REMOVE_FIELD(KSbeta_D2)
    REMOVE_FIELD(KSalpha)
  }
  
  printf("Initializing the fields using TOV and Kerr-Schild solution ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}

/* given the radius of NS and BH and their separation,
// create a grid with these properties.
// ->return value: grid of NS and BH in which inside of the BH excised. */
static Grid_T *creat_grid_TOV_KerrShild(const double R_NS_l,const double R_BH_r,const double a_BH)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  /* calculate the characteristics of this grid */
  const unsigned gn = grid->gn;
  const double C    = GetParameterD_E("BH_NS_separation");
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double *R_outermost = alloc_double(N_Outermost_Split);
  double box_size_l;
  unsigned nlb[3]/*left box*/,n;
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  char val[100] = {'\0'};
  const char *kind;
  unsigned i;
  
  /* finding the kind of grid */
  kind = GetParameterS_E("grid_kind");
  grid->kind = dup_s(kind);
  
  assert(GRT(C,0));
  assert(GRT(R_NS_l,0));
  assert(GRT(R_BH_r,0));
  assert(LSS(2*R_NS_l,C));
  assert(LSS(2*R_BH_r,C));
  
  /* making NS and BH surfaces function */
  NS_BH_surface_CubedSpherical_grid(grid,R_NS_l,R_BH_r,a_BH);
  
  box_size_l = GetParameterD_E("left_central_box_length_ratio")*R_NS_l;
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = GetParameterD_E(var);
    
    if (LSS(R_outermost[i],2*C))
      abortEr("the radius of outermost patches must be greater than twice of BBN distance.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        abortEr("The radius of outermost must be increasing.");
    
  }
  
  /* filling n */
  
  /* left box */
  nlb[0] = (unsigned)GetParameterI("n_a");
  nlb[1] = (unsigned)GetParameterI("n_b");
  nlb[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  /* adding the results to the parameter data base */
  
  /* n_a, n_b, n_c */
  /* left box */
  sprintf(par,"grid%u_left_centeral_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_centeral_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_left_centeral_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  /* surrounding box length */
  sprintf(par,"grid%u_surrounding_box_length",gn);
  add_parameter_double(par,C);
  
  /* R1 and R2 outermost */
  sprintf(par,"grid%u_outermost%u_R2",gn,0);
  add_parameter_double(par,R_outermost[0]);
    
  for (i = 1; i < N_Outermost_Split; i++)
  {
    /* R1: */
    sprintf(par,"grid%u_outermost%u_R1",gn,i);
    add_parameter_double(par,R_outermost[i-1]);
    
    /* R2: */
    sprintf(par,"grid%u_outermost%u_R2",gn,i);
    add_parameter_double(par,R_outermost[i]);
    
  }
  
  /* assuming the center of left NS at (0,-C/2,0) */
  sprintf(par,"grid%u_left_NS_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_left_NS_center_b",gn);
  add_parameter_double(par,-C/2);
  
  sprintf(par,"grid%u_left_NS_center_c",gn);
  add_parameter_double(par,0.0);
  
  /* assuming the center of right BH at (0,C/2,0) */
  sprintf(par,"grid%u_right_BH_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_BH_center_b",gn);
  add_parameter_double(par,C/2);
  
  sprintf(par,"grid%u_right_BH_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R_outermost);

  
  make_patches(grid);/* making patch(es) to cover the grid */
  realize_geometry(grid);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  
  return grid;
  
}

/* making  NS and BH surfaces function */
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,const double R_NS_l,const double R_BH_r,const double a_BH)
{
  double *R;
  char par[100] = {'\0'};
  unsigned N[3],n,i,j,k,N_total;
  Patch_T patch[1] = {0};
  struct Collocation_s coll_s[2] = {0};
  double X[2],r;
  
  /* left NS */
  
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("left_NS_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("left_NS_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("left_NS_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
    
  N_total = N[0]*N[1]*N[2];
  
  /* surface */
  R = alloc_double(N_total);
  for (i = 0; i < N[0]; ++i)
    for (j = 0; j < N[1]; ++j)
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = R_NS_l;
      
  sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
  add_parameter_array(par,R,N_total);
  
  free(R);
  
  /* right BH */
  
  /* filling min */
  patch->min[0] = -1;
  patch->min[1] = -1;
  patch->min[2] = 0;

  /* filling max */
  patch->max[0] = 1;
  patch->max[1] = 1;
  patch->max[2] = 1;
  
  /* collocation */
  patch->collocation[0] = Chebyshev_Extrema;
  patch->collocation[1] = Chebyshev_Extrema;
  patch->collocation[2] = Chebyshev_Extrema;

  /* basis */
  patch->basis[0] = Chebyshev_Tn_BASIS;
  patch->basis[1] = Chebyshev_Tn_BASIS;
  patch->basis[2] = Chebyshev_Tn_BASIS;
    
  /* filling N */
  N[0] = (unsigned)GetParameterI("n_a");
  N[1] = (unsigned)GetParameterI("n_b");
  N[2] = (unsigned)GetParameterI("n_c");
  
  /* check for override */
  n = (unsigned)GetParameterI("right_BH_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("right_BH_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("right_BH_n_c");
  if (n != INT_MAX)     N[2] = n;
  
  if(N[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(N[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(N[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
   
  patch->n[0] = N[0];
  patch->n[1] = N[1];
  patch->n[2] = N[2];
  
  initialize_collocation_struct(patch,&coll_s[0],0);
  initialize_collocation_struct(patch,&coll_s[1],1);
  
  N_total = N[0]*N[1]*N[2];
  
  R = alloc_double(N_total);
  
  /* surface up and down */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);
      r = sqrt(
               (1+SQR(X[0])+SQR(X[1]))/
               ((SQR(X[0])+SQR(X[1]))/(SQR(R_BH_r)+SQR(a_BH)) + 1/SQR(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_up",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_right_BH_surface_function_down",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface back */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = z/x */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = y/x */
      r = sqrt(
               (1+SQR(X[0])+SQR(X[1]))/
               (((1+SQR(X[1])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[0])/SQR(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_back",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface front */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = y/x */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = z/x */
      r = sqrt(
               (1+SQR(X[0])+SQR(X[1]))/
               (((1+SQR(X[0])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[1])/SQR(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_front",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface left */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = x/y */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = z/y */
      r = sqrt(
               (1+SQR(X[0])+SQR(X[1]))/
               (((1+SQR(X[0])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[1])/SQR(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_left",grid->gn);
  add_parameter_array(par,R,N_total);
  
  /* surface right */
  for (i = 0; i < N[0]; ++i)
  {
    X[0] = point_value(i,&coll_s[0]);/* a = z/y */
    for (j = 0; j < N[1]; ++j)
    {
      X[1] = point_value(j,&coll_s[1]);/* b = x/y */
      r = sqrt(
               (1+SQR(X[0])+SQR(X[1]))/
               (((1+SQR(X[1])))/(SQR(R_BH_r)+SQR(a_BH)) + SQR(X[0])/SQR(R_BH_r))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_right_BH_surface_function_right",grid->gn);
  add_parameter_array(par,R,N_total);
  
  free(R);
}

