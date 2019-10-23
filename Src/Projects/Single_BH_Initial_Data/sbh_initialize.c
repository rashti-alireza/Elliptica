/*
// Alireza Rashti
// October 2019
*/

#include "sbh_initialize.h"

/* initialize this system according to the previous grid. 
// ->return value: the next grid as a result of this initialization. */
Grid_T *sbh_initialize_next_grid(Grid_T *const grid_prev)
{
  Grid_T *grid_next = 0;
  
  if (!grid_prev)/* if grid is empty come up with an approximation */
  {
    /* if we use TOV and Kerr-Schil black hole approximation */
    if (strcmp_i(GetParameterS_E("BH_initialization"),"KerrShild"))
      grid_next = KerrShild_approximation();
    else
      abortEr(NO_OPTION);
  }
  else/* use previous grid to make the next one with new adjustments */
  {
    grid_next = make_next_grid_using_previous_grid(grid_prev);
  }
  
  return grid_next;   
}

/* finding different quantities and then make the next grid using previous grid
// ->return value: the next grid called 'grid_next' */
static Grid_T *make_next_grid_using_previous_grid(Grid_T *const grid_prev)
{
  Grid_T *grid_next = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of next grid */
  
  /* calculate ADM momenta */
  calculate_P_ADMs(grid_prev);
  
  /* find the BH radius to acquire the desired BH mass */
  //find_BH_radius(grid_prev);
  
  /* find the Omega_BH to acquire the desired BH spin */
  //find_BH_Omega(grid_prev);
  
  /* make new grid with new parameters */
  const double bh_chi  = GetParameterD_E("BH_X_U2");
  const double bh_mass = GetParameterD_E("BH_mass");
  const double bh_R    = bh_mass*(1+sqrt(1-SQR(bh_chi)));
  GridParams->R_BH = bh_R;
  GridParams->a_BH   = bh_chi*bh_mass;
  grid_next = creat_sbh_grid_CS(GridParams);
  
  /* fields: */
  /* creating all of the fields needed for construction of Initial Data */
  sbh_allocate_fields(grid_next);
  
  /* populating the free data part of initial data that we chose ourself */
  sbh_populate_free_data(grid_next);

  /* use previous grid to interpolate values of the fields for 
  // the next grid  and initialzing some other fields */
  interpolate_and_initialize_to_next_grid(grid_next,grid_prev);
  
  /* taking partial derivatives of the fields needed for equations */
  sbh_partial_derivatives_fields(grid_next);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  sbh_update_Aij(grid_next);
  
  /* make normal vectorn on BH horizon */
  make_normal_vector_on_BH_horizon(grid_next);
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  
  return grid_next;
}

/* use previous grid to interpolate values of the fields that will be solved for the next grid */
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev)
{
  const unsigned np = grid_next->np;
  unsigned p;
 
  /* the following fields are interpolated: */
  /* B0_U[0-2],psi,eta */
  
  /* to avoid race condition between threads write all coeffs */
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid_prev->patch[p];
    
    DECLARE_FIELD(B0_U0)
    DECLARE_FIELD(B0_U1)
    DECLARE_FIELD(B0_U2)
    DECLARE_FIELD(psi)
    DECLARE_FIELD(eta)
    make_coeffs_3d(B0_U0);
    make_coeffs_3d(B0_U1);
    make_coeffs_3d(B0_U2);
    make_coeffs_3d(psi);
    make_coeffs_3d(eta);
    
  }
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid_next->patch[p];
    unsigned nn = patch->nn;
    char hint[100],*root_name;
    unsigned ijk;
    
    root_name = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
    assert(root_name);
    root_name++;
    sprintf(hint,"%s",root_name);
    
    PREP_FIELD(B0_U0)
    PREP_FIELD(B0_U1)
    PREP_FIELD(B0_U2)
    PREP_FIELD(psi)
    PREP_FIELD(eta)
  
    for (ijk = 0; ijk < nn; ++ijk)
    {
      double x[3] = { patch->node[ijk]->x[0],
                      patch->node[ijk]->x[1],
                      patch->node[ijk]->x[2] };
      double Xp[3] = {0};/* (X,Y,Z)(x,y,z) in grid_prev */
      Patch_T *patchp = 0;/* patch in grid_prev contains (x,y,z) */
      
      /* finding X and patch in grid_prev, associated to x */
      find_X_and_patch(x,hint,grid_prev,Xp,&patchp);
      
      B0_U0[ijk] = interpolate_from_patch_prim("B0_U0",Xp,patchp);
      B0_U1[ijk] = interpolate_from_patch_prim("B0_U1",Xp,patchp);
      B0_U2[ijk] = interpolate_from_patch_prim("B0_U2",Xp,patchp);
      psi[ijk]   = interpolate_from_patch_prim("psi",Xp,patchp);
      eta[ijk]   = interpolate_from_patch_prim("eta",Xp,patchp);
    }
  }/* end of for (p = 0; p < np; ++p) */
  
  /* initializing some other fields: */
  /* Beta_U[0-2],B1_U[0-2] */
  //const double Omega_BHNS = GetParameterD_E("BH_orbital_angular_velocity");
  //const double Vr   = GetParameterD_E("BH_infall_velocity");
  //const double y_CM = GetParameterD_E("y_CM");  
  //const double D    = GetParameterD_E("BH_NS_separation");
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid_next->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    GET_FIELD(B0_U0)
    GET_FIELD(B0_U1)
    GET_FIELD(B0_U2)
    PREP_FIELD(B1_U0)
    PREP_FIELD(B1_U1)
    PREP_FIELD(B1_U2)
    PREP_FIELD(Beta_U0)
    PREP_FIELD(Beta_U1)
    PREP_FIELD(Beta_U2)
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      //double x = patch->node[ijk]->x[0];
      //double y = patch->node[ijk]->x[1];
       
      /* B1 */
      B1_U0[ijk] = 0;//Omega_BHNS*(-y+y_CM)+Vr*x/D;
      B1_U1[ijk] = 0;//Omega_BHNS*x+Vr*(y-y_CM)/D;
      B1_U2[ijk] = 0;
       
      /* Beta */
      Beta_U0[ijk] = B0_U0[ijk]+B1_U0[ijk];
      Beta_U1[ijk] = B0_U1[ijk]+B1_U1[ijk];
      Beta_U2[ijk] = B0_U2[ijk]+B1_U2[ijk];
    }
    
  }/* end of for (p = 0; p < np; ++p) */
  
}

/* given field name, X and patch, finds the value of the field in X  
// using interpolation.
// ->return value: f(X) */
static double interpolate_from_patch_prim(const char *const field,const double *const X,Patch_T *const patch)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  
  interp_s->field = patch->pool[Ind(field)];
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* given a cartesian point x on the grid, it finds the corresponding X and patch 
// on which this x takes place. 
// hint, is the name of the patch that potentially has the given x */
static void find_X_and_patch(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch)
{
  Interface_T **face;
  SubFace_T *subf;
  Needle_T *needle = alloc_needle();
  unsigned *found;
  unsigned p,f,sf;
  
  needle->grid = grid;
  needle->x    = x;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!strstr(patch->name,hint))
    {
      continue;
    }
    else
    {
      needle_in(needle,patch);
      
      /* find all neighbors of this patch */
      face = patch->interface;
      /* for all faces */
      FOR_ALL(f,face)
      {
        /* for all subfaces */
        for (sf = 0; sf < face[f]->ns; ++sf)
        {
          subf = face[f]->subface[sf];
          if (subf->outerB || subf->innerB)
            continue;
            
          needle_in(needle,grid->patch[subf->adjPatch]);
        }
      }
      break;
    }/* end of else */
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  point_finder(needle);
  found = needle->ans;
  
  /* if it could not find X in neither the given hint patch nor its neighbors */
  if (!needle->Nans)
  {
    *ppatch = grid->patch[found[0]];
    X_of_x(X,x,*ppatch);
  }
  else/* if no patch found let's find it in the other patches */
  {
    needle->ex   = needle->in;
    needle->Nex  = needle->Nin;
    needle->in   = 0;
    needle->Nin  = 0;
    
    point_finder(needle);
    found = needle->ans;
    
    if (!needle->Nans)
      abortEr("It could not find the x at the given patches!\n");
      
    *ppatch = grid->patch[found[0]];
    X_of_x(X,x,*ppatch);
  }
  free_needle(needle);
}

/* Kerr-Schild black hole approximation.
// ->return value: resultant grid from this approximation */
static Grid_T *KerrShild_approximation(void)
{
  Grid_T *grid = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of grid */
  
  /* basics of Kerr Shild black hole located at right side of y axis */
  pr_line_custom('=');
  printf("Acquiring Black Hole properties ...\n");
  const double bh_chi  = GetParameterD_E("BH_X_U2");
  const double bh_mass = GetParameterD_E("BH_mass");
  const double bh_R    = bh_mass*(1+sqrt(1-SQR(bh_chi)));
  printf("BH properties:\n");
  printf("--> BH radius (Kerr-Schild Coords.) = %e\n",bh_R);
  printf("--> BH dimensionless spin (z comp.) = %e\n",bh_chi);
  printf("--> BH ADM mass                     = %e\n",bh_mass);
  printf("Acquiring Black Hole properties ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
 
  /* combining these two geometry to create the grid */
  GridParams->R_BH = bh_R;
  GridParams->a_BH = bh_chi*bh_mass;
  grid = creat_sbh_grid_CS(GridParams);
  
  /* creating all of the fields needed for construction of Initial Data */
  sbh_allocate_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  sbh_populate_free_data(grid);
  
  /* initialize the fields using TOV and Kerr-Shild solution */
  init_field_KerrSchild(grid,bh_chi*bh_mass,bh_mass);
  
  /* taking partial derivatives of the fields needed for equations */
  sbh_partial_derivatives_fields(grid);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  sbh_update_Aij(grid);
  
  /* make normal vectorn on BH horizon */
  make_normal_vector_on_BH_horizon(grid);
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  
  return grid;
}

/* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
// _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
static void sbh_update_Aij(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ...\n");
  unsigned p;

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    sbh_update_psi10A_UiUj(patch);
  }
  
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* make normal vectorn on BH horizon */
static void make_normal_vector_on_BH_horizon(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Making normal vector on BH horizon ...\n");
  
  unsigned p,nn,ijk;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItHorizonPatch(patch))
      continue;
      
    nn = patch->nn;
    
    GET_FIELD(_gamma_D2D2)
    GET_FIELD(_gamma_D0D2)
    GET_FIELD(_gamma_D0D0)
    GET_FIELD(_gamma_D0D1)
    GET_FIELD(_gamma_D1D2)
    GET_FIELD(_gamma_D1D1)
    
    /* normal vector on horizon */
    PREP_FIELD(_HS_U0);
    PREP_FIELD(_HS_U1);
    PREP_FIELD(_HS_U2);
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      /* minus sign to point outside the black hole */
      _HS_U0[ijk] = dq2_dq1(patch,_c_,_x_,ijk);
      _HS_U1[ijk] = dq2_dq1(patch,_c_,_y_,ijk);
      _HS_U2[ijk] = dq2_dq1(patch,_c_,_z_,ijk);
      
      double N2 = 
pow(_HS_U0[ijk], 2)*_gamma_D0D0[ijk] + 2.0*_HS_U0[ijk]*_HS_U1[ijk]*
_gamma_D0D1[ijk] + 2.0*_HS_U0[ijk]*_HS_U2[ijk]*_gamma_D0D2[ijk] +
pow(_HS_U1[ijk], 2)*_gamma_D1D1[ijk] + 2.0*_HS_U1[ijk]*_HS_U2[ijk]*
_gamma_D1D2[ijk] + pow(_HS_U2[ijk], 2)*_gamma_D2D2[ijk];
        
      double N = sqrt(N2);
      
      /* normalizing */
      _HS_U0[ijk] /= N;
      _HS_U1[ijk] /= N;
      _HS_U2[ijk] /= N;
      
    }
     
  }
  
  printf("Making normal vector on BH horizon ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
  
}

/* initialize the fields using TOV and Kerr-Shild solution.
// the idea is to superimpose two fields of each solutions. */
static void init_field_KerrSchild(Grid_T *const grid,const double a_BH, const double M_BH)
{
  pr_line_custom('=');
  printf("Initializing the fields using Kerr-Schild solution ...\n");

  const double C_BH  = 0;
  const double a2_BH = SQR(a_BH);/* spin vector of BH */
  //const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  //const double Vr = GetParameterD_E("BH_NS_infall_velocity");
  unsigned p;
  
  //add_parameter_double("y_CM",y_CM);
  
  /* black hole parts */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    PREP_FIELD(Beta_U0)
    PREP_FIELD(Beta_U1)
    PREP_FIELD(Beta_U2)
    PREP_FIELD(_gammaI_U0U2)
    PREP_FIELD(_gammaI_U0U0)
    PREP_FIELD(_gammaI_U0U1)
    PREP_FIELD(_gammaI_U1U2)
    PREP_FIELD(_gammaI_U1U1)
    PREP_FIELD(_gammaI_U2U2)

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
  
  /* initialization psi, eta: */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    PREP_FIELD(psi)
    PREP_FIELD(eta)
    PREP_FIELD(KSalpha)
    
    for (ijk = 0; ijk < nn; ++ijk)
    {
      /* note that we naturally using isotropic for our coordiante, 
      // so bar in rbar is dropped */
      double alpha;
      
      /* psi */
      psi[ijk] = 1;
      /* + 1 for KS but we won't add and so we won't subtrac 1 either */
      
      /* eta */
      alpha = KSalpha[ijk];
      eta[ijk] = psi[ijk]*alpha;
      
    }
    
  }/* end of initialization psi, eta and matter fields */
  
  /* initializing Beta and B's */
  FOR_ALL_PATCHES(p,grid)
  {
     Patch_T *patch = grid->patch[p];
     unsigned nn = patch->nn;
     unsigned ijk;
     double psim4;/* psi^-4 */
      
     PREP_FIELD(psi)
     PREP_FIELD(B0_U0)
     PREP_FIELD(B0_U1)
     PREP_FIELD(B0_U2)
     PREP_FIELD(B1_U0)
     PREP_FIELD(B1_U1)
     PREP_FIELD(B1_U2)
     PREP_FIELD(Beta_U0)
     PREP_FIELD(Beta_U1)
     PREP_FIELD(Beta_U2)
    
     for (ijk = 0; ijk < nn; ++ijk)
     {
       //double x   = patch->node[ijk]->x[0];
       //double y   = patch->node[ijk]->x[1];
       
       psim4 = pow(psi[ijk],-4);
       
       /* Beta */
       Beta_U0[ijk] *= psim4;
       Beta_U1[ijk] *= psim4;
       Beta_U2[ijk] *= psim4;
       
       /* B1 */
       B1_U0[ijk] = 0;//Omega_BHNS*(-y+y_CM)+Vr*x/D;
       B1_U1[ijk] = 0;//Omega_BHNS*x+Vr*(y-y_CM)/D;
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
  
  printf("Initializing the fields using Kerr-Schild solution ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}

/* given the of radius and spin of BH ,
// as the parameters, create a grid with these properties.
// NOTE: WE assume we are using cubed spherical grid.
// ->return value: grid of single BH in which inside of the BH excised. */
static Grid_T *creat_sbh_grid_CS(struct Grid_Params_S *const GridParams)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  /* calculate the characteristics of this grid */
  const double R_BH = GridParams->R_BH;
  const unsigned gn = grid->gn;
  const unsigned N_Outermost_Split = (unsigned)GetParameterI_E("Number_of_Outermost_Split"); 
  double *R_outermost = alloc_double(N_Outermost_Split);
  char var[100] = {'\0'};
  char par[100] = {'\0'};
  const char *kind;
  unsigned i;
  
  /* finding the kind of grid */
  kind = GetParameterS_E("grid_kind");
  if (!strcmp_i(kind,"SBH_CubedSpherical_grid"))
    abortEr("This function only works with cubed spherical grid.\n");
    
  grid->kind = dup_s(kind);
  
  assert(GRT(R_BH,0));
  
  /* making BH surfaces function */
  BH_surface_CubedSpherical_grid(grid,GridParams);
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = GetParameterD_E(var);
    
    if (LSS(R_outermost[i],4*R_BH))
      abortEr("the radius of outermost patches must be greater than fourth times of of BH radius.");
    
    if (i > 0)
      if (LSSEQL(R_outermost[i],R_outermost[i-1]))
        abortEr("The radius of outermost must be increasing.");
    
  }
  
  /* adding the results to the parameter data base */
  
  /* surrounding box length */
  sprintf(par,"grid%u_surrounding_box_length",gn);
  add_parameter_double(par,2*R_BH);

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
  
  /* assuming the center of BH at (0,0,0) */
  sprintf(par,"grid%u_BH_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_BH_center_b",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_BH_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R_outermost);
  
  make_patches(grid);/* making patch(es) to cover the grid */
  realize_geometry(grid);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  
  return grid;
}

/* making BH surfaces function */
static void BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  const double R_BH = GridParams->R_BH;
  const double a_BH = GridParams->a_BH;
  double *R;
  char par[1000] = {'\0'};
  unsigned N[3],n,i,j,k,N_total;
  Patch_T patch[1] = {0};
  struct Collocation_s coll_s[2] = {0};
  double X[3],r;
  
  /* right BH: */
  
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
  n = (unsigned)GetParameterI("BH_n_a");
  if (n != INT_MAX)     N[0] = n;
  n = (unsigned)GetParameterI("BH_n_b");
  if (n != INT_MAX)     N[1] = n;
  n = (unsigned)GetParameterI("BH_n_c");
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
               ((SQR(X[0])+SQR(X[1]))/(SQR(R_BH)+SQR(a_BH)) + 1/SQR(R_BH))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_BH_surface_function_up",grid->gn);
  add_parameter_array(par,R,N_total);
  sprintf(par,"grid%u_BH_surface_function_down",grid->gn);
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
               (((1+SQR(X[1])))/(SQR(R_BH)+SQR(a_BH)) + SQR(X[0])/SQR(R_BH))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_BH_surface_function_back",grid->gn);
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
               (((1+SQR(X[0])))/(SQR(R_BH)+SQR(a_BH)) + SQR(X[1])/SQR(R_BH))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_BH_surface_function_front",grid->gn);
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
               (((1+SQR(X[0])))/(SQR(R_BH)+SQR(a_BH)) + SQR(X[1])/SQR(R_BH))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_BH_surface_function_left",grid->gn);
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
               (((1+SQR(X[1])))/(SQR(R_BH)+SQR(a_BH)) + SQR(X[0])/SQR(R_BH))
              );
      for (k = 0; k < N[2]; ++k)
        R[L(N,i,j,k)] = r;
    }
  }
  sprintf(par,"grid%u_BH_surface_function_right",grid->gn);
  add_parameter_array(par,R,N_total);
  
  free(R);
}

/* initialize Grid_Params struct */
static struct Grid_Params_S *init_GridParams(void)
{
  struct Grid_Params_S *par = calloc(1,sizeof(*par));
  return par;
}

/* free Grid_Params struct */
static void free_Grid_Params_S(struct Grid_Params_S *par)
{
  //_free(par->NS_R_Ylm->realClm);
  //_free(par->NS_R_Ylm->imagClm);
  _free(par);
}

/* calculating ADM momenta */
static void calculate_P_ADMs(Grid_T *const grid)
{
  Observable_T *obs   = init_observable(grid);
  
  obs->quantity = "ADM_momentums";
  plan_observable(obs);
  
  const double Px_ADM = obs->Px_ADM(obs);
  const double Py_ADM = obs->Py_ADM(obs);
  const double Pz_ADM = obs->Pz_ADM(obs);
  
  printf("ADM momenta:\n"
         "(Px,Py,Pz) = (%e,%e,%e)\n",Px_ADM,Py_ADM,Pz_ADM);
  
  free_observable(obs);
}

