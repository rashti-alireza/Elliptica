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
  abortEr(NO_JOB);
  Grid_T *grid_next = 0;
  struct Grid_Params_S GridParams[1] = {0};/* adjust some pars for construction of next grid */
  
  /* find Euler equation constant to meet NS baryonic mass */
  find_Euler_eq_const(grid_prev);
  /* find BH center to make P_ADM zero */
  /* find the BH radius to acquire the desired BH mass */
  /* find the Omega_BH to acquire the desired BH spin */
  /* find y_CM using force balance equation */
  
  if (!strcmp_i(grid_prev->kind,"BBN_CubedSpherical_grid"))
  {
    /* extrapolate fluid fields outside of NS in case their value needed. */
    extrapolate_fluid_fields_outsideNS_CS(grid_prev);
    
    GridParams->NS_R_type = GetParameterS_E("NS_surface_collocation");
    
    /* find NS surface using cubed spherical points */
    if (strstr_i(GridParams->NS_R_type,"CubedSpherical"))
      find_NS_surface_CS_method_CS(grid_prev,GridParams);
      
    /* find NS surface using spherical harmonic points */
    else if (strstr_i(GridParams->NS_R_type,"SphericalHarmonic"))
      find_NS_surface_Ylm_method_CS(grid_prev,GridParams);
    else
      abortEr(NO_OPTION);
  }
  else
    abortEr(NO_OPTION);
    
  /* make new grid with new parameters */
  const double bh_chi  = GetParameterD_E("BH_X_U2");
  const double bh_mass = GetParameterD_E("BH_mass");
  const double bh_R    = bh_mass*(1+sqrt(1-SQR(bh_chi)));
  GridParams->R_BH_r = bh_R;
  GridParams->a_BH   = bh_chi*bh_mass;
  grid_next = creat_grid_CS(GridParams);
  
  /* fields: */
  /* creating all of the fields needed for construction of Initial Data */
  bbn_allocate_fields(grid_next);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid_next);

  /* use previous grid to interpolate values of the fields for 
  // the next grid  and initialzing some other fields */
  interpolate_and_initialize_to_next_grid(grid_next,grid_prev);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid_next);
  
  /* update u0, _J^i, _E and _S */
  Tij_IF_CTS_psi6Sources(grid_next);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid_next);
  
  /* make normal vectorn on BH horizon */
  make_normal_vector_on_BH_horizon(grid_next);
  
  return grid_next;
}

/* find Euler equation constant to meet NS baryonic mass */
static void find_Euler_eq_const(Grid_T *const grid)
{
  Root_Finder_T *root = init_root_finder(1);
  double *Euler_const = 0;
  double guess[1];/* initial guess for Euler const */
  struct Params_S
  {
    Grid_T *grid;
    double NS_baryonic_mass;
  }params[1];
  
  params->grid = grid;
  params->NS_baryonic_mass = GetParameterD_E("NS_baryonic_mass");
  guess[0] = GetParameterD_E("Euler_equation_constant");
  
  root->description = "Finding Euler equation constant using NS baryonic mass:";
  root->type        = GetParameterS_E("RootFinder_Method");
  root->tolerance   = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_rootfinder_eq;
  plan_root_finder(root);
  Euler_const       = execute_root_finder(root);
  
  update_parameter_double_format("Euler_equation_constant",Euler_const[0]);
  
  free(Euler_const);
  free_root_finder(root);
}

/* root finder eqution for Euler equation constant */
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x)
{
  struct Params_S
  {
    Grid_T *grid;
    double NS_baryonic_mass;
  }*par = params;
  
  return bbn_NS_baryonic_mass(par->grid,x[0]) - par->NS_baryonic_mass;
}

/* use previous grid to interpolate values of the fields that will be solved for the next grid */
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev)
{
  const unsigned np = grid_next->np;
  unsigned p;
 
  /* test for race condition, test for residual after interpolation in grid_next */
  printf("THIS function has not been tested deeply %s,%d\n",__FILE__,__LINE__);
  
  /* the following fields are interpolated: */
  /* B0_U[0-2],psi,eta,phi,enthalpy */
  
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
    
    if (IsItNSPatch(patch))
    {
      DECLARE_FIELD(phi)
      DECLARE_FIELD(enthalpy)
      
      make_coeffs_3d(B0_U0);
      make_coeffs_3d(B0_U1);
      make_coeffs_3d(B0_U2);
      make_coeffs_3d(psi);
      make_coeffs_3d(eta);
      make_coeffs_3d(phi);
      make_coeffs_3d(enthalpy);
    }
    else
    {
      make_coeffs_3d(B0_U0);
      make_coeffs_3d(B0_U1);
      make_coeffs_3d(B0_U2);
      make_coeffs_3d(psi);
      make_coeffs_3d(eta);
    }
    
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
    sprintf(hint,"%s",root_name);
    
    PREP_FIELD(B0_U0)
    PREP_FIELD(B0_U1)
    PREP_FIELD(B0_U2)
    PREP_FIELD(psi)
    PREP_FIELD(eta)
    
    if (IsItNSPatch(patch))
    {
      PREP_FIELD(phi)
      PREP_FIELD(enthalpy)
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x[3] = { patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2] };
        double Xp[3] = {0};/* (X,Y,Z)(x,y,z) in grid_prev */
        Patch_T *patchp = 0;/* patch in grid_prev contains (x,y,z) */
        
        /* finding X and patch in grid_prev, associated to x */
        find_Xp_and_patchp(x,hint,grid_prev,Xp,&patchp);
        
        B0_U0[ijk]    = interpolate_from_prev_grid("B0_U0",Xp,patchp);
        B0_U1[ijk]    = interpolate_from_prev_grid("B0_U1",Xp,patchp);
        B0_U2[ijk]    = interpolate_from_prev_grid("B0_U2",Xp,patchp);
        psi[ijk]      = interpolate_from_prev_grid("psi",Xp,patchp);
        eta[ijk]      = interpolate_from_prev_grid("eta",Xp,patchp);
        phi[ijk]      = interpolate_from_prev_grid("phi",Xp,patchp);
        enthalpy[ijk] = interpolate_from_prev_grid("enthalpy",Xp,patchp);
      }
    }
    else
    {
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x[3] = { patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2] };
        double Xp[3] = {0};/* (X,Y,Z)(x,y,z) in grid_prev */
        Patch_T *patchp = 0;/* patch in grid_prev contains (x,y,z) */
        
        /* finding X and patch in grid_prev, associated to x */
        find_Xp_and_patchp(x,hint,grid_prev,Xp,&patchp);
        
        B0_U0[ijk] = interpolate_from_prev_grid("B0_U0",Xp,patchp);
        B0_U1[ijk] = interpolate_from_prev_grid("B0_U1",Xp,patchp);
        B0_U2[ijk] = interpolate_from_prev_grid("B0_U2",Xp,patchp);
        psi[ijk]   = interpolate_from_prev_grid("psi",Xp,patchp);
        eta[ijk]   = interpolate_from_prev_grid("eta",Xp,patchp);
      }
    }
  }/* end of FOR_ALL_PATCHES(p,grid_next) */
  
  /* initializing some other fields: */
  /* rho0,W_U[0-2],Beta_U[0-2],B1_U[0-2] */
  const double Omega_BHNS = GetParameterD_E("BH_NS_orbital_angular_velocity");
  const double Omega_NS_x = GetParameterD_E("NS_Omega_U0");
  const double Omega_NS_y = GetParameterD_E("NS_Omega_U1");
  const double Omega_NS_z = GetParameterD_E("NS_Omega_U2");
  const double Vr   = GetParameterD_E("BH_NS_infall_velocity");
  const double y_CM = GetParameterDoubleF_E("y_CM");  
  const double D    = GetParameterD_E("BH_NS_separation");
  const double C_NS = GetParameterDoubleF_E("NS_Center");
  
  FOR_ALL_PATCHES(p,grid_next)
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
      double x     = patch->node[ijk]->x[0];
      double y     = patch->node[ijk]->x[1];
       
      /* B1 */
      B1_U0[ijk] = Omega_BHNS*(-y+y_CM)+Vr*x/D;
      B1_U1[ijk] = Omega_BHNS*x+Vr*(y-y_CM)/D;
      B1_U2[ijk] = 0;
       
      /* Beta */
      Beta_U0[ijk] = B0_U0[ijk]+B1_U0[ijk];
      Beta_U1[ijk] = B0_U1[ijk]+B1_U1[ijk];
      Beta_U2[ijk] = B0_U2[ijk]+B1_U2[ijk];
    }
    
    if (IsItNSPatch(patch))
    {
      EoS_T *eos = initialize_EoS();
      
      GET_FIELD(enthalpy)
      PREP_FIELD(rho0)
      PREP_FIELD(W_U0)
      PREP_FIELD(W_U1)
      PREP_FIELD(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1]-C_NS;
        double z = patch->node[ijk]->x[2];
        
        /* rho0 */
        eos->h    = enthalpy[ijk];
        rho0[ijk] = eos->rest_mass_density(eos);
        
        /* spin part */
        W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
        W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
        W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
      }
      free_EoS(eos);
    }/* end of if (IsItNSPatch(patch)) */
    
  }/* end of FOR_ALL_PATCHES(p,grid_next) */
  
}

/* given field name, X and patch, finds the value of the field in X  
// using interpolation.
// ->return value: f(X) */
static double interpolate_from_prev_grid(const char *const field,const double *const X,Patch_T *const patch)
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
static void find_Xp_and_patchp(const double *const x,const char *const hint,Grid_T *const grid,double *const X,Patch_T **const ppatch)
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
  
  /* if it could not find X in neither the given hint patch nor its neighbors, raise a flag! */
  if (!needle->Nans)
    abortEr("It could not find the x at the given patches!\n");
  
  *ppatch = grid->patch[found[0]];
  
  X_of_x(X,x,*ppatch);
  
  free_needle(needle);
}

/* given the grid, find the NS surface on all cubed spherical points 
// using the fact that at the surface enthalpy = 1.
// it fills also NS radius attributes:
//   GridParams->Max_R_NS_l;
//
// we assumed cubed spherical grid */
static void find_NS_surface_CS_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  pr_line_custom('=');
  printf("Finding the surface of NS, cubed spherical ...\n");
  
  UNUSED(GridParams);
  abortEr("setup R_max and etc. in GridParams");
  
  /* the stucture for the root finder */
  struct Params_S
  {
    Patch_T *patch;
    double x0[3];/* (x,y,z) at the surface */
    double *N;/* the direction of increasing or decreasing of x = x0+N*d */
    double Euler_C;/* Euler equation const. */
  } par[1];
  char par_str[1000];
  unsigned p;
  
  par->Euler_C = GetParameterD_E("Euler_equation_constant");
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    /* finding at which patch h equals 1 takes place */
    if (!IsItNSPatch(patch))
      continue;
    if (strstr(patch->name,"left_centeral_box"))
      continue;
    
    const double *c = patch->c;
    const unsigned *n = patch->n;
    double *R_NS = alloc_double(patch->nn);
    double R0_NS;/* initial length of R_NS */
    double Rnew_NS;/* new R_NS */
    double *x_sol;
    double x,y,z;
    unsigned ijk,i,j,k,kp;
      
    GET_FIELD(enthalpy)
    
    /* go over NS surface and test the value of h */
    k = n[2]-1;/* this is where the NS surface takes place */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        Patch_T *h_patch = 0;
        Flag_T NS_patch_flg = NONE;
        
        ijk   = L(n,i,j,k);
        x     = patch->node[ijk]->x[0]-c[0];
        y     = patch->node[ijk]->x[1]-c[1];
        z     = patch->node[ijk]->x[2]-c[2];
        R0_NS = sqrt(SQR(x)+SQR(y)+SQR(z));
        
        if(EQL(enthalpy[ijk],1))/* if it takes place on the current NS surface */
        {
          for (kp = 0; kp < n[2]; ++kp)
            R_NS[L(n,i,j,kp)] = R0_NS;
        }
        else/* if NS surface is displaced */
        {
          Point_T point[1];
          double *N;/* normal vector */
            
          point->ind = ijk;
          point->patch = patch;
          point->face  = K_n2;
          
          if (LSS(enthalpy[ijk],1))
          {
            h_patch      = patch;
            NS_patch_flg = YES;
          }
          else/* which means h = 1 occures in neighboring patch */
          {
            NS_patch_flg = NO;
            /* find the neighbor of this patch at face K_n2 
            // using the maximum number of common points. */
            Interface_T *face = patch->interface[K_n2];
            unsigned s,max = 0;
            
            for (s = 0; s < face->ns; ++s)
            {
              SubFace_T *subf = face->subface[s];
              /* since this is a patch in the middle of the grid, we won't
              // expect to have outerB, moreover, since this is NS patch
              // we should not have innerB. */
              assert(!subf->outerB);
              assert(!subf->innerB);
              
              if (subf->np > max)
              {
                max = subf->np;
                h_patch = grid->patch[subf->adjPatch];
              }
            }
            assert(max);/* max = 0 is wrong */
          }/* end of else */
          
          /* having found h_patch, now find the the position of h = 1 */
          N          = normal_vec(point);/* normal vector pointing outwards */
          par->patch = h_patch;
          par->x0[0] = patch->node[ijk]->x[0];
          par->x0[1] = patch->node[ijk]->x[1];
          par->x0[2] = patch->node[ijk]->x[2];
          par->N     = N;
          if (NS_patch_flg == YES)
          {
            /* normal vector is toward the NS center */
            N[0] *= -1;
            N[1] *= -1;
            N[2] *= -1;
          }
          
          /* populate root finder */
          Root_Finder_T *root = init_root_finder(1);
          root->type      = GetParameterS_E("RootFinder_Method");
          plan_root_finder(root);
          root->tolerance = GetParameterD_E("RootFinder_Tolerance");
          root->MaxIter   = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
          root->x_gss     = &R0_NS;
          root->params    = par;
          root->f[0]      = bbn_NS_surface_enthalpy_eq;
          
          if (NS_patch_flg == YES)
          {
            root->FD_Left  = 1;
          }
          else
            root->FD_Right = 1;
          
          double xp[3];
          x_sol   = execute_root_finder(root);
          /*  new coords of R respect to the center of NS */
          xp[0] = x+x_sol[0]*N[0];
          xp[1] = y+x_sol[0]*N[1];
          xp[2] = z+x_sol[0]*N[2];
          
          Rnew_NS = rms(3,xp,0);
          free_root_finder(root);
          free(x_sol);
          
          /* having found new R now fill R_NS field */
          for (kp = 0; kp < n[2]; ++kp)
            R_NS[L(n,i,j,kp)] = Rnew_NS;
        }
      }/* end of for (j = 0; j < n[1]; ++j) */
    }/* end of for (i = 0; i < n[0]; ++i) */
    
    sprintf(par_str,"grid%u_%s_h=1",grid->gn+1/* since we want it for the next grid */
                                   ,patch->name);
    add_parameter_array(par_str,R_NS,patch->nn);
    free(R_NS);
  }/* end of FOR_ALL_PATCHES(p,grid) */
  
  printf("Finding the surface of NS, cubed spherical method ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

#define ij(i,j) (j)+Nphi*(i)
/* given the grid, find the NS surface on Ylm points (Legendre,EquiSpaced)
// using the fact that at the surface enthalpy = 1.
// it fills also NS radius attributes:
//   GridParams->Max_R_NS_l;
//   GridParams->NS_R_Ylm->realClm;
//   GridParams->NS_R_Ylm->imagClm;
//   GridParams->NS_R_Ylm->Lmax;
//
// we assumed cubed spherical grid */
static void find_NS_surface_Ylm_method_CS(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  pr_line_custom('=');
  printf("Finding the surface of NS, Ylm method ...\n");
  
  /* the stucture for the root finder */
  struct Params_S
  {
    Patch_T *patch;
    double x0[3];/* (x,y,z) at the surface */
    double *N;/* the direction of increasing or decreasing of x = x0+N*d */
    double Euler_C;/* Euler equation const. */
  } par[1];
  unsigned Ntheta,Nphi;/* total number of theta and phi points */
  const unsigned lmax = (unsigned)GetParameterI_E("NS_surface_Ylm_expansion_max_l");
  double theta,phi;
  double *Rnew_NS = 0;/* new R for NS */
  double Max_R_NS = 0;/* maximum radius of NS */
  double X[3],x[3],N[3];
  char stem[1000],*affix;
  Flag_T NS_patch_flg = NONE;
  unsigned i,j;
  
  par->Euler_C = GetParameterD_E("Euler_equation_constant");
  
  /* initialize tables */
  init_Legendre_root_function();
  
  Ntheta  = Nphi = 2*lmax+1;
  Rnew_NS = alloc_double(Ntheta*Nphi);
  
  /* for each points of Ylm find the surface of NS */
  for (i = 0; i < Ntheta; ++i)
  {
    theta = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      phi = j*2*M_PI/Nphi;
      NS_patch_flg = NONE;
      Patch_T *h_patch = 0,*patch = 0;
      double h,y[3] = {0},R0_NS,*dr;
      
      /* find patch and X,Y,Z at NS surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_NS_CS(X,&patch,theta,phi,grid);
      
      /* find enthalpy at the (X,Y,Z) */
      Interpolation_T *interp_h = init_interpolation();
      interp_h->field = patch->pool[Ind("enthalpy")];
      interp_h->XYZ_dir_flag = 1;
      interp_h->X            = X[0];
      interp_h->Y            = X[1];
      interp_h->Z            = X[2];
      plan_interpolation(interp_h);
      h = execute_interpolation(interp_h);/* enthalpy */
      free_interpolation(interp_h);
      
      /* finding x */
      x_of_X(x,X,patch);
      y[0] = x[0]-patch->c[0];
      y[1] = x[1]-patch->c[2];
      y[2] = x[2]-patch->c[1];
      R0_NS = rms(3,y,0);
      
      if(EQL(h,1))/* if it takes place on the current NS surface */
      {
        Rnew_NS[ij(i,j)] = R0_NS;
      }
      else/* if NS surface is displaced */
      {
        /* r^ = sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^ */
        N[0] = sin(theta)*cos(phi);
        N[1] = sin(theta)*sin(phi);
        N[2] = cos(theta);
        
        if (LSS(h,1))
        {
          h_patch      = patch;
          NS_patch_flg = YES;
        }
        else/* which means h = 1 occures in neighboring patch */
        {
          NS_patch_flg = NO;
          /* finding the side of the patch */
          affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);
          assert(affix);
          sprintf(stem,"left_NS_surrounding_%s",affix);
          free(affix);
          h_patch = GetPatch(stem,grid);
        }
        /* having found h_patch, now find the the position of h = 1 */
        par->patch = h_patch;
        par->x0[0] = x[0];
        par->x0[1] = x[1];
        par->x0[2] = x[2];
        par->N     = N;
        if (NS_patch_flg == YES)
        {
          /* normal vector is toward the NS center */
          N[0] *= -1;
          N[1] *= -1;
          N[2] *= -1;
        }
        
        /* populate root finder */
        Root_Finder_T *root = init_root_finder(1);
        root->type      = GetParameterS_E("RootFinder_Method");
        plan_root_finder(root);
        root->tolerance = GetParameterD_E("RootFinder_Tolerance");
        root->MaxIter   = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
        root->x_gss     = &R0_NS;
        root->params    = par;
        root->f[0]      = bbn_NS_surface_enthalpy_eq;
        
        if (NS_patch_flg == YES)
        {
          root->FD_Left  = 1;
        }
        else
          root->FD_Right = 1;
        
        dr    = execute_root_finder(root);
        /*  new coords of R respect to the center of NS */
        y[0] += N[0]*dr[0];
        y[1] += N[1]*dr[0];
        y[2] += N[2]*dr[0];
        
        Rnew_NS[ij(i,j)] = rms(3,y,0);
        if (Rnew_NS[ij(i,j)] > Max_R_NS)
          Max_R_NS = Rnew_NS[ij(i,j)];
          
        free_root_finder(root);
        free(dr);
        
      }/* end of else */
      
    }
  }
  
  /* adding maximum radius of NS to grid parameters */
  GridParams->Max_R_NS_l = Max_R_NS;
  
  /* making radius of NS parameter at each patch using Ylm interpolation */
  double *realClm = alloc_ClmYlm(lmax);
  double *imagClm = alloc_ClmYlm(lmax);
  
  init_Ylm();
  
  /* calculating coeffs */
  get_Ylm_coeffs(realClm,imagClm,Rnew_NS,Ntheta,Nphi,lmax);
  GridParams->NS_R_Ylm->realClm = realClm;
  GridParams->NS_R_Ylm->realClm = imagClm;
  GridParams->NS_R_Ylm->Lmax    = lmax;
  
  free(Rnew_NS);
  
  printf("Finding the surface of NS, Ylm method ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}
#ifdef ij
#undef ij
#endif

/* given theta, phi and the patch, it finds the corresponding
// patch and X,Y,Z coordinate. */
static void find_XYZ_and_patch_of_theta_phi_NS_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid)
{
  const double tgphi    = tan(phi);
  const double costheta = cos(theta);
  Flag_T found_flg = NO;
  unsigned p;
  
  /* check all of NS patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistence */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSPatch(patch))
      continue;
    if (strstr(patch->name,"left_centeral_box"))
      continue;

    Flag_T side = patch->CoordSysInfo->CubedSphericalCoord->side;
    const double *c = patch->c;
    double a,b;
    double a_sign,b_sign,c_sign;
    double x[3],phi2,theta2,r;
    
    X[2] = 1;/* since we are on NS surface */
    
    switch (side)
    {
      case UP:
        a = Sqrt(1 - Power(costheta,2))/
            Sqrt(Power(costheta,2) + Power(costheta,2)*Power(tgphi,2));
        b = (Sqrt(1 - Power(costheta,2))*tgphi)/
            Sqrt(Power(costheta,2)*(1 + Power(tgphi,2)));   
      break;
      case DOWN:
        b = Sqrt(1 - Power(costheta,2))/
            Sqrt(Power(costheta,2) + Power(costheta,2)*Power(tgphi,2));
        a = (Sqrt(1 - Power(costheta,2))*tgphi)/
            Sqrt(Power(costheta,2)*(1 + Power(tgphi,2)));
      break;
      case LEFT:
        a = 1/tgphi;
        b = Sqrt((-Power(costheta,2) - Power(costheta,2)*Power(tgphi,2))/
            ((-1 + Power(costheta,2))*Power(tgphi,2)));
      break;
      case RIGHT:
        b = 1/tgphi;
        a = Sqrt((-Power(costheta,2) - Power(costheta,2)*Power(tgphi,2))/
            ((-1 + Power(costheta,2))*Power(tgphi,2)));
      break;
      case BACK:
        a = Sqrt(-Power(costheta,2) - Power(costheta,2)*Power(tgphi,2))/
            Sqrt(-1 + Power(costheta,2));
        b = tgphi;
      break;
      case FRONT:
        b = Sqrt(-Power(costheta,2) - Power(costheta,2)*Power(tgphi,2))/
            Sqrt(-1 + Power(costheta,2));
        a = tgphi;
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* having found the magnitude of a and b, we need to find out the sign of them.
    // this is done by pay attention to side, signum(costheta) and range of tanphi */
    switch (side)
    {
      case UP:
        arctan_argument_signum(&b_sign,&a_sign,phi);
      break;
      case DOWN:
        arctan_argument_signum(&a_sign,&b_sign,phi);
      break;
      case LEFT:
        arctan_argument_signum(&c_sign,&a_sign,phi);
        assert(c_sign < 0);
        if (costheta > 0) b_sign = 1;
        else		b_sign = -1;
      break;
      case RIGHT:
        arctan_argument_signum(&c_sign,&b_sign,phi);
        assert(c_sign > 0);
        if (costheta > 0) a_sign = 1;
        else		a_sign = -1;
      break;
      case BACK:
        arctan_argument_signum(&b_sign,&c_sign,phi);
        assert(c_sign < 0);
        if (costheta > 0) a_sign = 1;
        else		a_sign = -1;
      break;
      case FRONT:
        arctan_argument_signum(&a_sign,&c_sign,phi);
        assert(c_sign > 0);
        if (costheta > 0) b_sign = 1;
        else		b_sign = -1;
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    X[0] = fabs(a)*a_sign;
    X[1] = fabs(b)*b_sign;
    
    /* check if x of X really gives you the correct angles */
    x_of_X(x,X,patch);
    x[0] -= c[0];
    x[1] -= c[1];
    x[2] -= c[2];
    r = rms(3,x,0);
    theta2 = acos(x[2]/r);
    phi2   = arctan(x[1],x[0]);
    if (EQL(theta2,theta) && EQL(phi2,phi))
    {
      found_flg = YES;
      *ppatch = patch;
      break;
    }
  }
  if (found_flg == NO)
    abortEr("(X,Y,Z) or patch could not be found.\n");
}

/* extrapolating phi, dphi and W in NS surrounding coords 
// in case they are needed for interpolation to the next grid
// or in calculation of enthalpy at NS surrounding patches.
// for extrapolation we demand:
// f(r) = f(r1)/(e^{-1}-e^{-r2/r1})*(e^{-r/r1}-e^{-r2/r1}), 
// where r1 is radius of NS surface and r2 is twice of r1
// note: f(r) = 0, if r >= r2 */
static void extrapolate_fluid_fields_outsideNS_CS(Grid_T *const grid)
{
  const double BN_NS_d = GetParameterD_E("BH_NS_separation");
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if (!IsItNSSurroundingPatch(patch))
      continue;
     
    /* irrotational part of fluid */
    ADD_FIELD(phi)
    GET_FIELD(phi)
    ADD_FIELD_NoMem(dphi_D2)
    ADD_FIELD_NoMem(dphi_D1)
    ADD_FIELD_NoMem(dphi_D0)
    DECLARE_AND_EMPTY_FIELD(dphi_D2)
    DECLARE_AND_EMPTY_FIELD(dphi_D1)
    DECLARE_AND_EMPTY_FIELD(dphi_D0)
    
    /* spin part of fluid W^i */
    ADD_FIELD(W_U0)
    ADD_FIELD(W_U1)
    ADD_FIELD(W_U2)
    GET_FIELD(W_U0)
    GET_FIELD(W_U1)
    GET_FIELD(W_U2)
    
    /* find the corresponding NS patch to be used for extrapolation */
    Patch_T *NS_patch;/* corresponding NS patch, to extrapolate out */
    char stem[1000];
    char *affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
    
    assert(affix);
    sprintf(stem,"left_NS%s",affix);
    free(affix);
    NS_patch = GetPatch(stem,grid);
    
    /* populating phi and W fields */
    Interpolation_T *interp_phi = init_interpolation();
    interp_phi->field = NS_patch->pool[LookUpField_E("phi",NS_patch)];
    interp_phi->XY_dir_flag = 1;
    interp_phi->K     = NS_patch->n[2]-1;
    plan_interpolation(interp_phi);
    
    Interpolation_T *interp_W_U0 = init_interpolation();
    interp_W_U0->field = NS_patch->pool[LookUpField_E("W_U0",NS_patch)];
    interp_W_U0->XY_dir_flag = 1;
    interp_W_U0->K     = NS_patch->n[2]-1;
    plan_interpolation(interp_W_U0);
    
    Interpolation_T *interp_W_U1 = init_interpolation();
    interp_W_U1->field = NS_patch->pool[LookUpField_E("W_U1",NS_patch)];
    interp_W_U1->XY_dir_flag = 1;
    interp_W_U1->K     = NS_patch->n[2]-1;
    plan_interpolation(interp_W_U1);
    
    Interpolation_T *interp_W_U2 = init_interpolation();
    interp_W_U2->field = NS_patch->pool[LookUpField_E("W_U2",NS_patch)];
    interp_W_U2->XY_dir_flag = 1;
    interp_W_U2->K     = NS_patch->n[2]-1;
    plan_interpolation(interp_W_U2);
    
    const unsigned *n = patch->n;
    unsigned ijk,i,j,k;
    double x[3],X[3];
    k = 0;/* at the surface of NS surface */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        ijk = L(n,i,j,k);
        /* find the value of phi at NS surface */
        X_of_x(X,patch->node[ijk]->x,NS_patch);
        interp_phi->X = X[0];
        interp_phi->Y = X[1];
        phi[ijk]  = execute_interpolation(interp_phi);
        W_U0[ijk] = execute_interpolation(interp_W_U0);
        W_U1[ijk] = execute_interpolation(interp_W_U1);
        W_U2[ijk] = execute_interpolation(interp_W_U2);
      }
    }
    
    /* since cubed spherical coords are radial in Z direction
    // and angular in X and Y direction we have: */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        ijk = L(n,i,j,0);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2];
        
        double phi0  = phi[ijk];
        double W0_U0 = W_U0[ijk];
        double W0_U1 = W_U1[ijk];
        double W0_U2 = W_U2[ijk];
        double r1 = rms(3,x,0);
        double r2 = 2*r1;
        assert(LSS(r2,BN_NS_d));
        
        for (k = 1; k < n[2]; ++k)
        {
         double r,fac;
         ijk = L(n,i,j,k);
         x[0] = patch->node[ijk]->x[0]-patch->c[0];
         x[1] = patch->node[ijk]->x[1]-patch->c[1];
         x[2] = patch->node[ijk]->x[2]-patch->c[2];
         r    = rms(3,x,0);
         
         if (GRTEQL(r,r2))
           fac     = 0;
         else
           fac     = (exp(-r/r1)-exp(-r2/r1))/(1./M_E-exp(-r2/r1));
           
         phi[ijk]  = phi0*fac;
         W_U0[ijk] = W0_U0*fac;
         W_U1[ijk] = W0_U1*fac;
         W_U2[ijk] = W0_U2*fac;
        }
      }
    }
    Field_T *phi_field = patch->pool[Ind("phi")];
    dphi_D2->v = Partial_Derivative(phi_field,"z");
    dphi_D1->v = Partial_Derivative(phi_field,"y");
    dphi_D0->v = Partial_Derivative(phi_field,"x");
  }
}
/* use TOV and Kerr-Schil black hole approximation.
// ->return value: resultant grid from this approximation */
static Grid_T *TOV_KerrShild_approximation(void)
{
  Grid_T *grid = 0;
  struct Grid_Params_S GridParams[1] = {0};/* adjust some pars for construction of grid */
  
  /* solve fields for a TOV star located at left side of y axis */
  TOV_T *tov = TOV_init();
  tov->bar_m = GetParameterD_E("NS_baryonic_mass");
  tov->description = "Estimating NS";
  tov = TOV_solution(tov);
  const double ns_R = tov->rbar[tov->N-1];
  
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
  GridParams->Max_R_NS_l = ns_R;
  GridParams->R_BH_r     = bh_R;
  GridParams->a_BH       = bh_chi*bh_mass;
  GridParams->NS_R_type  = "PefectSphere";
  grid = creat_grid_CS(GridParams);
  
  /* creating all of the fields needed for construction of Initial Data */
  bbn_allocate_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  bbn_populate_free_data(grid);
  
  /* initialize the fields using TOV and Kerr-Shild solution */
  init_field_TOV_plus_KerrSchild(grid,tov,bh_chi*bh_mass,bh_mass);
  
  /* taking partial derivatives of the fields needed for equations */
  bbn_partial_derivatives_fields(grid);
  
  /* update u0, _J^i, _E and _S */
  Tij_IF_CTS_psi6Sources(grid);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  bbn_update_Aij(grid);
  
  /* make normal vectorn on BH horizon */
  make_normal_vector_on_BH_horizon(grid);
  
  /* find Euler equation const using enthalpy of TOV star and other fields */
  find_Euler_eq_const_TOV_KerrSchild(grid);
  
  TOV_free(tov);
  
  return grid;
}

/* find Euler equation const using enthalpy of TOV star and 
// other fields that are found in TOV Kerr-Schild approximation.
// this is important, since otherwise when h gets updated, it may
// become 'nan' due to incorrect value of Euler equation constant.
//
// FORMULA I USED: 
// C+h/u0+D_{i}phi*(D^{i}phi+W^{i})/(h*u0)-Beta^{i}*D_{i}phi = 0.
// 
// CPI SCRIPT that I used:
 
Declare = 
{

 # enthalpy
 (obj = Field,name = enthalpy, rank = 0, C_macro);

 # conformal metric inverse
 (obj = Field,name = _gammaI, rank = UU, C_macro);

 # spin part of fluid
 (obj = Field,name = W, rank = U, C_macro);

 # d(phi)/d? for irrotional part of fluid
 (obj = Field,name = dphi, rank = D, C_macro);

 # Beta
 (obj = Field,name = Beta, rank = U, C_macro);

 # conformal factor
 (obj = Field,name = psi, rank = 0, C_macro);

 # u0
 (obj = Field,name = u0, rank = 0, C_macro);
}
# symmetries:
Symm[_gammaI(i,j)  = _gammaI(j,i)];

psim4 = psi**(-4);
dphiP = psim4*_gammaI(i,j)*dphi(-i)*dphi(-j)+dphi(-i)*W(i);
Euler_C = -enthalpy/u0 - dphiP/(enthalpy*u0) + Beta(i)*dphi(-i);

*/
static void find_Euler_eq_const_TOV_KerrSchild(Grid_T *const grid)
{
  unsigned p,ijk;
  double Euler_C;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if(!IsItNSPatch(patch))
      continue;
    
    ijk = 0;/* for an arbitrary point */
    GET_FIELD(enthalpy)
    GET_FIELD(_gammaI_U0U2)
    GET_FIELD(_gammaI_U0U0)
    GET_FIELD(_gammaI_U0U1)
    GET_FIELD(_gammaI_U1U2)
    GET_FIELD(_gammaI_U1U1)
    GET_FIELD(_gammaI_U2U2)
    GET_FIELD(W_U1)
    GET_FIELD(W_U0)
    GET_FIELD(W_U2)
    GET_FIELD(dphi_D2)
    GET_FIELD(dphi_D1)
    GET_FIELD(dphi_D0)
    GET_FIELD(Beta_U1)
    GET_FIELD(Beta_U0)
    GET_FIELD(Beta_U2)
    GET_FIELD(psi)
    GET_FIELD(u0)

    double psim4 = 
pow(psi[ijk], -4);

    double dphiP = 
W_U0[ijk]*dphi_D0[ijk] + W_U1[ijk]*dphi_D1[ijk] + W_U2[ijk]*
dphi_D2[ijk] + psim4*(_gammaI_U0U0[ijk]*pow(dphi_D0[ijk], 2) + 2.0*
_gammaI_U0U1[ijk]*dphi_D0[ijk]*dphi_D1[ijk] + 2.0*_gammaI_U0U2[ijk]*
dphi_D0[ijk]*dphi_D2[ijk] + _gammaI_U1U1[ijk]*pow(dphi_D1[ijk], 2) +
2.0*_gammaI_U1U2[ijk]*dphi_D1[ijk]*dphi_D2[ijk] + _gammaI_U2U2[ijk]*
pow(dphi_D2[ijk], 2));

    Euler_C = 
(-dphiP - pow(enthalpy[ijk], 2) + enthalpy[ijk]*u0[ijk]*(Beta_U0[ijk]*
dphi_D0[ijk] + Beta_U1[ijk]*dphi_D1[ijk] + Beta_U2[ijk]*dphi_D2[ijk]))/
(enthalpy[ijk]*u0[ijk]);
    
    break;/* only for 1 patch we find Euler constant */
  }
  update_parameter_double_format("Euler_equation_constant",Euler_C);
}

/* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
// _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
static void bbn_update_Aij(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ...\n");
  unsigned p;

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    bbn_update_psi10A_UiUj(patch);
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
  const double Omega_NS_x = GetParameterD_E("NS_Omega_U0");
  const double Omega_NS_y = GetParameterD_E("NS_Omega_U1");
  const double Omega_NS_z = GetParameterD_E("NS_Omega_U2");
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
  
  /* initialization psi, eta and matter fields: */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    PREP_FIELD(psi)
    PREP_FIELD(eta)
    PREP_FIELD(KSalpha)
    
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
      
      PREP_FIELD(enthalpy)
      PREP_FIELD(rho0)
      PREP_FIELD(phi)
      PREP_FIELD(W_U0)
      PREP_FIELD(W_U1)
      PREP_FIELD(W_U2)
      
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
        phi[ijk] = 0;//-Omega_BHNS*(y-y_CM)*x;
        
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

/* given the max of radius of NS, radius of BH and their separation,
// and BH spin as the parameters, create a grid with these properties.
// NOTE: WE assume we are using cubed spherical grid.
// ->return value: grid of NS and BH in which inside of the BH excised. */
static Grid_T *creat_grid_CS(struct Grid_Params_S *const GridParams)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  /* calculate the characteristics of this grid */
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const double R_BH_r     = GridParams->R_BH_r;
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
  if (!strcmp_i(kind,"BBN_CubedSpherical_grid"))
    abortEr("This function only works with cubed spherical grid.\n");
    
  grid->kind = dup_s(kind);
  
  assert(GRT(C,0));
  assert(GRT(Max_R_NS_l,0));
  assert(GRT(R_BH_r,0));
  assert(LSS(2*Max_R_NS_l,C));
  assert(LSS(2*R_BH_r,C));
  
  /* making NS and BH surfaces function */
  NS_BH_surface_CubedSpherical_grid(grid,GridParams);
  
  box_size_l = GetParameterD_E("left_central_box_length_ratio")*Max_R_NS_l;
  
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
static void NS_BH_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const double R_BH_r     = GridParams->R_BH_r;
  const double a_BH       = GridParams->a_BH;
  double *R;
  char par[1000] = {'\0'};
  unsigned N[3],n,i,j,k,N_total;
  Patch_T patch[1] = {0};
  struct Collocation_s coll_s[2] = {0};
  double X[3],r;
  
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
  /* if NS is perfect sphere like TOV*/
  if (strcmp_i(GridParams->NS_R_type,"PrefectSphere"))
  {
    for (i = 0; i < N[0]; ++i)
      for (j = 0; j < N[1]; ++j)
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = Max_R_NS_l;
  }
  /* if NS radius is varied and we know its expansion in Ylm bases */
  else if (strcmp_i(GridParams->NS_R_type,"SphericalHarmonic"))
  {
    /* we need interpolation */
    double *const realClm = GridParams->NS_R_Ylm->realClm;
    double *const imagClm = GridParams->NS_R_Ylm->imagClm;
    const unsigned Lmax   = GridParams->NS_R_Ylm->Lmax;
    double theta,phi;
    
    X[3] = 1;/* since we are on the NS surface */
    
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
      
    patch->n[0] = N[0];
    patch->n[1] = N[1];
    patch->n[2] = N[2];
  
    initialize_collocation_struct(patch,&coll_s[0],0);
    initialize_collocation_struct(patch,&coll_s[1],1);
    
    /* surface up */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,UP);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_up",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface down */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,DOWN);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_down",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface back */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,BACK);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_back",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface front */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,FRONT);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_front",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface left */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,LEFT);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_left",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* surface right */
    for (i = 0; i < N[0]; ++i)
    {
      X[0] = point_value(i,&coll_s[0]);
      for (j = 0; j < N[1]; ++j)
      {
        X[1] = point_value(j,&coll_s[1]);
        find_theta_phi_of_XYZ_NS_CS(&theta,&phi,X,RIGHT);
        r = interpolation_Ylm(realClm,imagClm,Lmax,theta,phi);
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = r;
      }
    }
    sprintf(par,"grid%u_left_NS_surface_function_right",grid->gn);
    add_parameter_array(par,R,N_total);
    
    /* free coeffs */
    free(realClm);
    free(imagClm);
  }
  else
    abortEr(NO_OPTION);
  
  free(R);
  
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

/* given (X,Y,Z) in the specified slice of NS in cubed spherical coords
// it finds the associated polar and azimuthal angels */
static void find_theta_phi_of_XYZ_NS_CS(double *const theta,double *const phi,const double *const X,const Flag_T side)
{
  const double a = X[0];
  const double b = X[1];
  const double d = sqrt(1+SQR(a)+SQR(b));
  
  switch (side)
  {
    case UP:
      *phi   = arctan(b,a);
      *theta = acos(1/d);
    break;
    case DOWN:
      *phi   = arctan(a,b);
      *theta = acos(-1/d);
    break;
    case LEFT:
      *phi   = arctan(-1,a);
      *theta = acos(b/d);
    break;
    case RIGHT:
      *phi   = arctan(1,b);
      *theta = acos(a/d);
    break;
    case BACK:
      *phi   = arctan(b,-1);
      *theta = acos(a/d);
    break;
    case FRONT:
      *phi   = arctan(a,1);
      *theta = acos(b/d);
    break;
    default:
      abortEr(NO_OPTION);
  }
  
}
