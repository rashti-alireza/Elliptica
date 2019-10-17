/*
// Alireza Rashti
// October 2019
*/

#include "sns_initialize.h"

/* initialize this system according to the previous grid. 
// ->return value: the next grid as a result of this initialization. */
Grid_T *sns_initialize_next_grid(Grid_T *const grid_prev)
{
  Grid_T *grid_next = 0;
  
  if (!grid_prev)/* if grid is empty come up with an approximation */
  {
    /* if we use TOV and Kerr-Schil black hole approximation */
    if (strcmp_i(GetParameterS_E("NS_initialization"),"TOV"))
      grid_next = TOV_approximation();
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
  
  /* updating enthalpy */
  sns_update_enthalpy_and_denthalpy(grid_prev);
  
  if (strcmp_i(grid_prev->kind,"SNS_CubedSpherical+Box_grid"))
  {
    /* extrapolate fluid fields outside of NS in case their value needed. */
    extrapolate_fluid_fields_outsideNS_CS(grid_prev);
  }
  else
    abortEr(NO_OPTION);
    
  /* find the NS center */
  find_NS_center(grid_prev);
  
  /* if needed, drag the NS to its designated point */
  adjust_NS_center(grid_prev);
  
  /* find Euler equation constant to meet NS baryonic mass */
  find_Euler_eq_const(grid_prev);
    
  /* find NS surface */
  if (strcmp_i(grid_prev->kind,"SNS_CubedSpherical+Box_grid"))
  {
    
    GridParams->NS_R_type = GetParameterS_E("NS_surface_finder_method");
    /* find NS surface using spherical harmonic points */
    if (strstr_i(GridParams->NS_R_type,"SphericalHarmonic"))
      find_NS_surface_Ylm_method_CS(grid_prev,GridParams);
    else
      abortEr(NO_OPTION);
  }
  else
    abortEr(NO_OPTION);
    
  /* make new grid */  
  grid_next = creat_sns_grid_CS(GridParams);
  
  /* fields: */
  /* creating all of the fields needed for construction of Initial Data */
  sns_allocate_fields(grid_next);
  
  /* populating the free data part of initial data that we chose ourself */
  sns_populate_free_data(grid_next);

  /* use previous grid to interpolate values of the fields for 
  // the next grid  and initialzing some other fields */
  interpolate_and_initialize_to_next_grid(grid_next,grid_prev);
  
  /* taking partial derivatives of the fields needed for equations */
  sns_partial_derivatives_fields(grid_next);
  //sns_study_initial_data(grid_next);
  
  /* update enthalpy,denthalpy,rho0, drho0, u0, _J^i, _E and _S */
  sns_update_matter_fields(grid_next);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  sns_update_Aij(grid_next);  
  
  /* freeing */
  free_Grid_Params_S(GridParams);
  
  return grid_next;
}

/* find Euler equation constant to meet NS baryonic mass */
static void find_Euler_eq_const(Grid_T *const grid)
{
  Root_Finder_T *root = init_root_finder(1);
  double *Euler_const = 0;
  double guess[1];/* initial guess for Euler const */
  struct Euler_eq_const_RootFinder_S params[1];
  
  params->grid = grid;
  params->NS_baryonic_mass = GetParameterD_E("NS_baryonic_mass");
  guess[0] = GetParameterD_E("Euler_equation_constant");
  
  root->description = "Finding Euler equation constant using NS baryonic mass";
  root->type        = GetParameterS_E("RootFinder_Method");
  root->tolerance   = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->f[0]        = Euler_eq_const_rootfinder_eq;
  root->verbose     = 1;
  plan_root_finder(root);
  Euler_const       = execute_root_finder(root);
  update_parameter_double_format("Euler_equation_constant",Euler_const[0]);
  printf("Euler Equation const. updated: %g -> %g\n",guess[0],Euler_const[0]);
  free(Euler_const);
  free_root_finder(root);
}

/* root finder eqution for Euler equation constant */
static double Euler_eq_const_rootfinder_eq(void *params,const double *const x)
{
  struct Euler_eq_const_RootFinder_S *const par = params;
  
  return sns_NS_baryonic_mass(par->grid,x[0]) - par->NS_baryonic_mass;
}

/* use previous grid to interpolate values of the fields that will be solved for the next grid */
static void interpolate_and_initialize_to_next_grid(Grid_T *const grid_next,Grid_T *const grid_prev)
{
  const unsigned np = grid_next->np;
  unsigned p;
 
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
    make_coeffs_3d(B0_U0);
    make_coeffs_3d(B0_U1);
    make_coeffs_3d(B0_U2);
    make_coeffs_3d(psi);
    make_coeffs_3d(eta);
    
    if (IsItNSPatch(patch))
    {
      DECLARE_FIELD(phi)
      DECLARE_FIELD(enthalpy)
      make_coeffs_3d(phi);
      make_coeffs_3d(enthalpy);
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
    root_name++;
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
        
        B0_U0[ijk]    = interpolate_from_patch_prim("B0_U0",Xp,patchp);
        B0_U1[ijk]    = interpolate_from_patch_prim("B0_U1",Xp,patchp);
        B0_U2[ijk]    = interpolate_from_patch_prim("B0_U2",Xp,patchp);
        psi[ijk]      = interpolate_from_patch_prim("psi",Xp,patchp);
        eta[ijk]      = interpolate_from_patch_prim("eta",Xp,patchp);
        phi[ijk]      = interpolate_from_patch_prim("phi",Xp,patchp);
        enthalpy[ijk] = interpolate_from_patch_prim("enthalpy",Xp,patchp);
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
        
        B0_U0[ijk] = interpolate_from_patch_prim("B0_U0",Xp,patchp);
        B0_U1[ijk] = interpolate_from_patch_prim("B0_U1",Xp,patchp);
        B0_U2[ijk] = interpolate_from_patch_prim("B0_U2",Xp,patchp);
        psi[ijk]   = interpolate_from_patch_prim("psi",Xp,patchp);
        eta[ijk]   = interpolate_from_patch_prim("eta",Xp,patchp);
      }
    }
  }/* end of for (p = 0; p < np; ++p) */
  
  /* initializing some other fields: */
  /* W_U[0-2],Beta_U[0-2],B1_U[0-2] */
  const double Omega_NS_x = GetParameterD_E("NS_Omega_U0");
  const double Omega_NS_y = GetParameterD_E("NS_Omega_U1");
  const double Omega_NS_z = GetParameterD_E("NS_Omega_U2");
  const double C_NS = GetParameterD_E("NS_Center");

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
       
      /* B1 */
      B1_U0[ijk] = 0;
      B1_U1[ijk] = 0;
      B1_U2[ijk] = 0;
       
      /* Beta */
      Beta_U0[ijk] = B0_U0[ijk]+B1_U0[ijk];
      Beta_U1[ijk] = B0_U1[ijk]+B1_U1[ijk];
      Beta_U2[ijk] = B0_U2[ijk]+B1_U2[ijk];
    }
    
    if (IsItNSPatch(patch))
    {
      
      PREP_FIELD(W_U0)
      PREP_FIELD(W_U1)
      PREP_FIELD(W_U2)
      
      for (ijk = 0; ijk < nn; ++ijk)
      {
        double x = patch->node[ijk]->x[0];
        double y = patch->node[ijk]->x[1]-C_NS;
        double z = patch->node[ijk]->x[2];
        
        /* spin part */
        W_U0[ijk] = Omega_NS_y*z-Omega_NS_z*y;
        W_U1[ijk] = Omega_NS_z*x-Omega_NS_x*z;
        W_U2[ijk] = Omega_NS_x*y-Omega_NS_y*x;
      }
    }/* end of if (IsItNSPatch(patch)) */
    
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
  if (needle->Nans)
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

#define ij(i,j) ((j)+Nphi*(i))
/* given the grid, find the NS surface on Ylm points 
// i.e (theta,phi) collocations are = (Legendre,EquiSpaced),
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
  struct NS_surface_RootFinder_S par[1];
  unsigned Ntheta,Nphi;/* total number of theta and phi points */
  const unsigned lmax = (unsigned)GetParameterI_E("NS_surface_Ylm_expansion_max_l");
  double theta,phi;
  double *Rnew_NS = 0;/* new R for NS */
  double Max_R_NS = 0;/* maximum radius of NS */
  double guess    = 0;
  const double maxR_out     = (1./3.)*(-2*GetParameterD_E("NS_center_-y_axis"));
  const double SMALL_FACTOR = 0.8;/* initial r = SMALL_FACTOR*R0_NS */
  const double RESIDUAL     = 1E3*GetParameterD_E("RootFinder_Tolerance");
  const double HowBigFac    = 1.5;
  const double HowSmallFac  = 0.5;
  const double INCREMENT    = 1.2;/* %20 increase */
  const double DECREMENT    = 0.8;/* %20 decrease */
  double X[3],x[3],N[3];
  char stem[1000],*affix;
  //Flag_T NS_patch_flg = NONE;
  unsigned i,j;
  
  par->Euler_C = GetParameterD_E("Euler_equation_constant");
  par->scale   = 1E-2;
  
  /* populate root finder */
  Root_Finder_T *root = init_root_finder(1);
  root->type      = GetParameterS_E("RootFinder_Method");
  plan_root_finder(root);
  root->tolerance = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter   = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss     = &guess;
  root->params    = par;
  root->f[0]      = sns_NS_surface_enthalpy_eq;
  root->FD_Right  = 1;
  //root->verbose  = 1;
  par->root_finder = root;
  
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
      
      //NS_patch_flg = NONE;
      Patch_T *h_patch = 0,*patch = 0;
      double y1[3] = {0},y2[3] = {0};
      double h,R0_NS,*dr;
      
      /* find patch and X,Y,Z at NS surface in which theta and phi take place */
      find_XYZ_and_patch_of_theta_phi_NS_CS(X,&patch,theta,phi,grid);
      
      /* find enthalpy at the (X,Y,Z) */
      Interpolation_T *interp_h = init_interpolation();
      interp_h->field = patch->pool[Ind("enthalpy")];
      interp_h->XY_dir_flag  = 1;
      interp_h->X            = X[0];
      interp_h->Y            = X[1];
      interp_h->K            = patch->n[2]-1;
      plan_interpolation(interp_h);
      h = execute_interpolation(interp_h);/* enthalpy */
      free_interpolation(interp_h);
      
      /* finding x */
      x_of_X(x,X,patch);
      
      /* r^ = sin(theta)cos(phi)x^+sin(theta)sin(phi)y^+cos(theta)z^ */
      N[0] = sin(theta)*cos(phi);
      N[1] = sin(theta)*sin(phi);
      N[2] = cos(theta);
        
      /* y is where the previous radius was located */
      y1[0] = x[0]-patch->c[0];
      y1[1] = x[1]-patch->c[1];
      y1[2] = x[2]-patch->c[2];
      R0_NS = rms(3,y1,0);
      
      assert(h > 0);
      root->interrupt = 0;
      if(LSSEQL(h,1))/* if it takes place at NS patch */
      {
        //NS_patch_flg = YES;
        h_patch      = patch;
        
        y2[0] = SMALL_FACTOR*R0_NS*N[0];
        y2[1] = SMALL_FACTOR*R0_NS*N[1];
        y2[2] = SMALL_FACTOR*R0_NS*N[2];
        
        par->x0[0] = y2[0]+h_patch->c[0];
        par->x0[1] = y2[1]+h_patch->c[1];
        par->x0[2] = y2[2]+h_patch->c[2];
        par->maxR = R0_NS;
        
      }
      else/* which means h = 1 occures in neighboring patch */
      {
        //NS_patch_flg = NO;
        /* finding the juxtapose patch of this NS patch */
        affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);
        assert(affix);
        sprintf(stem,"left_NS_surrounding%s",affix);
        free(affix);
        
        h_patch = GetPatch(stem,grid);
        
        y2[0] = x[0]-patch->c[0];
        y2[1] = x[1]-patch->c[1];
        y2[2] = x[2]-patch->c[2];
        
        par->x0[0] = x[0];
        par->x0[1] = x[1];
        par->x0[2] = x[2];
        par->maxR = maxR_out;
        guess = 1;
      }
      /* having found h_patch, now find the the position of h = 1 */
      par->patch = h_patch;
      par->N     = N;
      
      dr = execute_root_finder(root);
      //printf("---- %s\n",h_patch->name);
      //print_root_finder_exit_status(root);
      /* if root finder went wrong for some reason */
      if (root->interrupt || GRT(root->residual,RESIDUAL))
      {
        //print_root_finder_exit_status(root);
        Rnew_NS[ij(i,j)] = R0_NS;
        printf("---- %s\n",h_patch->name);
        print_root_finder_exit_status(root);
      }
      else/* if root finder was successful */
      {
        /*  new coords of R respect to the center of NS */
        y2[0] += N[0]*dr[0]*par->scale;
        y2[1] += N[1]*dr[0]*par->scale;
        y2[2] += N[2]*dr[0]*par->scale;
        Rnew_NS[ij(i,j)] = rms(3,y2,0);
      }
      
      /* if we have big difference in NS radius, ameliorate it */
      if (GRTEQL(Rnew_NS[ij(i,j)],HowBigFac*R0_NS))
        Rnew_NS[ij(i,j)] = R0_NS*INCREMENT;
      else if (LSSEQL(Rnew_NS[ij(i,j)],HowSmallFac*R0_NS))
        Rnew_NS[ij(i,j)] = R0_NS*DECREMENT;
      
      if (Rnew_NS[ij(i,j)] > Max_R_NS)
        Max_R_NS = Rnew_NS[ij(i,j)];
        
      free(dr);
        
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  /* adding maximum radius of NS to grid parameters */
  GridParams->Max_R_NS_l = Max_R_NS;
  
  /* making radius of NS parameter at each patch using Ylm interpolation */
  double *realClm = alloc_ClmYlm(lmax);
  double *imagClm = alloc_ClmYlm(lmax);
  
  /* calculating coeffs */
  get_Ylm_coeffs(realClm,imagClm,Rnew_NS,Ntheta,Nphi,lmax);
  GridParams->NS_R_Ylm->realClm = realClm;
  GridParams->NS_R_Ylm->imagClm = imagClm;
  GridParams->NS_R_Ylm->Lmax    = lmax;
  
  free(Rnew_NS);
  free_root_finder(root);
  
  printf("Finding the surface of NS, Ylm method ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}
#ifdef ij
#undef ij
#endif

/* given theta, phi and knowing the fact that they are on NS surface, 
// it finds the corresponding patch and X,Y,Z coordinate. */
static void find_XYZ_and_patch_of_theta_phi_NS_CS(double *const X,Patch_T **const ppatch,const double theta,const double phi,Grid_T *const grid)
{
  const double tan_phi    = tan(phi);
  const double cos_theta  = cos(theta);
  const double tan_phi2   = SQR(tan_phi);
  const double cos_theta2 = SQR(cos_theta);
  Flag_T found_flg = NO;
  unsigned p;
  
  X[2] = 1;/* since we are on NS surface */
  
  /* check all of NS patches in which (x,y,z) and 
  // (X,Y,Z) and (theta,phi) are consistent */
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
    
    /* we know that theta = 0 or Pi occures only at UP or DOWN patches
    // so don't bother to follow algorithm all the way down.
    // furthermore, this prevent 0 in denominator of unrelated patches. */
    if (EQL(theta,0) || EQL(theta,M_PI))
    {
      if (side == LEFT || side == RIGHT || 
          side == BACK || side == FRONT   )
        continue;
    }
    
    /* first calculate the magnetitude of a and b 
    // which are related to X[0] and X[1] with a sign */
    switch (side)
    {
      case UP:
        a = Sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        b = tan_phi*Sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case DOWN:
        b = Sqrt((1 - cos_theta2)/(cos_theta2 + cos_theta2*tan_phi2));
        a = tan_phi*Sqrt((1 - cos_theta2)/(cos_theta2*(1 + tan_phi2)));
      break;
      case LEFT:
        a = 1/tan_phi;
        b = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case RIGHT:
        b = 1/tan_phi;
        a = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/((1. - cos_theta2)*tan_phi2));
      break;
      case BACK:
        a = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        b = tan_phi;
      break;
      case FRONT:
        b = Sqrt((cos_theta2 + cos_theta2*tan_phi2)/(1 - cos_theta2));
        a = tan_phi;
      break;
      default:
        abortEr(NO_OPTION);
    }
    
    /* having found the magnitude of a and b, we need to find out the sign of them.
    // this is done by paying attention to side, signum(cos_theta) and range of tanphi */
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
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
      break;
      case RIGHT:
        arctan_argument_signum(&c_sign,&b_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case BACK:
        arctan_argument_signum(&b_sign,&c_sign,phi);
        if (cos_theta > 0) a_sign = 1;
        else		   a_sign = -1;
      break;
      case FRONT:
        arctan_argument_signum(&a_sign,&c_sign,phi);
        if (cos_theta > 0) b_sign = 1;
        else		   b_sign = -1;
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
// the fields spread out in the same fashion as it is changing toward 
// the NS sarface. what we have :
// f(r_out) = (f(r_in) + df)*exp(g(r)), in which df is f(r2)-f(r1), 
// r1 = FACTOR*r2 and r2 is the radius of NS surface, and g(r) is 
// a function of r_out/r2 to control the radial trend of fieldd.
// a = (r_max-r2)/(r2-r1)
// b = (r1*r_max-SQR(r2))/(r1-r2)
// r_out = a*r_in + b, where (r_in) r_out is r (inside)outside NS and
// r_max is the max radius in NS surrounding patch.
// in effect, it means we emulate the trend of the fields inside of NS
// from r1 to r2 in NS surrounding, from r2 to r_max. */
static void extrapolate_fluid_fields_outsideNS_CS(Grid_T *const grid)
{
  const double FACTOR = 0.8;/* r1 = FACTOR*r2 */
  //const double EXP    = 1;/* (r_out/r2)^EXP */
  unsigned p;
  
  FOR_ALL_PATCHES(p,grid)
  {
    /* surrounding patch */
    Patch_T *patch = grid->patch[p];
    if (!IsItNSSurroundingPatch(patch))
      continue;
     
    /* add fields: */
    /* enthalpy */
    ADD_FIELD(enthalpy)
    
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
    
    Patch_T *NS_patch;/* corresponding NS patch, to extrapolate out */
    const unsigned *n = patch->n;
    double r1,r2,r_max,r_in,r_out,a,b;
    double phii,W_U0i,W_U1i,W_U2i;/* fields at r1 */
    double phif,W_U0f,W_U1f,W_U2f;/* fields at r2 */
    double dphi,dW_U0,dW_U1,dW_U2;
    //double sphi,sW_U0,sW_U1,sW_U2;/* signs */
    double x[3],X[3];
    double THETA,PHI;
    unsigned ijk,i,j,k;

    /* find the corresponding NS patch to be used for extrapolation */
    char stem[1000];
    char *affix = regex_find("_[[:alpha:]]{2,5}$",patch->name);/* finding the side of the patch */
    assert(affix);
    sprintf(stem,"left_NS%s",affix);
    free(affix);
    NS_patch = GetPatch(stem,grid);
    
    /* prepare interpolation arguments */
    Interpolation_T *interp_phi = init_interpolation();
    interp_phi->field = NS_patch->pool[LookUpField_E("phi",NS_patch)];
    interp_phi->XYZ_dir_flag = 1;
    plan_interpolation(interp_phi);

    Interpolation_T *interp_W_U0 = init_interpolation();
    interp_W_U0->field = NS_patch->pool[LookUpField_E("W_U0",NS_patch)];
    interp_W_U0->XYZ_dir_flag = 1;
    plan_interpolation(interp_W_U0);

    Interpolation_T *interp_W_U1 = init_interpolation();
    interp_W_U1->field = NS_patch->pool[LookUpField_E("W_U1",NS_patch)];
    interp_W_U1->XYZ_dir_flag = 1;
    plan_interpolation(interp_W_U1);

    Interpolation_T *interp_W_U2 = init_interpolation();
    interp_W_U2->field = NS_patch->pool[LookUpField_E("W_U2",NS_patch)];
    interp_W_U2->XYZ_dir_flag = 1;
    plan_interpolation(interp_W_U2);

    /* spreading out the fields value using interpolation */
    for (i = 0; i < n[0]; ++i)
    {
      for (j = 0; j < n[1]; ++j)
      {
        /* calculate r_max at NS surrounding patch */
        ijk = L(n,i,j,n[2]-1);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2]; 
        r_max = rms(3,x,0);

        /* calculate r2 at NS surface */
        ijk = L(n,i,j,0);
        x[0] = patch->node[ijk]->x[0]-patch->c[0];
        x[1] = patch->node[ijk]->x[1]-patch->c[1];
        x[2] = patch->node[ijk]->x[2]-patch->c[2];
        r2   = rms(3,x,0);

        /* calulate r1,a,b */
        r1 = FACTOR*r2;/* small scale to define r1 */
        a  = (r_max-r2)/(r2-r1);
        b  = (r1*r_max-SQR(r2))/(r1-r2);

        /* find the value of phi and W^{i} at r2 */
        THETA = acos(x[2]/r2);
        PHI   = arctan(x[1],x[0]);

        x[0] = r2*sin(THETA)*cos(PHI)+NS_patch->c[0];
        x[1] = r2*sin(THETA)*sin(PHI)+NS_patch->c[1];
        x[2] = r2*cos(THETA)+NS_patch->c[2];

        X_of_x(X,x,NS_patch);

        interp_phi->X  = X[0];
        interp_phi->Y  = X[1];
        interp_phi->Z  = X[2];

        interp_W_U0->X = X[0];
        interp_W_U0->Y = X[1];
        interp_W_U0->Z = X[2];

        interp_W_U1->X = X[0];
        interp_W_U1->Y = X[1];
        interp_W_U1->Z = X[2];

        interp_W_U2->X = X[0];
        interp_W_U2->Y = X[1];
        interp_W_U2->Z = X[2];

        phif  = execute_interpolation(interp_phi);
        W_U0f = execute_interpolation(interp_W_U0);
        W_U1f = execute_interpolation(interp_W_U1);
        W_U2f = execute_interpolation(interp_W_U2);

        /* find the value of phi and W^{i} at r1 */
        x[0] = r1*sin(THETA)*cos(PHI)+NS_patch->c[0];
        x[1] = r1*sin(THETA)*sin(PHI)+NS_patch->c[1];
        x[2] = r1*cos(THETA)+NS_patch->c[2];

        X_of_x(X,x,NS_patch);
        
        interp_phi->X  = X[0];
        interp_phi->Y  = X[1];
        interp_phi->Z  = X[2];

        interp_W_U0->X = X[0];
        interp_W_U0->Y = X[1];
        interp_W_U0->Z = X[2];

        interp_W_U1->X = X[0];
        interp_W_U1->Y = X[1];
        interp_W_U1->Z = X[2];

        interp_W_U2->X = X[0];
        interp_W_U2->Y = X[1];
        interp_W_U2->Z = X[2];

        phii  = execute_interpolation(interp_phi);
        W_U0i = execute_interpolation(interp_W_U0);
        W_U1i = execute_interpolation(interp_W_U1);
        W_U2i = execute_interpolation(interp_W_U2);

        dphi  = phif  - phii;
        dW_U0 = W_U0f - W_U0i;
        dW_U1 = W_U1f - W_U1i;
        dW_U2 = W_U2f - W_U2i;
        //sphi  = dphi  > 0 ?  -1 : 1;
        //sW_U0 = dW_U0 > 0 ?  -1 : 1;
        //sW_U1 = dW_U1 > 0 ?  -1 : 1;
        //sW_U2 = dW_U2 > 0 ?  -1 : 1;
        
        for (k = 0; k < n[2]; ++k)
        {
          ijk = L(n,i,j,k);
          x[0] = patch->node[ijk]->x[0]-patch->c[0];
          x[1] = patch->node[ijk]->x[1]-patch->c[1];
          x[2] = patch->node[ijk]->x[2]-patch->c[2];
          
          r_out = rms(3,x,0);
          r_in  = (r_out-b)/a;

          x[0] = r_in*sin(THETA)*cos(PHI)+NS_patch->c[0];
          x[1] = r_in*sin(THETA)*sin(PHI)+NS_patch->c[1];
          x[2] = r_in*cos(THETA)+NS_patch->c[2];

          /* find the value of phi and W^{i} at NS surface */
          X_of_x(X,x,NS_patch);

          interp_phi->X  = X[0];
          interp_phi->Y  = X[1];
          interp_phi->Z  = X[2];

          interp_W_U0->X = X[0];
          interp_W_U0->Y = X[1];
          interp_W_U0->Z = X[2];

          interp_W_U1->X = X[0];
          interp_W_U1->Y = X[1];
          interp_W_U1->Z = X[2];

          interp_W_U2->X = X[0];
          interp_W_U2->Y = X[1];
          interp_W_U2->Z = X[2];
          
          //double e0 = pow(r_out/r2,EXP);
          
          phi[ijk]  = execute_interpolation(interp_phi)+dphi;
          //phi[ijk] *= exp(sphi*e0)*pow(r_out/r2,sphi);
          W_U0[ijk] = execute_interpolation(interp_W_U0)+dW_U0;
          //W_U0[ijk] *= exp(sW_U0*e0)*pow(r_out/r2,sW_U0);
          W_U1[ijk] = execute_interpolation(interp_W_U1)+dW_U1;
          //W_U1[ijk] *= exp(sW_U1*e0)*pow(r_out/r2,sW_U1);
          W_U2[ijk] = execute_interpolation(interp_W_U2)+dW_U2;
          //W_U2[ijk] *= exp(sW_U2*e0)*pow(r_out/r2,sW_U2);
          
        }/* end of for (k = 0 ; k < n[2]; ++k) */
      }/* end of for (j = 0; j < n[1]; ++j) */
    }/* end of for (i = 0; i < n[0]; ++i) */
    free_interpolation(interp_phi);
    free_interpolation(interp_W_U0);
    free_interpolation(interp_W_U1);
    free_interpolation(interp_W_U2);

    Field_T *phi_field = patch->pool[Ind("phi")];
    dphi_D2->v = Partial_Derivative(phi_field,"z");
    dphi_D1->v = Partial_Derivative(phi_field,"y");
    dphi_D0->v = Partial_Derivative(phi_field,"x");
    
    /* populating enthalpy in NS surroundings */
    Tij_IF_CTS_enthalpy(patch);
    
  }/* end of FOR_ALL_PATCHES(p,grid) */
}

/* use TOV and Kerr-Schil black hole approximation.
// ->return value: resultant grid from this approximation */
static Grid_T *TOV_approximation(void)
{
  Grid_T *grid = 0;
  struct Grid_Params_S *GridParams = init_GridParams();/* adjust some pars for construction of grid */
  
  /* solve fields for a TOV star located at left side of y axis */
  TOV_T *tov = TOV_init();
  tov->bar_m = GetParameterD_E("NS_baryonic_mass");
  tov->description = "Make a TOV Neutron Star";
  tov = TOV_solution(tov);
  const double ns_R = tov->rbar[tov->N-1];
  
  /* combining these two geometry to create the grid */
  GridParams->Max_R_NS_l = ns_R;
  GridParams->NS_R_type  = "PerfectSphere";
  grid = creat_sns_grid_CS(GridParams);
  
  /* creating all of the fields needed for construction of Initial Data */
  sns_allocate_fields(grid);
  
  /* populating the free data part of initial data that we chose ourself */
  sns_populate_free_data(grid);
  
  /* initialize the fields using TOV solution */
  init_field_TOV(grid,tov);
  
  /* taking partial derivatives of the fields needed for equations */
  sns_partial_derivatives_fields(grid);
  
  /* update u0, _J^i, _E and _S */
  Tij_IF_CTS_psi6Sources(grid);
  
  /* update _Aij in K^{ij} = A^{ij}+1/3*gamma^{ij}*K and 
  // _A^{ij} = gamma^10*A^{ij} and _dA^{ij} */
  sns_update_Aij(grid);
  
  /* find Euler equation const using enthalpy of TOV star and other fields */
  find_Euler_eq_const_TOV(grid);

  /* freeing */
  free_Grid_Params_S(GridParams);

  TOV_free(tov);
  
  return grid;
}

/* find Euler equation const using enthalpy of TOV star and 
// other fields that are found in TOV approximation.
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
static void find_Euler_eq_const_TOV(Grid_T *const grid)
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
static void sns_update_Aij(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ...\n");
  unsigned p;

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    sns_update_psi10A_UiUj(patch);
  }
  
  printf("Updating _A^{ij}, _dA^{ij} and _A^{ij}*A_{ij} ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* initialize the fields using TOV solution. */
static void init_field_TOV(Grid_T *const grid,const TOV_T *const tov)
{
  pr_line_custom('=');
  printf("Initializing the fields using TOV solution ...\n");

  const double M_NS = tov->ADM_m;/* NS adm mass */
  const double C_NS = GetParameterD_E("NS_center_-y_axis");/* center of NS it's on -y axis*/
  const double R_Schwar = tov->r[tov->N-1];/* NS's Schwarzchild radius */
  const double Omega_NS_x = GetParameterD_E("NS_Omega_U0");
  const double Omega_NS_y = GetParameterD_E("NS_Omega_U1");
  const double Omega_NS_z = GetParameterD_E("NS_Omega_U2");
  unsigned p;
  
  add_parameter_double("NS_Center",C_NS);
  
  /* initialization psi, eta and matter fields: */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    unsigned nn = patch->nn;
    unsigned ijk;
    
    PREP_FIELD(psi)
    PREP_FIELD(eta)
    PREP_FIELD(ETA)
    
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
        alpha = sqrt(1-2*M_NS/R_Schwar)/enthalpy_h; 
        ETA[ijk] = eta[ijk] = psi[ijk]*alpha;
        
        /* enthalpy */
        enthalpy[ijk] = enthalpy_h;
        
        /* rho0 */
        eos->h = enthalpy_h;
        rho0[ijk] = eos->rest_mass_density(eos);
        
        /* phi Newtonian approximation */
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
        
        /* eta */
        alpha = (1-0.5*M_NS/r)/(1+0.5*M_NS/r);
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
       
       /* Beta */
       Beta_U0[ijk] = 0;
       Beta_U1[ijk] = 0;
       Beta_U2[ijk] = 0;
       
       /* B1 */
       B1_U0[ijk] = 0;
       B1_U1[ijk] = 0;
       B1_U2[ijk] = 0;
       
       /* B0 */
       B0_U0[ijk] = Beta_U0[ijk]-B1_U0[ijk];
       B0_U1[ijk] = Beta_U1[ijk]-B1_U1[ijk];
       B0_U2[ijk] = Beta_U2[ijk]-B1_U2[ijk];
   }
      
  }/* end of * initializing Beta and B */
  
  printf("Initializing the fields using TOV solution ==> Done.\n");
  pr_clock();
  pr_line_custom('=');

}

/* given the max of radius of NS and its ceneter,
// create a grid with these properties.
// NOTE: WE assume we are using cubed spherical grid.
// ->return value: grid for single NS. */
static Grid_T *creat_sns_grid_CS(struct Grid_Params_S *const GridParams)
{
  Grid_T *grid = alloc_grid();/* adding a new grid */
  /* calculate the characteristics of this grid */
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
  const unsigned gn = grid->gn;
  const double NS_C = GetParameterD_E("NS_center_-y_axis");
  const double C    = -2*NS_C;
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
  if (!strcmp_i(kind,"SNS_CubedSpherical+Box_grid"))
    abortEr("This function only works with cubed spherical + box grid.\n");
    
  grid->kind = dup_s(kind);
  
  assert(GRT(C,0));
  assert(GRT(Max_R_NS_l,0));
  assert(LSS(2*Max_R_NS_l,C));
  
  /* making NS surface function */
  NS_surface_CubedSpherical_grid(grid,GridParams);
  
  box_size_l = GetParameterD_E("left_central_box_length_ratio")*Max_R_NS_l;
  
  for (i = 0; i < N_Outermost_Split; i++)
  {
    sprintf(var,"Outermost%u_radius",i);
    R_outermost[i] = GetParameterD_E(var);
    
    if (LSS(R_outermost[i],2*C))
      abortEr("the radius of outermost patches must be greater than twice of NS center.");
    
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
  
  /* right box */
  
  nlb[0] = (unsigned)GetParameterI("n_a");
  nlb[1] = (unsigned)GetParameterI("n_b");
  nlb[2] = (unsigned)GetParameterI("n_c");
  /* check for override */
  n = (unsigned)GetParameterI("right_box_n_a");
  if (n != INT_MAX)   nlb[0] = n;
  n = (unsigned)GetParameterI("right_box_n_b");
  if (n != INT_MAX)   nlb[1] = n;
  n = (unsigned)GetParameterI("right_box_n_c");
  if (n != INT_MAX)   nlb[2] = n;
    
  if(nlb[0] == INT_MAX)
    abortEr("n_a could not be set.\n");
  if(nlb[1] == INT_MAX)
    abortEr("n_b could not be set.\n");
  if(nlb[2] == INT_MAX)
    abortEr("n_c could not be set.\n");
  
  sprintf(par,"grid%u_right_box_n_a",gn);
  sprintf(val,"%u",nlb[0]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_box_n_b",gn);
  sprintf(val,"%u",nlb[1]);
  add_parameter_string(par,val);
  
  sprintf(par,"grid%u_right_box_n_c",gn);
  sprintf(val,"%u",nlb[2]);
  add_parameter_string(par,val);
  
  /* size a,b,c */
  sprintf(par,"grid%u_left_centeral_box_size_a",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_b",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_left_centeral_box_size_c",gn);
  add_parameter_double(par,box_size_l);
  
  sprintf(par,"grid%u_right_box_size_a",gn);
  add_parameter_double(par,C);
  
  sprintf(par,"grid%u_right_box_size_b",gn);
  add_parameter_double(par,C);
  
  sprintf(par,"grid%u_right_box_size_c",gn);
  add_parameter_double(par,C);
  
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
  
  /* assuming the center of right box at (0,C/2,0) */
  sprintf(par,"grid%u_right_box_center_a",gn);
  add_parameter_double(par,0.0);
  
  sprintf(par,"grid%u_right_box_center_b",gn);
  add_parameter_double(par,C/2);
  
  sprintf(par,"grid%u_right_box_center_c",gn);
  add_parameter_double(par,0.0);
  
  free(R_outermost);

  
  make_patches(grid);/* making patch(es) to cover the grid */
  realize_geometry(grid);/* realizing the geometry of whole grid
                     // including the way patches have been sewed,
                     // normal to the boundary, 
                     // outer-boundary, inner boundary and etc. */
  
  return grid;
  
}

/* making NS surfaces function */
static void NS_surface_CubedSpherical_grid(Grid_T *const grid,struct Grid_Params_S *const GridParams)
{
  const double Max_R_NS_l = GridParams->Max_R_NS_l;/* maximum radius of NS */
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
  if (strcmp_i(GridParams->NS_R_type,"PerfectSphere"))
  {
    for (i = 0; i < N[0]; ++i)
      for (j = 0; j < N[1]; ++j)
        for (k = 0; k < N[2]; ++k)
          R[L(N,i,j,k)] = Max_R_NS_l;
    
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
  }
  /* if NS radius is varied and we know its expansion in Ylm bases */
  else if (strcmp_i(GridParams->NS_R_type,"SphericalHarmonic"))
  {
    /* we need interpolation */
    const double *const realClm = GridParams->NS_R_Ylm->realClm;
    const double *const imagClm = GridParams->NS_R_Ylm->imagClm;
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
  }
  else
    abortEr(NO_OPTION);
  
  free(R);
}

/* given (X,Y,Z) in the specified slice of NS in cubed spherical coords
// it finds the associated polar and azimuthal angels on the surface of NS */
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

/* initialize Grid_Params struct */
static struct Grid_Params_S *init_GridParams(void)
{
  struct Grid_Params_S *par = calloc(1,sizeof(*par));
  return par;
}

/* free Grid_Params struct */
static void free_Grid_Params_S(struct Grid_Params_S *par)
{
  _free(par->NS_R_Ylm->realClm);
  _free(par->NS_R_Ylm->imagClm);
  _free(par);
}

/* find the NS center using d(enthalpy)/dx^i = 0 */
static void find_NS_center(Grid_T *const grid)
{
  double *NS_center;
  Root_Finder_T *root = init_root_finder(3);
  struct NC_Center_RootFinder_S params[1];
  double guess[3];/* initial guess for root finder */
  char par_name[1000];
  unsigned p;
  
  guess[0] = guess[2] = 0;
  guess[1] = GetParameterD_E("NS_center_-y_axis");
  params->root_finder = root;
  root->description = "Finding NS center:";
  root->type        = GetParameterS_E("RootFinder_Method");
  root->tolerance   = GetParameterD_E("RootFinder_Tolerance");
  root->MaxIter     = (unsigned)GetParameterI_E("RootFinder_Max_Number_of_Iteration");
  root->x_gss       = guess;
  root->params      = params;
  root->verbose     = 1;
  root->f[0]        = dh_dx0_root_finder_eq;
  root->f[1]        = dh_dx1_root_finder_eq;
  root->f[2]        = dh_dx2_root_finder_eq;    
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSPatch(patch))
      continue;

    root->interrupt = 0;
    params->patch = patch;
    plan_root_finder(root);
    NS_center     = execute_root_finder(root);
    
    /* if root finder was successful */
    if (root->residual < 1E-6)
    {
      //test
      printf("(%f,%f,%f)\n",NS_center[0],NS_center[1],NS_center[2]);
      //end
      
      /* save the position of NS center */
      sprintf(par_name,"grid%u_NS_center_x",grid->gn);
      add_parameter_array(par_name,NS_center,3);
      
      /* save the patch stem where the NS center takes place */
      sprintf(par_name,"grid%u_NS_center_patch",grid->gn);
      char *stem = strstr(patch->name,"_");
      assert(stem);
      stem++;
      add_parameter(par_name,stem);
      free(NS_center);
      break;
    }
    free(NS_center);
  }
  if (root->exit_status != ROOT_FINDER_OK)
  {
    print_root_finder_exit_status(root);
    abortEr("NS center could not be found.\n");
  }
  free_root_finder(root);
}

/* dh/dx^0 = 0 */
static double dh_dx0_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T *const patch = par->patch;
  DECLARE_FIELD(denthalpy_D0);
  Interpolation_T *interp_s;
  Needle_T *needle = alloc_needle();
  unsigned flg;
  double interp,X[3];
  
  needle->grid = patch->grid;
  needle->x    = x;
  needle_in(needle,patch);
  point_finder(needle);
  flg = needle->Nans;
  free_needle(needle);
  
  /* if this point is out of this patch, exit */
  if (flg == 0)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
  X_of_x(X,x,patch);
  interp_s = init_interpolation();
  interp_s->field = denthalpy_D0;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* dh/dx^1 = 0 */
static double dh_dx1_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T *const patch = par->patch;
  DECLARE_FIELD(denthalpy_D1);
  Interpolation_T *interp_s;
  Needle_T *needle = alloc_needle();
  unsigned flg;
  double interp,X[3];
  
  needle->grid = patch->grid;
  needle->x    = x;
  needle_in(needle,patch);
  point_finder(needle);
  flg = needle->Nans;
  free_needle(needle);
  
  /* if this point is out of this patch, exit */
  if (flg == 0)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
  X_of_x(X,x,patch);
  interp_s = init_interpolation();
  interp_s->field = denthalpy_D1;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* dh/dx^2 = 0 */
static double dh_dx2_root_finder_eq(void *params,const double *const x)
{
  struct NC_Center_RootFinder_S *const par = params;
  Patch_T *const patch = par->patch;
  DECLARE_FIELD(denthalpy_D2);
  Interpolation_T *interp_s;
  Needle_T *needle = alloc_needle();
  unsigned flg;
  double interp,X[3];
  
  needle->grid = patch->grid;
  needle->x    = x;
  needle_in(needle,patch);
  point_finder(needle);
  flg = needle->Nans;
  free_needle(needle);
  
  /* if this point is out of this patch, exit */
  if (flg == 0)
  {
    par->root_finder->interrupt = 1;
    return 0;
  }
  
  X_of_x(X,x,patch);
  interp_s = init_interpolation();
  interp_s->field = denthalpy_D2;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* adjust the center of NS at the designated point, in case it moved.
// we need only to draw enthalpy to (0,NS_C,0).
// to do so, we demand shifted_enthalpy(r) = enthalpy(R+r), in which R 
// is the amount the center is displaced from (0,-D/2,0);
// and finally update the enthalpy and its derivatives. */
static void adjust_NS_center(Grid_T *const grid)
{
  char par_name[1000];
  double *NS_center = 0;
  const double C    = GetParameterD_E("NS_center_-y_axis");
  sprintf(par_name,"grid%u_NS_center_x",grid->gn);
  NS_center = GetParameterArrayF_E(par_name);
  const double R[3] = {NS_center[0],NS_center[1]-C,NS_center[2]};
  unsigned p,ijk;
  
  /* if it is already located at the designted point */
  if (EQL(0,R[0]) && EQL(0,R[1]) && EQL(0,R[2]))
    return;
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (!IsItNSPatch(patch))
      continue;
    
    /* now shift enthalpy */
    Patch_T *patchp = 0;
    char *stem, hint[1000];
    DECLARE_FIELD(enthalpy);
    ADD_FIELD(shifted_enthalpy);
    DECLARE_FIELD(shifted_enthalpy);
    
    make_coeffs_3d(enthalpy);
    
    stem = strstr(patch->name,"_");/* the patch->name convention is grid\d?_root */
    assert(stem);
    stem++;
    sprintf(hint,"%s",stem);
    unsigned nn = patch->nn;
    double x[3],Xp[3];
    for (ijk = 0; ijk < nn; ++ijk)
    {
      x[0] = patch->node[ijk]->x[0]+R[0];
      x[1] = patch->node[ijk]->x[1]+R[1];
      x[2] = patch->node[ijk]->x[2]+R[2];
      find_Xp_and_patchp(x,hint,grid,Xp,&patchp);
      shifted_enthalpy->v[ijk] = interpolate_from_patch_prim("enthalpy",Xp,patchp);
    }
    /* now clean enthalpy and copy new value and remove extras */
    free_coeffs(enthalpy);
    free(enthalpy->v);
    enthalpy->v = shifted_enthalpy->v;
    shifted_enthalpy->v = 0;
    REMOVE_FIELD(shifted_enthalpy);
  }
  
  /* update enthalpy derivatives */
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    if(!IsItNSPatch(patch))
      continue;
      
    sns_update_derivative_enthalpy(patch);  
  }
  
  /* check if it works */
  if(1)
  {
    struct NC_Center_RootFinder_S par[1] = {0};
    Root_Finder_T root_finder[1] = {0};
    
    par->patch = GetPatch("left_centeral_box",grid);
    par->root_finder = root_finder;
    
    printf("dh/dx(%g,%g,%g)|NS center = %g\n",
      dh_dx0_root_finder_eq(par,NS_center),
      NS_center[0],NS_center[1],NS_center[2]);
    printf("dh/dy(%g,%g,%g)|NS center = %g\n",
      dh_dx1_root_finder_eq(par,NS_center),
      NS_center[0],NS_center[1],NS_center[2]);
    printf("dh/dz(%g,%g,%g)|NS center = %g\n",
      dh_dx2_root_finder_eq(par,NS_center),
      NS_center[0],NS_center[1],NS_center[2]);
  }
}
