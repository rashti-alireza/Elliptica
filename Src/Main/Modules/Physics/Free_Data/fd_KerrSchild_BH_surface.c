#include "fd_KerrSchild_BH_surface.h"

/* ->: alloc and populate radius for surface of BH
// find EVENT HORIZON intersected with the slice of a Kerr-Schild BH
// with arbitrary boost and spin. In principle it is different with AH.
// NOTE: R measured from the center of BH. */
double *fd_find_KerrSchild_BH_surface(Physics_T *const phys)
{
  FUNC_TIC
  
  const double M = Getd("Christodoulou_mass");
  const double a = Getd("spin_a");
  const double R = M+sqrt(Pow2(M)+Pow2(a));/* R non-rotated and boosted */
  struct BH_surface_RootFinder_S par[1] = {0};
  const Uint lmax   = (Uint)Geti("surface_Ylm_max_l");
  const Uint Ntheta = Ntheta_Ylm(lmax);
  const Uint Nphi   = Nphi_Ylm(lmax);
  const Uint Ntot   = Ntotal_Ylm(lmax);
  double *const R_BH = alloc_double(Ntot);/* R surface */
  double *R_res = alloc_double(Ntot);
  Uint i,j;
  
  /* set KS parameters */
  fd_KerrSchild_set_params(phys);
  
  /* populate root finder */
  Root_Finder_T *const root = init_root_finder(1);
  root->type      = "Steepest_Descent";
  root->tolerance = Getd("RootFinder_Tolerance");
  root->MaxIter   = (Uint)Geti("RootFinder_Iteration");
  root->params    = par;
  root->f[0]      = BH_surface_root_finder_eq;
  root->x_gss     = &R;
  root->verbose   = 0;//strstr_i(Gets("RootFinder_verbose"),"yes");
  plan_root_finder(root);
  
  /* initialize tables */
  init_Legendre_root_function();
  
  /* R surface in non-rotated and non-boosted coords. */
  par->R  = R;
  
  /* for each points of Ylm find the surface of BH */
  for (i = 0; i < Ntheta; ++i)
  {
    par->th = acos(-Legendre_root_function(i,Ntheta));
    for (j = 0; j < Nphi; ++j)
    {
      par->ph = j*2*M_PI/Nphi;
      
      execute_root_finder(root);
      
      R_BH[IJ_Ylm(i,j,Nphi)]  = root->x_sol[0];
      R_res[IJ_Ylm(i,j,Nphi)] = root->residual;
    }/* end of for (j = 0; j < Nphi; ++j) */
  }/* end of for (i = 0; i < Ntheta; ++i) */
  
  printf(Pretty0"L2 norm of BH surface finder = %e\n",
         L2_norm(Ntot,R_res,0));
  
  Free(R_res);
  free_root_finder(root);
  
  FUNC_TOC
  return R_BH;
}

/* find r such that f(r)-R = 0, where R is the radius of BH 
// in inertial and non rotated coordinates. */
static double BH_surface_root_finder_eq(void *params,const double *const r)
{
  const struct BH_surface_RootFinder_S *const pars = params;
  const double R  = pars->R;
  const double th = pars->th;
  const double ph = pars->ph;
  const double x  = r[0]*sin(th)*cos(ph);
  const double y  = r[0]*sin(th)*sin(ph);
  const double z  = r[0]*cos(th);
  KS_set_args
  
  return (R-farg->R);
}
