/*
// Alireza Rashti
// September 2019
*/

#include "root_finders.h"

/*
// synopsis:
// =========
//
// "Steepest_Descent":
// ===================
// note: this function always finds the roots of the equations
// if they are smooth, if they be only C^0 the root might not be found.
// in this case we recommand bisect method.
//
// Root_Finder_T *root = init_root_finder(n); # n is number of equations or (equivalently variables)
// root->type          = "Steepest_Descent";
// plan_root_finder(root);
// root->description = "solving f = 0";
// root->verbose     = 1;# to prints every steps of solution
// root->tolerance   = 10E-10;
// root->MaxIter     = 10;
// root->x_gss       = x0; # initial guess (OPTIONAL)
// root->params      = params; # parameters used for evaluation of fs
// root->f[0]        = f0; # f0 = 0 equation
// root->f[1]        = f1; # f1 = 0 equation
// # also if df_dx's are available (OPTIONAL):
// # note if df_dx's not specified it automatically uses finite difference method
// root->df_dx[0]    = df0_dx;
// root->df_dx[1]    = df1_dx;
// # to force use right side ot left side finite difference:
// root->FD_Right = 1; # right side f.d. method
// # or
// root->FD_Left  = 1; # left side f.d. method
// execute_root_finder(root);
// double *x_sol  = root->x_sol;
// # some checks:
// * to check the exit status of root finder do: *
// print_root_finder_exit_status(root);
// * to find out the residual *
// printf("%e",root->residual);
// * note one can interrupt the root finder by setting 
//   root->interupt = 1 inside the params of equation. *
// root->x_sol = 0; ** if you DON'T want to free the solution
// free_root_finder(root); # free struct root
//
//
// "Bisect_Single":
// ================
// Root_Finder_T *root = init_root_finder(1); # only 1 equation
// root->type          = "Bisect_Single";
// plan_root_finder(root);
// root->description = "solving f = 0";
// root->verbose     = 1;# to prints every steps of solution
// root->tolerance   = 10E-10;
// root->MaxIter     = 1000;
// root->params      = params; # parameters used for evaluation of fs
// root->f[0]        = f0; # f0 = 0 equation
// root->a_bisect    = a; # f(x) must change sign for x in [a,b]
// root->b_bisect    = b; # f(x) must change sign for x in [a,b]
// execute_root_finder(root);
// double *x_sol = root->x_sol;
// # some checks:
// * to check the exit status of root finder do: *
// print_root_finder_exit_status(root);
// * to find out the residual *
// printf("%e",root->residual);
// * note one can interrupt the root finder by setting 
//   root->interupt = 1 inside the params of equation. *
// root->x_sol = 0; ** if you DON'T want to free the solution
// free_root_finder(root); # free struct root
*/

/* initializing a Root_Finder_T struct with calloc.
// n is the number of equations or equivalently number of variables
// ->return value: a pristine struct */
Root_Finder_T *init_root_finder(const Uint n)
{
  Root_Finder_T *root = calloc(1,sizeof(*root));
  IsNull(root);
  
  root->x_sol = calloc(n,sizeof(*root->f));
  IsNull(root->x_sol);
  
  root->f = calloc(n+1,sizeof(*root->f));
  IsNull(root->f);
  root->f[n] = 0;
  
  root->df_dx = calloc(n+1,sizeof(*root->df_dx));
  IsNull(root->df_dx);
  root->df_dx[n] = 0;
  
  root->n = n;
  root->eq_number = UINT_MAX;
  return root;
}

/* ->return value: root of system of equations {f0 = 0, f1 = 0, ...} */
double *execute_root_finder(Root_Finder_T *const root)
{
  if (!root->root_finder_func)
    Error0("No plan for the root finder.\n");
    
  return root->root_finder_func(root);
}

/* given the information it decides how to perform root finding */
void plan_root_finder(Root_Finder_T *const root)
{
  if (strcmp_i(root->type,"Steepest_Descent"))
  {
    root->root_finder_func = root_finder_steepest_descent;
  }
  else if (strcmp_i(root->type,"Bisect_Single"))
  {
    root->root_finder_func = root_finder_bisect_single;
  }
  
  else
    Error0(NO_OPTION);
}

/* free the root_finder struct */
void free_root_finder(Root_Finder_T *root)
{
  if (!root)
    return;
  
  Free(root->x_sol);
  Free(root->f);
  Free(root->df_dx);
  
  free(root);
}

/* using bisect method, find the root of function f(x)
// and then return the solution.
// note: this is only works for a single equation.
// note: it allocates memory for the solution.
// note: f(x) must change sign for x in [a,b].
// ->return value: x solution that makes f(x) = 0 */
static double *root_finder_bisect_single(Root_Finder_T *const root)
{
  const double tic = get_time_sec();
  const Uint MaxIter = root->MaxIter;
  const Uint n = root->n;
  void *params     = root->params;
  const double TOL = root->tolerance;
  double (**f)(void *params,const double *const x) = root->f;
  const double A = root->a_bisect;
  const double B = root->b_bisect;
  const char *const desc = root->description;
  double *const x = root->x_sol;
  double a = A,b = B;
  double FA,FP,p[1] = {0},d;
  Uint i;
  
  /* some checks */
  if (!f)
    Error0("\n~> No equation has been given.\n");
  if (n != 1)
    Error0("This bisect method only works for single equation.");
  if (f[0](params,&A)*f[0](params,&B) > 0)
    Error0("f(x) must change sign for x in [a,b].\n");
  
  if (desc)
    printf("%s\n",desc);
  if (root->verbose)
  {
    printf("\n{ Root finder ...\n\n");
    printf("|--> Root Finder  = Bisect Single Method\n");
    printf("|--> Num. of Eqs. = %u\n"
           "|--> Tolerance    = %e\n"
           "|--> Max. Iter.   = %u\n",
            n,TOL,MaxIter);
  }
  
  /* initialize */
  root->exit_status = ROOT_FINDER_UNDEF;
  root->residual    = DBL_MAX;
  root->interrupt   = 0;
  
  i = 1;
  FA = f[0](params,&a);
  while(i <= MaxIter)
  {
    d    = 0.5*(b-a);
    p[0] = a+d;
    FP   = f[0](params,p);
    root->residual = FP;
    
    /* some checks */
    if (root->interrupt != 0)
    {
      root->exit_status = ROOT_FINDER_INTERRUPTED;
      if (root->verbose)
      {
        printf("\n~> Root finder was interrupted!\n");
      }
      break;
    }
    if (!isfinite(root->residual))
    {
      root->exit_status = ROOT_FINDER_NAN;
      if (root->verbose)
      {
        printf("\n~> Residual is abnormal; Residual = %e\n",root->residual);
      }
      break;
    }
    /* verbose */
    if (root->verbose)
    {
      printf("|--> Step[%02u]: Residual{f(x) = 0} = %+e\n",i,root->residual);
      fflush(stdout);
    }
    
    if (EQL(FP,0) || LSS(fabs(d),TOL))
    {
      x[0] = p[0];
      if (root->verbose)
        printf("|--> Step[%02u]: Residual{f(x) = 0} = %+e\n",i,root->residual);
      root->exit_status = ROOT_FINDER_OK;
      break;
    }
    
    i++;
    
    if (FA*FP > 0)
    {
      a  = p[0];
      FA = FP;
    }
    else
    {
      b = p[0];
    }
    
  }/* end of while(i <= MaxIter) */
  if (i > MaxIter)
  { 
    root->exit_status = ROOT_FINDER_MAX_ITER; 
    if (root->verbose)
    {
      printf("|--> Step[%02u]: Residual[f(x) = 0] = %+e\n",i-1,root->residual);
      printf("\n~> Exceeds maximum number of iterations => Residual = %e\n",root->residual);
    }
  }
  
  root->x_sol = x;
  
  if (root->verbose)
    printf("\n} Root finder --> Done. ( Wall-Clock = %.2fs )\n\n",
           get_time_sec()-tic);
  fflush(stdout);
  
  return x;
}

/* using steepest descent method, find the root of function f(x)
// and then return the solution.
// note: it allocates memory for the solution.
// ->return value: x solution that makes f(x) = 0 */
static double *root_finder_steepest_descent(Root_Finder_T *const root)
{
  const double tic = get_time_sec();
  const Uint MaxIter = root->MaxIter;
  const Uint n = root->n;
  void *params = root->params;
  const double *const x_gss = root->x_gss;
  const double TOL = root->tolerance;
  double (**f)(void *params,const double *const x) = root->f;
  double (**df_dx)(void *params,const double *const x,const Uint dir) = root->df_dx;
  double (*dg_dx)(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root) = 0;
  const char *const desc = root->description;
  double *const x = root->x_sol;
  double z[n],y[n];
  double g0,g1,g2,g3,g,h1,h2,h3,alpha0,alpha,alpha2,alpha3,z0;
  double res;
  Flag_T small_alpha3_flg = NO;
  Uint i,k;
  
  if (!f)
    Error0("\n~> No equation has been given.\n");
  
  if (desc)
    printf("%s\n",desc);
  if (root->verbose)
  {
    printf("\n{ Root finder ...\n\n");
    printf("|--> Root Finder  = Steepst Descent Method\n");
    printf("|--> Num. of Eqs. = %u\n"
           "|--> Tolerance    = %e\n"
           "|--> Max. Iter.   = %u\n",
            n,TOL,MaxIter);
  }
  
  /* setup differentials */
  if (df_dx[0])/* if differentials are given */
    dg_dx = dg_dx_of_df_dx_SD;
  else if (root->FD_Left)/* left side stencil */
    dg_dx = dg_dx_FD3L_SD;
  else if (root->FD_Right)/* right side stencil */
    dg_dx = dg_dx_FD3R_SD;
  else/* if not specified */
    dg_dx = dg_dx_FD3M_SD;
  
  if (x_gss)/* if no guess is given, x has been initialized to 0 */
    for (i = 0; i < n; ++i)
      x[i] = x_gss[i];
  
  /* initialize */
  root->exit_status = ROOT_FINDER_UNDEF;
  root->residual    = DBL_MAX;
  root->interrupt   = 0;
  
  k = 1;
  while (k <= MaxIter)
  {
    g1 = g_SD(f,params,x,root);
    if (root->interrupt != 0)
    {
      root->exit_status = ROOT_FINDER_INTERRUPTED;
      if (root->verbose)
      {
        printf("\n~> Root finder was interrupted!\n");
      }
      break;
    }
    if (root->verbose)
    {
      res = sqrt(g1);
      if (!isfinite(res))
      {
        root->exit_status = ROOT_FINDER_NAN;
        break;
      }
      else
        printf("|--> Step[%02u]: Residual{f(x) = 0} = %+e\n",k-1,res);
    }
    
    for (i = 0; i < n; i++)
      z[i] = dg_dx(params,x,i,f,df_dx,root);
      
    z0 = root_square(n,z,0);
    if (EQL(z0,0.))
    {
      root->residual = sqrt(g1);
      if (!isfinite(root->residual))
        root->exit_status = ROOT_FINDER_NAN;
      else
        root->exit_status = ROOT_FINDER_EXTREMA;
        
      if (root->verbose)
      {
        if (root->exit_status == ROOT_FINDER_NAN)
        {
          printf("\n~> Residual is abnormal; Residual = %e\n",root->residual);
        }
        else
        {
          printf("\n~> Zero gradient, it hit an extrema; Residual = %e\n",root->residual);
        }
      }
      
      break;
    }
    
    for (i = 0; i < n; i++)
      z[i] /= z0;
    alpha3 = 1.;
    for (i = 0; i < n; i++)
      y[i] = x[i] - alpha3*z[i];
    g3 = g_SD(f,params,y,root);
    if (root->interrupt != 0)
    {
      root->exit_status = ROOT_FINDER_INTERRUPTED;
      if (root->verbose)
      {
        printf("\n~> Root finder was interrupted!\n");
      }
      break;
    }
    
    while (g3 >= g1)
    {
      alpha3 /= 2.;
      for (i = 0; i < n; i++)
        y[i] = x[i] - alpha3*z[i];
      g3 = g_SD(f,params,y,root);
      if (root->interrupt != 0)
      {
        root->exit_status = ROOT_FINDER_INTERRUPTED;
        if (root->verbose)
        {
          printf("\n~> Root finder was interrupted!\n");
        }
        break;
      }
      
      if(alpha3 < 0.5*TOL)
      {
        root->residual = sqrt(g3);
        if (!isfinite(root->residual))
          root->exit_status = ROOT_FINDER_NAN;
        else
          root->exit_status = ROOT_FINDER_NO_IMPROVEMENT;
        
        if (root->verbose)
        {
          if (root->exit_status == ROOT_FINDER_NAN)
          {
            printf("\n~> Residual is abnormal; Residual = %e\n",root->residual);
          }
          else
          {
            printf("\n~> No likely improvement; Residual = %e\n",root->residual);
          }
        }
        
        for (i = 0; i < n; i++)
          x[i] = y[i];
          
        small_alpha3_flg = YES;
        break;
      }
    }
    if (small_alpha3_flg == YES)
      break;
      
    alpha2 = alpha3*0.5;
    for (i = 0; i < n; i++)
      y[i] = x[i] - alpha2*z[i];
    g2 = g_SD(f,params,y,root);
    if (root->interrupt != 0)
    {
      root->exit_status = ROOT_FINDER_INTERRUPTED;
      if (root->verbose)
      {
        printf("\n~> Root finder was interrupted!\n");
      }
      break;
    }
    h1 = (g2-g1)/alpha2;
    h2 = (g3-g2)/(alpha3-alpha2);
    h3 = (h2-h1)/alpha3;
    alpha0 = 0.5*(alpha2 - h1/h3);
    for (i = 0; i < n; i++)
      y[i] = x[i] - alpha0*z[i];
    g0 = g_SD(f,params,y,root);
    if (root->interrupt != 0)
    {
      root->exit_status = ROOT_FINDER_INTERRUPTED;
      if (root->verbose)
      {
        printf("\n~> Root finder was interrupted!\n");
      }
      break;
    }
    
    alpha = g0 < g3 ? alpha0 : alpha3;
    
    for (i = 0; i < n; i++)
      x[i] -= alpha*z[i];
    
    g = g0 < g3 ? g0 : g3;
    root->residual = sqrt(g);
    if (!isfinite(root->residual))
    {
      root->exit_status = ROOT_FINDER_NAN;
      if (root->verbose)
      {
        printf("\n~> Residual is abnormal; Residual = %e\n",root->residual);
      }
      break;
    }
    if (fabs(g-g1) < TOL)
    { 
      root->exit_status = ROOT_FINDER_OK;
      
      if (root->verbose)
      {
        printf("|--> Step[%02u]: Residual[f(x) = 0] = %+e\n",k,root->residual);
        printf("\n~> The root(s) are found => Residual = %e\n",root->residual);
      }
      break;
    }
    
    k++;
    if (k == MaxIter+1)
    { 
      root->exit_status = ROOT_FINDER_MAX_ITER; 
      if (root->verbose)
      {
        printf("|--> Step[%02u]: Residual[f(x) = 0] = %+e\n",k-1,root->residual);
        printf("\n~> Exceeds maximum number of iterations => Residual = %e\n",root->residual);
      }
      break;
    }
  }
  
  fflush(stdout);
  root->x_sol = x;
  
  if (root->verbose)
    printf("\n} Root finder --> Done. ( Wall-Clock = %.2fs )\n\n",
           get_time_sec()-tic);
  
  return x;
}

/* printing the status of root finder */
void print_root_finder_exit_status(const Root_Finder_T *const root)
{
  if (root->description)
    printf("%s\n",root->description);
    
  switch(root->exit_status)
  {
    case ROOT_FINDER_OK:
      printf(Pretty0"Root finder found the root(s) up to the specified tolerance.\n");
      printf(Pretty0"Root finder residual = %e.\n",root->residual);
    break;
    case ROOT_FINDER_EXTREMA:
      printf(Pretty0"Root finder hit an extrema.\n");
      printf(Pretty0"Root finder residual = %e.\n",root->residual);
    break;
    case ROOT_FINDER_MAX_ITER:
      printf(Pretty0"Root finder exceeded the maximum number of iteration.\n");
      printf(Pretty0"Root finder residual = %e.\n",root->residual);
    break;
    case ROOT_FINDER_NO_IMPROVEMENT:
      printf(Pretty0"Root finder cannot improve the solution further.\n");
      printf(Pretty0"Root finder residual = %e.\n",root->residual);
    break;
    case ROOT_FINDER_INTERRUPTED:
      printf(Pretty0"Root finder was interrupted.\n");
    break;
    case ROOT_FINDER_NAN:
      printf(Pretty0"Root finder failed with an abnormal residual.\n");
    break;
    default:
      printf(Pretty0"The status is not defined.\n");
  }
  fflush(stdout);
}

/* ->return value: f0^2(params,x)+f1^2(params,x) + ... */
static double g_SD(double (**f)(void *params,const double *const x),void *params,const double *const x,Root_Finder_T *const root)
{
  Uint i = 0;
  double g = 0;
  
  while (f[i])
  {
    root->eq_number = i;
    g += pow(f[i](params,x),2);
    i++;
  }
  root->eq_number = UINT_MAX;/* catch error */
  return g;
}

/* ->return value: derivative of e.g. f0^2(x0,x1)+f1^2(x0,x1) given df0_dx, df1_dx etc. 
// at point x with respect to x^{dir} using the given df_dx's */
static double dg_dx_of_df_dx_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root)
{
  Uint i = 0;
  double dg_dx = 0;
  
  while (f[i])
  {
    root->eq_number = i;
    dg_dx += f[i](params,x)*df_dx[i](params,x,dir);
    i++;
  }
  dg_dx *= 2.;
  root->eq_number = UINT_MAX;/* catch error */
  return dg_dx;
}

/* ->return value: derivative of e.g. f0^2(x0,x1)+f1^2(x0,x1) 
// at point x with respect to x^{dir} using finite difference (three-point midpoint formula)*/
static double dg_dx_FD3M_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root)
{
  const double eps    = 10E-5;
  const double fabsx  = fabs(x[dir])*eps;
  const double h      = fabsx > eps ? fabsx : eps;/* just in case fabsx is 0, so h won't get 0 */
  double gl,gr,dg_dx;
  
  x[dir] -= h;/* x' = x-h */
  gl      = g_SD(f,params,x,root);/* g(x-h) */
  
  x[dir] += 2.*h;/* x' + 2h => x'' = x+h */
  gr      = g_SD(f,params,x,root);/* g(x+h) */
  
  dg_dx = 0.5*(gr-gl)/h;
  
  x[dir] -= h;/* x''-h => x */
  
  return dg_dx;
  UNUSED(df_dx);
}

/* ->return value: derivative of e.g. f0^2(x0,x1)+f1^2(x0,x1)
// at point x with respect to x^{dir} using finite difference (three-point right end formula)*/
static double dg_dx_FD3R_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root)
{
  const double eps    = 10E-5;
  const double fabsx  = fabs(x[dir])*eps;
  const double h      = fabsx > eps ? fabsx : eps;/* just in case fabsx is 0, so h won't get 0 */
  double g,gr1,gr2,dg_dx;
  
  g       = g_SD(f,params,x,root);/* g(x) */
  
  x[dir] += h;/* x' = x+h */
  gr1     = g_SD(f,params,x,root);/* g(x+h) */
  
  x[dir] += h;/* x' + h => x'' = x+2h */
  gr2     = g_SD(f,params,x,root);/* g(x+2h) */
  
  dg_dx = 0.5*(-3*g+4*gr1-gr2)/h;
  
  x[dir] -= 2.*h;/* x''-2h => x */
  
  return dg_dx;
  UNUSED(df_dx);
}

/* ->return value: derivative of e.g. f0^2(x0,x1)+f1^2(x0,x1)
// at point x with respect to x^{dir} using finite difference (three-point left end formula)*/
static double dg_dx_FD3L_SD(void *params,double *const x,const Uint dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const Uint dir),Root_Finder_T *const root)
{
  const double eps    = 10E-5;
  const double fabsx  = fabs(x[dir])*eps;
  const double h      = fabsx > eps ? fabsx : eps;/* just in case fabsx is 0, so h won't get 0 */
  double g,gl1,gl2,dg_dx;
  
  g       = g_SD(f,params,x,root);/* g(x) */
  
  x[dir] -= h;/* x' = x-h */
  gl1     = g_SD(f,params,x,root);/* g(x-h) */
  
  x[dir] -= h;/* x' - h => x'' = x-2h */
  gl2     = g_SD(f,params,x,root);/* g(x-2h) */
  
  dg_dx = 0.5*(3*g-4*gl1+gl2)/h;
  
  x[dir] += 2.*h;/* x''+2h => x */
  
  return dg_dx;
  UNUSED(df_dx);
}

