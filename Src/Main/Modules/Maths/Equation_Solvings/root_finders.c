/*
// Alireza Rashti
// September 2019
*/

#include "root_finders.h"

/* initializing a Root_Finder_T struct with calloc.
// n is the number of equations or equivalently number of variables
// ->return value: a pristine struct */
Root_Finder_T *init_root_finder(const unsigned n)
{
  Root_Finder_T *root = calloc(1,sizeof(*root));
  pointerEr(root);
  
  root->f = calloc(n+1,sizeof(*root->f));
  pointerEr(root->f);
  root->f[n] = 0;
  
  root->df_dx = calloc(n+1,sizeof(*root->df_dx));
  pointerEr(root->df_dx);
  root->df_dx[n] = 0;
  
  root->n = n;
  
  return root;
}

/* ->return value: root of system of equations {f0 = 0, f1 = 0, ...} */
double *execute_root_finder(Root_Finder_T *const root)
{
  return root->root_finder_func(root);
}

/* given the information it decides how to perform root finding */
void plan_root_finder(Root_Finder_T *const root)
{
  if (strcmp_i(root->type,"Steepest_Descent"))
  {
    root->root_finder_func = root_finder_steepest_descent;
  }
  else
    abortEr(NO_OPTION);
}

/* free the root_finder struct */
void free_root_finder(Root_Finder_T *root)
{
  if (!root)
    return;
    
  if (root->f)
   free(root->f);
   
  if (root->df_dx)
   free(root->df_dx);
   
  free(root);
}


/* using steepest descent method, find the root of function f(x)
// and then return the solution.
// note: it allocates memory for the solution.
//
// synopsis:
// =========
//
// Root_Finder_T *root = init_root_finder(n); # n is number of equations or (equivalently variables)
// root->type          = "Steepest_Descent";
// plan_root_finder(root);
// root->description = "solving f = 0";
// root->tolerance   = 10E-10;
// root->MaxIter     = 10;
// root->x_gss       = x0; # initial guess (OPTIONAL)
// root->params      = params; # parameters used for evaluation of fs
// root->f[0]        = f0; # f0 = 0 equation
// root->f[1]        = f1; # f1 = 0 equation
// # also if df_dx's are available (OPTIONAL):
// root->df_dx[0]    = df0_dx;
// root->df_dx[1]    = df1_dx;
// double *x_sol     = execute_root_finder(root);
// free_root_finder(root); # free struct root
// free(x_sol);
//
// ->return value: x solution that makes f(x) = 0 */
static double *root_finder_steepest_descent(Root_Finder_T *const root)
{
  const unsigned MaxIter = root->MaxIter;
  const unsigned n = root->n;
  void *params     = root->params;
  const double *const x_gss = root->x_gss;
  const double TOL = root->tolerance;
  double (**f)(void *params,const double *const x) = root->f;
  double (**df_dx)(void *params,const double *const x,const unsigned dir) = root->df_dx;
  double (*dg_dx)(void *params,double *const x,const unsigned dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const unsigned dir)) = 0;
  const char *const desc = root->description;
  double *const x = alloc_double(n);
  double z[n],y[n];
  double g0,g1,g2,g3,g,h1,h2,h3,alpha0,alpha,alpha2,alpha3,z0;
  Flag_T small_alpha3_flg = NO;
  unsigned i,k;
  
  if (desc)
    printf("%s:\n",desc);
  else
    printf("Finding root of {f(x) = 0}:\n");
  printf("Num. of Eqs. = %u, Tolerance = %e, Max. Num. of Iter. = %u\n",n,TOL,MaxIter);
  
  /* setup differentials */
  if (df_dx[0])/* if differentials are given */
    dg_dx = dg_dx_of_df_dx_SD;
  else
    dg_dx = dg_dx_FD_SD;
  
  if (x_gss)/* if no guess is given, x has been initialized to 0 */
    for (i = 0; i < n; ++i)
      x[i] = x_gss[i];
    
  k = 1;
  while (k <= MaxIter)
  {
    g1 = g_SD(f,params,x);
    printf(".. Step[%02u]: Residual{f(x) = 0} = %+e\n",k-1,sqrt(g1));
    
    for (i = 0; i < n; i++)
      z[i] = dg_dx(params,x,i,f,df_dx);
      
    z0 = L2_norm(n,z,0);
    if (EQL(z0,0.))
    {
      root->residual = sqrt(g1);
      printf("Root Finder -> Steepest Descent Method:\n"
             "Zero gradient thus an extrema; Residual = %e\n",root->residual);
      break;
    }
    
    for (i = 0; i < n; i++)
      z[i] /= z0;
    alpha3 = 1.;
    for (i = 0; i < n; i++)
      y[i] = x[i] - alpha3*z[i];
    g3 = g_SD(f,params,y);
    
    while (g3 >= g1)
    {
      alpha3 /= 2.;
      for (i = 0; i < n; i++)
        y[i] = x[i] - alpha3*z[i];
      g3 = g_SD(f,params,y);
      
      if(alpha3 < 0.5*TOL)
      {
        root->residual = sqrt(g3);
        printf("Root Finder -> Steepest Descent Method:\n"
             "No likely improvement; Residual = %e\n",root->residual);
        
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
    g2 = g_SD(f,params,y);
    h1 = (g2-g1)/alpha2;
    h2 = (g3-g2)/(alpha3-alpha2);
    h3 = (h2-h1)/alpha3;
    alpha0 = 0.5*(alpha2 - h1/h3);
    for (i = 0; i < n; i++)
      y[i] = x[i] - alpha0*z[i];
    g0 = g_SD(f,params,y);
    
    alpha = g0 < g3 ? alpha0 : alpha3;
    
    for (i = 0; i < n; i++)
      x[i] -= alpha*z[i];
    
    g = g0 < g3 ? g0 : g3;
    root->residual = sqrt(g);
    if (fabs(g-g1) < TOL)
    { 
      printf(".. Step[%02u]: Residual[f(x) = 0] = %+e\n",k,root->residual);
      printf("Root Finder -> Steepest Descent Method:\n"
             "The root(s) are found => Residual = %e\n",root->residual);
      break;
    }
    
    k++;
    if (k == MaxIter+1)
    {  
      printf(".. Step[%02u]: Residual[f(x) = 0] = %+e\n",k-1,root->residual);
      printf("Root Finder -> Steepest Descent Method:\n"
             "Exceeds maximum number of iterations => Residual = %e\n",root->residual);
      break;
    }
  }
  
  fflush(stdout);
  root->x_sol = x;
  return x;
}

/* ->return value: f0^2(params,x)+f1^2(params,x) + ... */
static double g_SD(double (**f)(void *params,const double *const x),void *params,const double *const x)
{
  unsigned i = 0;
  double g = 0;
  
  while (f[i])
  {
    g += pow(f[i](params,x),2);
    i++;
  }
  
  return g;
}

/* ->return value: derivative of e.g. f0^2(x0,x1)+f1^2(x0,x1) given df0_dx, df1_dx etc. 
// at point x with respect to x^{dir} using the given df_dx's */
static double dg_dx_of_df_dx_SD(void *params,double *const x,const unsigned dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const unsigned dir))
{
  unsigned i = 0;
  double dg_dx = 0;
  
  while (f[i])
  {
    dg_dx += f[i](params,x)*df_dx[i](params,x,dir);
    i++;
  }
  dg_dx *= 2.;
  
  return dg_dx;
}

/* ->return value: derivative of e.g. f0^2(x0,x1)+f1^2(x0,x1) 
// at point x with respect to x^{dir} using finite difference (three-point midpoint formula)*/
static double dg_dx_FD_SD(void *params,double *const x,const unsigned dir,double (**f)(void *params,const double *const x),double (**df_dx)(void *params,const double *const x,const unsigned dir))
{
  const double eps    = 10E-10;
  const double fabsx  = fabs(x[dir])*eps;
  const double h      = fabsx > eps ? fabsx : eps;/* just in case fabsx is 0, so h won't get 0 */
  double gl,gr,dg_dx;
  
  x[dir] -= h;/* x' = x-h */
  gl      = g_SD(f,params,x);/* g(x-h) */
  
  x[dir] += 2.*h;/* x' + 2h => x' = x+h */
  gr      = g_SD(f,params,x);/* g(x+h) */
  
  dg_dx = 0.5*(gr-gl)/h;
  
  x[dir] -= h;/* x'-h => x */
  
  return dg_dx;
  UNUSED(df_dx);
}
