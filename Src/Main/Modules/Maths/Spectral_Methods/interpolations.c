/*
// Alireza Rashti
// August 2018
*/

#include "interpolations.h"

/* synopsis:
// =========
//
// ** filling the interpolation struct **
// Interpolation_T *interp_s = init_interpolation();
// interp_s->field = alpha;
// interp_s->X = p1;
// interp_s->Y = p2;
// * NOTE: X and Y and Z must be given in coords which patch uses *
// interp_s->Z_dir_flag = 1; * means it only interpolates along Z direction *
// interp_s->I = i0; * since the interpolation takes place for const i *
// interp_s->J = j0; * since the interpolation takes place for const j *
// ** planning the appropriate function for interpolation **
// plan_interpolation(interp_s);
//
// ** evaluating interpolation **
// value = execute_interpolation(interp_s);
// 
// ** freeing **
// free_interpolation(interp_s);
*/

/* initializing an Interpolation_T struct with calloc */
Interpolation_T *init_interpolation(void)
{
  Interpolation_T *s = calloc(1,sizeof(*s));
  IsNull(s);
  
  return s;
}

/* interpolation function */
double execute_interpolation(Interpolation_T *const interp_struct)
{
  return interp_struct->interpolation_func(interp_struct);
}

/* interpolation function for derivatives */
double execute_derivative_interpolation(Interpolation_T *const interp_struct)
{
  return interp_struct->interpolation_derivative_func(interp_struct);
}

// Assigns pointer such as interp_s->f so that spline interpolation
// functions can be used irrespective of spline method.
void assign_interpolation_ptrs(Interpolation_T *const interp_s)
{
  if (strstr_i(interp_s->method,"Natural_Cubic_Spline_1D"))
  {
    interp_s->f = &interp_s->N_cubic_spline_1d->f;
    interp_s->x = &interp_s->N_cubic_spline_1d->x;
    interp_s->h = &interp_s->N_cubic_spline_1d->h;
    interp_s->N = &interp_s->N_cubic_spline_1d->N;
  }
  else if (strstr_i(interp_s->method,"Hermite_Cubic_Spline") || strstr_i(interp_s->method,"Hermite_Spline")
          || strstr_i(interp_s->method,"Hermite"))
  {
    interp_s->f = &interp_s->Hermite_spline_1d->f;
    interp_s->x = &interp_s->Hermite_spline_1d->x;
    interp_s->h = &interp_s->Hermite_spline_1d->h;
    interp_s->N = &interp_s->Hermite_spline_1d->N;
  }
}

/* planning interpolation function based on 
// input file and patch and flags in interp_s. */
void plan_interpolation(Interpolation_T *const interp_s)
{
  Uint i;
  Patch_T *patch;
  /* choose method of interpolation based on inputs and patch
  // and narrow down the options. */
  if (strstr_i(interp_s->method,"Neville_1D"))
  {
    interp_s->interpolation_func = interpolation_Neville_1d;
  }
  else if (strstr_i(interp_s->method,"Natural_Cubic_Spline_1D"))
  {
    order_arrays_natural_cubic_spline_1d(interp_s);
    find_coeffs_natural_cubic_spline_1d(interp_s);
    interp_s->interpolation_func = interpolation_natural_cubic_spline_1d;
    interp_s->interpolation_derivative_func = interpolation_NCS_derivative;
  }
  else if (strstr_i(interp_s->method,"Hermite"))
  {
    order_arrays_spline_1d(interp_s);
    
    if (!&interp_s->finite_diff_order) { interp_s->finite_diff_order = (Uint)Pgeti("Interpolation_finite_diff_order"); }
    if (!&interp_s->Spline_Order) { interp_s->Spline_Order = (Uint)Pgeti("Interpolation_spline_order"); }
    find_coeffs_Hermite_spline(interp_s);
    interp_s->interpolation_func = interpolation_Hermite_spline;
    interp_s->interpolation_derivative_func = interpolation_finite_difference;
  }
  else if (strstr_i(interp_s->method,"Spectral") || 
       strstr_i(PgetsEZ("Interpolation_Method"),"Spectral"))
  {
    fPick_Func_T *func = 0;
    if (!interp_s->field)
      Error0("Bad argument:No field specified for the interpolation.\n");
    patch = interp_s->field->patch;
  
    for (i = 0; i < 3; ++i)
    {
      if (patch->basis[i]       == Chebyshev_Tn_BASIS &&
          patch->collocation[i] == Chebyshev_Extrema    )
      {
        func = interpolation_Chebyshev_Tn;
      }
      else
         Error0(INCOMPLETE_FUNC);

    }/* end of for (i = 0; i < 3; ++i) */
    
    /* having chosen the type of interpolation, 
    // find out in which direction(s) interpolation happens
    // and then assign which function should be called.
    */
    interp_s->interpolation_func = func(interp_s);
  }
  else
    Error0(INCOMPLETE_FUNC);
    
  // If finite difference derivative is manually assigned
  if (PgetsEZ("interpolation_derivative_method") 
      && strstr_i(PgetsEZ("interpolation_derivative_method"),
                    "finite_difference"))
  {
    interp_s->interpolation_derivative_func = interpolation_finite_difference;
    interp_s->finite_diff_order = (Uint)Pgeti("Interpolation_finite_diff_order");
  }
  if (PgetsEZ("interpolation_derivative_method")
      && strstr_i(PgetsEZ("interpolation_derivative_method"), "manual_FDM"))
  {
    interp_s->finite_diff_order = (Uint)Pgeti("Interpolation_finite_diff_order");
    interp_s->interpolation_derivative_func = interpolation_manual_FDM;
  }
}

/* xi's must be in increasing order, make sure this happens. */
static void order_arrays_natural_cubic_spline_1d(Interpolation_T *const interp_s)
{
  double *x = interp_s->N_cubic_spline_1d->x;
  double *f = interp_s->N_cubic_spline_1d->f;
  double *y;/* ordered x */
  double *g;/* ordered f */
  const Uint N = interp_s->N_cubic_spline_1d->N;
  Uint i;
  
  if (EQL(x[0],x[N-1]))
  {
    Error0("Periodic case has not been considered.\n");
  }
  if (GRT(x[0],x[N-1]))/* if the x's are in decreasing order */
  {
    interp_s->N_cubic_spline_1d->Alloc_Mem = 1;
    y = alloc_double(N);
    g = alloc_double(N);
    
    for (i = 0; i < N; ++i)
    {
      y[i] = x[N-1-i];
      g[i] = f[N-1-i];
    }
    interp_s->N_cubic_spline_1d->x = y;
    interp_s->N_cubic_spline_1d->f = g;
  }
  interp_s->N_cubic_spline_1d->Order = 1;
}

static void order_arrays_spline_1d(Interpolation_T *const interp_s)
{
  double *x = *interp_s->x;
  double *f = *interp_s->f;
  double *y;/* ordered x */
  double *g;/* ordered f */
  const Uint N = *interp_s->N;
  Uint i;
  
  if (EQL(x[0],x[N-1]))
  {
    Error0("Spline arrays: Periodic case has not been considered.\n");
  }
  if (GRT(x[0],x[N-1]))/* if the x's are in decreasing order */
  {
    y = alloc_double(N);
    g = alloc_double(N);
    
    for (i = 0; i < N; ++i)
    {
      y[i] = x[N-1-i];
      g[i] = f[N-1-i];
    }
    *interp_s->x = y;
    *interp_s->f = g;
  }
  set_interp_order_flag(interp_s, 1);
  set_interp_alloc_mem_flag(interp_s, 1);
}

///////////////////////////////////////Finite difference approximation////////////////////////
// Wrapper for FDM_Fornberg in derivatives.c (Fornberg finite difference method)
static double interpolation_finite_difference(Interpolation_T *const interp_s)
{
  const double *const x = *interp_s->x;
  const double *const f = *interp_s->f;
  const double h = *interp_s->h;
  const Uint N = *interp_s->N;
  const Uint n = interp_s->finite_diff_order;
  const Uint M = interp_s->FDM_derivative;
  
  return FDM_Fornberg(x, f, h, M, n, N);
}

/*
static double interpolation_finite_difference(Interpolation_T *const interp_s)
{
  // Approximates M-th derivative of f(x)|x=h by finite difference method.
  const double *const x = *interp_s->x;
  double *const f = *interp_s->f;
  const double h = *interp_s->h;
  const Uint N = *interp_s->N;
  Uint i = 0;
  const Uint n = interp_s->finite_diff_order;
  const Uint M = interp_s->FDM_derivative;
  double ret = DBL_MAX;/
  Flag_T flg = NONE;
  
  // Checks if we have enough data points for given order.
  if (N <= n+M)
  {
    printf("Points: %i\n", N);
    printf("Order of accuracy: %i\n", n);
    Error0("Finite difference error: Not enough points for desired accuracy.\n");
  }
  if (M >= n)
  {
    printf("Degree of derivative: %i\n", M);
    printf("Order of accuracy: %i\n", n);
    Error0("Finite difference error: Degree of derivative must not exceed degree of accuracy.\n");
  }
  
  // Finds the data segment
  for (i = 0; i < N-1; ++i)
  {
    if (GRTEQL(h,x[i]) && LSSEQL(h,x[i+1]))
    {
      flg = FOUND;
      break;
    }
  }
  
  if (flg != FOUND)
  {
    if (!interp_s->N_cubic_spline_1d->No_Warn)
      Warning("The given point for the interpolation is out of the domain.\n");
    
    return ret;
  }
  
  // If we have enough points on either side of the desired point h,
  // we can use a central finite difference. Otherwise, we must
  // use a shifted finite difference (e.g. the forward finite difference
  // if x[0] < = h < x[1]).
  // Note: If the desired accuracy order is an odd integer, we assign the
  // 'extra' point to the left of the point h.
  
  // Selects subset of 'x' array to use for finite difference method,
  int left_pt = (int)i - (int)floor((n+M)/2)+1;
  int right_pt = (int)i + (int)ceil((n+M)/2)+1;
  
  // Shifts right if sub-array is too far left:
  while (left_pt < 0)
  {
    left_pt++;
    right_pt++;
  }
  while ((Uint)right_pt >= N)
  {
    left_pt--;
    right_pt--;
  }
      
  // Allocates memory for delta coefficients
  // Stores (M+1)x(n+m)x(n+m) 3D array 
  // 'deltas' in linear format.
  Uint p[3] = { 0, M+1, n+M };
  double* deltas = alloc_double((M+1)*(n+M)*(n+M));
  
  // Fornberg algorithm calculates the deltas.
  deltas[i_j_k_to_ijk(p,0,0,0)] = 1;
  double c1 = 1;
  double c2 = 1;
  double c3 = 1;
  
  for (Uint l=1; l<n+M; l++)
  {
    c2 = 1;
    for (Uint v=0; v<l; v++)
    {
        c3 = x[(Uint)left_pt+l] - x[(Uint)left_pt+v];
        c2 = c2*c3;
        if (l <= M)
        { deltas[i_j_k_to_ijk(p,l,l-1,v)] = 0; }
        for (Uint m=0; m<=FDM_min(l,M); m++)
        { deltas[i_j_k_to_ijk(p,m,l,v)] = (((x[(Uint)left_pt+l] - h)
          *deltas[i_j_k_to_ijk(p,m,l-1,v)]
          - (m ? m*deltas[i_j_k_to_ijk(p,m-1,l-1,v)] : 0))/c3); }
    }
    for (Uint m=0; m<=FDM_min(l,M); m++)
    { deltas[i_j_k_to_ijk(p,m,l,l)] = ((c1/c2)*
      ((m ? m*deltas[i_j_k_to_ijk(p,m-1,l-1,l-1)] : 0)
      - (x[(Uint)left_pt+l-1]-h)*deltas[i_j_k_to_ijk(p,m,l-1,l-1)])); }
    c1 = c2;
  }
  
  // Approximates f^(m)|x=h using deltas.
  ret = 0;
  for (Uint v=0; v<=n; v++)
  { ret += deltas[i_j_k_to_ijk(p,M,n,v)] * f[left_pt+(int)v]; }
  
  free(deltas);
  return ret;
}
*/

// Manually-coded finite difference methods,
// only available for certain combinations of
// derivative and accuracy orders.

// Checks if method is available for given M, n.
static Uint interpolation_check_manual_FDM(Interpolation_T *const interp_s)
{
  if (!interp_s->finite_diff_order)
  { Error0("Interpolation finite difference error: Order of accuracy not assigned."); }
  if (!interp_s->FDM_derivative)
  { Error0("Interpolation finite difference error: Degree of derivative not assigned."); }
  const Uint n = interp_s->finite_diff_order;
  const Uint M = interp_s->FDM_derivative;
  const Uint N = *interp_s->N;
  
  // Checks if we have enough data points for given order.
  if (N <= n+M)
  {
    printf("Points: %i\n", N);
    printf("Order of accuracy: %i\n", n);
    Error0("Finite difference error: Not enough points for desired accuracy.\n");
    return 0;
  }
  if (M >= n)
  {
    printf("Degree of derivative: %i\n", M);
    printf("Order of accuracy: %i\n", n);
    Error0("Finite difference error: Degree of derivative must not exceed degree of accuracy.\n");
    return 0;
  }
  
  // Return success iff method available.
  if (M == 1)
  {
    if (n == 3)
    {
      return 1;
    }
  }
  
  printf("Manual FDM methods: %i-th derivative to order of accuracy %i not available.\n", M, n);
  Error0("Manual FDM invalid parameters.");
  return 0;
}

static double interpolation_manual_FDM(Interpolation_T *const interp_s)
{
  // Approximates M-th derivative of f(x)|x=h by finite difference method,
  // to order of accuracy n.
  const double *const x = *interp_s->x;
  const double *const f = *interp_s->f;
  const double h = *interp_s->h;
  const Uint N_total = *interp_s->N;
  Uint i = 0;
  const Uint n = interp_s->finite_diff_order;
  const Uint M = interp_s->FDM_derivative;
  double ret = DBL_MAX;
  Flag_T flg = NONE;
  
  if (!interpolation_check_manual_FDM(interp_s))
  {
    Error0("Interpolation: manual finite difference method failed.");
    return 0;
  }
  
  // First derivatives
  if (M == 1)
  {
    if (n == 3)
    {
      int N = 4;
      double d[16];
      double p[16];
      double q[4];
      ret = 0;
      
      // Finds the data segment
      for (i = 0; i < N_total-1; ++i)
      {
        if (GRTEQL(h,x[i]) && LSSEQL(h,x[i+1]))
        {
          flg = FOUND;
          break;
        }
      }
      if (flg != FOUND)
      {
        if (!interp_s->N_cubic_spline_1d->No_Warn)
          Warning("The given point for the interpolation is out of the domain.\n");
        
        return ret;
      }
      
      // Excise sub-arrays
      // Selects subset of 'x' array to use for finite difference method,
      // i.e. x[left_pt : right_pt].
      int left_pt = (int)i - (int)floor((n+M)/2)+1;
      int right_pt = (int)i + (int)ceil((n+M)/2)+1;
      
      // Shifts right if sub-array is too far left:
      while (left_pt < 0)
      {
        left_pt++;
        right_pt++;
      }
      while (right_pt >= (int)N_total)
      {
        left_pt--;
        right_pt--;
      }
      
      // 2D arrays are stored in linear format.
      // Array indices (j,k)
      int j;
      int k;
      int b;
      
      // Calculate d[j,k] array
      for (j = 0; j < N; j++)
      {
        for (k = 0; k < N; k++)
        { d[i_j_to_ij(N,j,k)] = x[left_pt + j] - x[left_pt + k]; }
      }
      
      // Calculate p[j,k] array
      for (j = 0; j < N; j++)
      {
        for (k = 0; k < N; k++)
        {
          p[i_j_to_ij(N,j,k)] = 1;
          for (b = 0; b < N; b++)
          {
            if (j != k && b != k && b != j)
            { p[i_j_to_ij(N,j,k)] *= (h - x[left_pt+b]) / d[i_j_to_ij(N,j,b)]; }
          }
        }
      }
      
      // Calculate l'_j(h) (= h[j]) for each j
      for (j = 0; j < N; j++)
      {
        q[j] = 0;
        for (k = 0; k < N; k++)
        {
          if (k != j)
          { q[j] += p[i_j_to_ij(N,j,k)] / d[i_j_to_ij(N,j,k)]; }
        }
      }
      
      // Calculate f'(h)
      ret = 0;
      for (j = 0; j < N; j++)
      ret += f[left_pt + j] * q[j];
      
      return ret;
    }
  }
  
  Error0("Manual FDM error: no method executed.");
  return ret;
}
     
///////////////////////////////////////Natural cubic spline///////////////////////////////
/* find a, b, c and d coeffs */
static void find_coeffs_natural_cubic_spline_1d(Interpolation_T *const interp_s)
{
  const Uint n = interp_s->N_cubic_spline_1d->N-1;
  const double *const x = interp_s->N_cubic_spline_1d->x;
  double *const a = interp_s->N_cubic_spline_1d->f;
  double *const b = alloc_double(n);
  double *const c = alloc_double(n+1);
  double *const d = alloc_double(n);
  double *h  = alloc_double(n),
         *l  = alloc_double(n+1),
         *z  = alloc_double(n+1),
         *al = alloc_double(n),
         *mu = alloc_double(n);
  Uint i;
  
  for (i = 0; i < n; ++i)
    h[i] = x[i+1]-x[i];
  for (i = 1; i < n; ++i)
    al[i] = 3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1];
  
  l[0]  = 1;
  mu[0] = 0;
  z[0]  = 0;
  for (i = 1; i < n; ++i)
  {
    l[i] = 2*(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
    mu[i] = h[i]/l[i];
    z[i] = (al[i]-h[i-1]*z[i-1])/l[i];
  }
  l[n] = 1;
  z[n] = 0;
  c[n] = 0;
  
  for (i = n-1; i >= 1; --i)
  {
    c[i] = z[i]-mu[i]*c[i+1];
    b[i] = (a[i+1]-a[i])/h[i]-h[i]*(c[i+1]+2.*c[i])/3.;
    d[i] = (c[i+1]-c[i])/3./h[i];
  }
  i = 0;
  c[i] = z[i]-mu[i]*c[i+1];
  b[i] = (a[i+1]-a[i])/h[i]-h[i]*(c[i+1]+2.*c[i])/3.;
  d[i] = (c[i+1]-c[i])/3./h[i];
  
  interp_s->N_cubic_spline_1d->a = a;
  interp_s->N_cubic_spline_1d->b = b;
  interp_s->N_cubic_spline_1d->c = c;
  interp_s->N_cubic_spline_1d->d = d;
  
  free(h);
  free(l);
  free(z);
  free(mu);
  free(al);
}

/* tutorial:
// =========
// given the values of the function f(x) for the points x_0, x_1, ... , and x_n,
// it calculates the value of f(h) for an h in [x_0, x_n]
// using Natural Cubic Spline method in 1-d
//
// NOTE: this algorithm doesn't work for an arbitrary order of x_i's.
// so ensure x_i's are in a decreasing or an increasing order.
//
// NOTE: Cubic spline doesn't use the values of df(x)/dx at x_i's. 
// hence, so in general, the interpolant derivatives don't match with the 
// true value of the derivatives(df/dx).
// if we need to include the derivative information, for instance,
// in EoS where d(epsl)/d(rho) should be consistent with pressure at T = 0, 
// we may want to use other methods such as Hermite which uses derivative.
// 
// ** filling the interpolation struct **
// Interpolation_T *interp_s = init_interpolation();
// interp_s->method          = "Natural_Cubic_Spline_1D";
// interp_s->N_cubic_spline_1d->f = f;
// interp_s->N_cubic_spline_1d->x = x;
// interp_s->N_cubic_spline_1d->h = h;// the point that we wanna interpolate f, i.e., f(h)
// interp_s->N_cubic_spline_1d->N = 20;// the dimension of the array f and x
// ** planning the appropriate function for interpolation **
// plan_interpolation(interp_s);
//
// ** evaluating interpolation **
// value = execute_interpolation(interp_s);
// ** freeing **
// free_interpolation(interp_s);
*/
static double interpolation_natural_cubic_spline_1d(Interpolation_T *const interp_s)
{
  if (!interp_s->N_cubic_spline_1d->Order)
    order_arrays_natural_cubic_spline_1d(interp_s);
    
  const double *const x = interp_s->N_cubic_spline_1d->x;
  const double *const a = interp_s->N_cubic_spline_1d->a;
  const double *const b = interp_s->N_cubic_spline_1d->b;
  const double *const c = interp_s->N_cubic_spline_1d->c;
  const double *const d = interp_s->N_cubic_spline_1d->d;
  const double h = interp_s->N_cubic_spline_1d->h;
  const double N = interp_s->N_cubic_spline_1d->N;
  double ret = DBL_MAX;/* it's important to be max double */
  Uint i = 0;
  Flag_T flg = NONE;
  
  /* find the segment */
  for (i = 0; i < N-1; ++i)
  {
    if (GRTEQL(h,x[i]) && LSSEQL(h,x[i+1]))
    {
      flg = FOUND;
      break;
    }
  }

  if (flg != FOUND)
  {
    if (!interp_s->N_cubic_spline_1d->No_Warn)
      Warning("The given point for the interpolation is out of the domain.\n");
    
    return ret;
  }
  
  ret = a[i]+b[i]*(h-x[i])+c[i]*Pow2(h-x[i])+d[i]*Pow3(h-x[i]);
  
  return ret; 
}

// Interpolates f'(x) using derivative of natural cubic spline.
static double interpolation_NCS_derivative(Interpolation_T *const interp_s)
{
  if (!interp_s->N_cubic_spline_1d->Order)
  order_arrays_natural_cubic_spline_1d(interp_s);
  
  const double *const x = interp_s->N_cubic_spline_1d->x;
  const double *const b = interp_s->N_cubic_spline_1d->b;
  const double *const c = interp_s->N_cubic_spline_1d->c;
  const double *const d = interp_s->N_cubic_spline_1d->d;
  const double h = interp_s->N_cubic_spline_1d->h;
  const double N = interp_s->N_cubic_spline_1d->N;
  double ret = DBL_MAX;/* it's important to be max double */
  Uint i = 0;
  Flag_T flg = NONE;
  
  /* find the segment */
  for (i = 0; i < N-1; ++i)
  {
    if (GRTEQL(h,x[i]) && LSSEQL(h,x[i+1]))
    {
      flg = FOUND;
      break;
    }
  }

  if (flg != FOUND)
  {
    if (!interp_s->N_cubic_spline_1d->No_Warn)
      Warning("The given point for the interpolation is out of the domain.\n");
    
    return ret;
  }
  
  ret = b[i] + 2*c[i]*(h-x[i]) + 3*d[i]*Pow2(h-x[i]);
  return ret; 
}

////////////////////////////////////////////Hermite Spline///////////////////////////////////////////////////////////
// Interpolates function f(x) over sub-intervals [a,b] of x domain.
// Given arrays for f(x) and x, each of length N,
// an n-th order Hermite spline will split the x domain into (N-1)overlapping
// sub-intervals of (n+1) points each, and construct an n-th degree interpolating
// polynomial Hj(x) on each sub-interval j (i.e. interval such that
// x[j] <= x < x[j+1]). Then f(x[j] <= h < x[j+1]) ~ Hj(h),
// and f(h) = Hj(h) exactly whenever h = x[j].
// The Hermite spline also interpolates the first derivative, i.e.
// f'(h) = Hj'(h) exactly whenever h = x[j].
//
// While generating the spline, 
// f'(x) can be estimated using finite difference methods.
//
// The spline order n must be greater than 2 (i.e. the lowest available
// spline is cubic). A very high n (n > 8 or so) is likely to result in
// oscillatory behavior and poor accuracy.
// A higher-order spline takes more time and memory to generate and store,
// and slightly longer to access once generated.

static void find_coeffs_Hermite_spline(Interpolation_T *const interp_s)
{
  // Finds values of the Hermite spline coefficients for each sub-interval,
  // For a spline of order G,
  // Hj(h) = a[0,j] + a[1,j]*(h-x0) + a[2,j]*(h-x0)^2 + a[3,j]*(h-x0)*(h-x1) + ...
  const double *const x = interp_s->Hermite_spline_1d->x;
  const double *const f = interp_s->Hermite_spline_1d->f;
  double *const fp = alloc_double(*interp_s->N);
  const int N = (int)interp_s->Hermite_spline_1d->N;
  const int n = (int)interp_s->Spline_Order - 1;
  double *const a = alloc_double(((Uint)N-1)*(2*(Uint)n+2));
  // a[] is a linearized 2D array where a[i,j] is
  // the j-th spline coefficient over interval i.
  
  // Error if not enough data points for desired spline order.
  if (n+1 > N)
  { 
    printf("Hermite spline order: %i\n", n+1);
    printf("Hermite spline data points: %i\n", *interp_s->N);
    Error0("Hermite spline error: Not enough data points for desired spline order.");
  }
  
  // Estimates f'(x) for each point in the domain x
  // using finite difference method.
  // Stores them for later use in generating the Hermite spline.
  interp_s->FDM_derivative = 1;
  for (int i = 0; i < N; i++)
  {
    interp_s->Hermite_spline_1d->h = x[i];
    fp[i] = interpolation_finite_difference(interp_s);
  }
  interp_s->Hermite_spline_1d->fp = fp;
  
  // Calculates spline coefficients on each interval x[j] <= x < x[j+1].
  // Temporary arrays
  double *Q = alloc_double((2*(Uint)n+2)*(2*(Uint)n+2));
  double *z = alloc_double(2*(Uint)n+2);
  int l; // Left-most sub-array point
  for (int j = 0; j < N; j++)
  {
    for (Uint c = 0; c < (2*(Uint)n+2)*(2*(Uint)n+2); c++)
    { Q[c] = 0; }
    // Excise sub-array from x[]. Shift left if out of bounds.
    l = j;
    while (l >= N - (n+1))
    { l--; }
    
    //Step 1
    for (int i = 0; i <= n; i++)
    {
      //Step 2
      z[2*i]                      = x[l+i];
      z[2*i+1]                    = x[l+i];
      Q[i_j_to_ij(2*n+1,2*i,0)]   = f[l+i];
      Q[i_j_to_ij(2*n+1,2*i+1,0)] = f[l+i];
      Q[i_j_to_ij(2*n+1,2*i+1,1)] = fp[l+i];

      //Step 3
      if (i)
      { Q[i_j_to_ij(2*n+1,2*i,1)] = (Q[i_j_to_ij(2*n+1,2*i,0)]
      - Q[i_j_to_ij(2*n+1,2*i-1,0)]) / (z[2*i] - z[2*i-1]); }
    }
    
    //Step 4
    for (int b = 2; b <= 2*n+1; b++)
    {
      for (int k = 2; k <= b; k++)
      { Q[i_j_to_ij(2*n+1,b,k)] = (Q[i_j_to_ij(2*n+1,b,k-1)]
      - Q[i_j_to_ij(2*n+1,b-1,k-1)]) / (z[b] - z[b-k]); }
    }

    //Step 5
    for (int b = 0; b <= 2*n+1; b++)
    { a[i_j_to_ij(2*n+1,j,b)] = Q[i_j_to_ij(2*n+1,b,b)]; }
  }
  
  interp_s->Hermite_spline_1d->Alloc_Mem = 1;
  interp_s->Hermite_spline_1d->a = a;
  free(Q);
  free(z);
}

static double interpolation_Hermite_spline(Interpolation_T *const interp_s)
{
  // Executes Hermite spline interpolation:
  // f(h) ~ H(h) = Q00 + Q11(h-x0) + Q22(h-x0)^2
  //        + Q33(h-x0)^2(h-x1) + ... + Q2n+1,2n+1(h-x0)^2...(h-xn).
  if (!interp_s->Hermite_spline_1d->Order)
  order_arrays_spline_1d(interp_s);
    
  const double *const x = interp_s->Hermite_spline_1d->x;
  const double *const a = interp_s->Hermite_spline_1d->a;
  const double h = interp_s->Hermite_spline_1d->h;
  const int N = (int)interp_s->Hermite_spline_1d->N;
  const int n = (int)interp_s->Spline_Order - 1;
  double ret = DBL_MAX;// it's important to be max double
  int i = 0;
  Flag_T flg = NONE;
  
  // Finds the spline interval
  for (i = 0; i < N-1; ++i)
  {
    if (GRTEQL(h,x[i]) && LSSEQL(h,x[i+1]))
    {
      flg = FOUND;
      break;
    }
  }

  if (flg != FOUND)
  {
    if (!interp_s->Hermite_spline_1d->No_Warn)
    {
      printf("Domain: [%E, %E]\n", x[0], x[N-1]);
      printf("Point for interpolation: %E\n", h);
      Warning("The given point for the interpolation is out of the domain.\n");
    }
    
    return ret;
  }
 
  int l = i;
  while (l >= N - (n+1))
  { l--; }
  
  ret = 0;
  double ret_term = 1;
  // ADD COMMENT
  for (int j = 0; j <= 2*n+1; j++)
  {
    ret += a[i_j_to_ij(2*n+1,l,j)]*ret_term;
    ret_term *= (h - x[l+(int)floor(j/2)]);
  }
  
  return ret;
}

/////////////////////////////////////////Neville interpolation
/* tutorial:
// =========
// given function f(xi)'s and points xi's it calculates value of f(h)
// using Neville iterate method in 1-d
// 
// ** filling the interpolation struct **
// Interpolation_T *interp_s = init_interpolation();
// interp_s->method         = "Neville_1D";
// interp_s->Neville_1d->f   = f;
// interp_s->Neville_1d->x   = x;
// interp_s->Neville_1d->h   = h;// the point that we wanna interpolate f i.e. f(h)
// interp_s->Neville_1d->N   = 20;// dimension of array f and x
// interp_s->Neville_1d->max = 10;// 10 out of N used for interpolation
// ** planning the appropriate function for interpolation **
// plan_interpolation(interp_s);
//
// ** evaluating interpolation **
// value = execute_interpolation(interp_s);
// ** freeing **
// free_interpolation(interp_s);
*/
static double interpolation_Neville_1d(Interpolation_T *const interp_s)
{
  if (!interp_s->Neville_1d->f)
    Error0("No function is specified.\n");
  
  const double *const f = interp_s->Neville_1d->f;
  double *x;
  const double h = interp_s->Neville_1d->h;
  double **q;
  double interp = 0;
  Matrix_T *Q   = 0;
  Uint N,i,j;
  Flag_T flg = NO;
  
  if (interp_s->Neville_1d->max    			  && 
      interp_s->Neville_1d->N > interp_s->Neville_1d->max &&
      /* for small N we might get x[i]-x[i-j] = 0 thus: */
      interp_s->Neville_1d->N > 10)
  {
    flg = YES;
    N = interp_s->Neville_1d->max;
    Uint s = interp_s->Neville_1d->N/N;
    assert(s);
    Q = alloc_matrix(REG_SF,N,N);
    q = Q->reg->A;
    
    x = alloc_double(N);
    x[0]      = interp_s->Neville_1d->x[0];
    x[N-1]    = interp_s->Neville_1d->x[N-1];
    q[0][0]   = f[0];
    q[N-1][0] = f[N-1];
    for (i = 1; i < N-1; i++)
    {
      q[i][0] = f[i*s];
      x[i] = interp_s->Neville_1d->x[i*s];
    }
  }
  else
  {
    N = interp_s->Neville_1d->N;
    Q = alloc_matrix(REG_SF,N,N);
    q = Q->reg->A;
    x = interp_s->Neville_1d->x;
    
    for (i = 0; i < N; i++)
      q[i][0] = f[i];
  }

  for (i = 1; i < N; ++i)
    for (j = 1; j <= i; ++j)
      q[i][j] = ((h-x[i-j])*q[i][j-1]-(h-x[i])*q[i-1][j-1])/(x[i]-x[i-j]);

  interp = q[N-1][N-1];

  free_matrix(Q);
  if (flg == YES)
    free(x);
  
  return interp;
}

/* calculating interpolation using Chebyshev spectral method and
// Chebyshev Extrema collocation points.
// ->return value: interpolation value
*/
static fInterpolation_T *interpolation_Chebyshev_Tn(Interpolation_T *const interp_s)
{
  fInterpolation_T *func = 0;
  Uint count = 0;
  
  if (interp_s->XYZ_dir_flag) count++;
  if (interp_s->XY_dir_flag)  count++;
  if (interp_s->XZ_dir_flag)  count++;
  if (interp_s->YZ_dir_flag)  count++;
  if (interp_s->X_dir_flag)   count++;
  if (interp_s->Y_dir_flag)   count++;
  if (interp_s->Z_dir_flag)   count++;
  
  if (count > 1)
    Error0("For spectral interpolation, more than one direction was flagged.\n");
  
  if (interp_s->XYZ_dir_flag)
    func = interpolation_Chebyshev_Tn_XYZ;
    
  else if (interp_s->XY_dir_flag)
    func = interpolation_Chebyshev_Tn_XY;
    
  else if (interp_s->XZ_dir_flag)
    func = interpolation_Chebyshev_Tn_XZ;
    
  else if (interp_s->YZ_dir_flag)
    func = interpolation_Chebyshev_Tn_YZ;
    
  else if (interp_s->X_dir_flag)
    func = interpolation_Chebyshev_Tn_X;
    
  else if (interp_s->Y_dir_flag)
    func = interpolation_Chebyshev_Tn_Y;
    
  else if (interp_s->Z_dir_flag)
    func = interpolation_Chebyshev_Tn_Z;
    
  else
    Error0("The flags in interpolation structure were not set correctly.\n");
  
  return func;
}

/* interpolation in X direction keeping the other two directions constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_X(Interpolation_T *const interp_s)
{
  //printf("\tinterpolation_Chebyshev_Tn_X:\n");//////////
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = General2ChebyshevExtrema(interp_s->X,0,field->patch);/* X coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const Uint J = interp_s->J;/* index of const coords */
  const Uint K = interp_s->K;/* index of const coords */
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint i;
  
  make_coeffs_1d(field,0);
  C = field->v2;
  
  for (i = 0; i < n[0]; ++i)
  {
    interp_v += C[i_j_k_to_ijk(n,i,J,K)]*Tx(i,X);
    //printf("\t\tC[%i][%i][%i] == %E\n", i, J, K, C[i]);//////////
    //printf("\t\tTx(%i,X) == %E\n", i, Tx(i,X));//////////
    //printf("\t\tinterp == %E\n", interp_v);//////////
  }
  
  return interp_v;
}

/* interpolation in Y direction keeping the other two directions constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_Y(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double Y = General2ChebyshevExtrema(interp_s->Y,1,field->patch);/* Y coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const Uint I = interp_s->I;/* index of const coords */
  const Uint K = interp_s->K;/* index of const coords */
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint j;
  
  make_coeffs_1d(field,1);
  C = field->v2;
  
  for (j = 0; j < n[1]; ++j)
    interp_v += C[i_j_k_to_ijk(n,I,j,K)]*Ty(j,Y);
  
  return interp_v;
}

/* interpolation in Z direction keeping the other two directions constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_Z(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double Z = General2ChebyshevExtrema(interp_s->Z,2,field->patch);/* Z coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const Uint I = interp_s->I;/* index of const coords */
  const Uint J = interp_s->J;/* index of const coords */
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint k;
  
  make_coeffs_1d(field,2);
  C = field->v2;
  
  for (k = 0; k < n[2]; ++k)
    interp_v += C[i_j_k_to_ijk(n,I,J,k)]*Tz(k,Z);
  
  return interp_v;
}

/* interpolation in X&Y directions keeping the other direction constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_XY(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = General2ChebyshevExtrema(interp_s->X,0,field->patch);/* X coord of the interesting point */
  const double Y = General2ChebyshevExtrema(interp_s->Y,1,field->patch);/* Y coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const Uint K = interp_s->K;/* index of const coords */
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint i,j;
  
  make_coeffs_2d(field,0,1);
  C = field->v2;
  
  for(i = 0; i < n[0]; ++i)
    for (j = 0; j < n[1]; ++j)
      interp_v += C[i_j_k_to_ijk(n,i,j,K)]*Tx(i,X)*Ty(j,Y);
      
  return interp_v;
}

/* interpolation in X&Z directions keeping the other direction constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_XZ(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = General2ChebyshevExtrema(interp_s->X,0,field->patch);/* X coord of the interesting point */
  const double Z = General2ChebyshevExtrema(interp_s->Z,2,field->patch);/* Z coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const Uint J = interp_s->J;/* index of const coords. */
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint i,k;
  
  make_coeffs_2d(field,0,2);
  C = field->v2;
  
  for(i = 0; i < n[0]; ++i)
    for (k = 0; k < n[2]; ++k)
      interp_v += C[i_j_k_to_ijk(n,i,J,k)]*Tx(i,X)*Tz(k,Z);
      
  return interp_v;
}

/* interpolation in Y&Z direction keeping the other direction constant.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_YZ(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double Y = General2ChebyshevExtrema(interp_s->Y,1,field->patch);/* Y coord of the interesting point */
  const double Z = General2ChebyshevExtrema(interp_s->Z,2,field->patch);/* Z coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const Uint I = interp_s->I;/* index of const coords. */
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint j,k;
  
  make_coeffs_2d(field,1,2);
  C = field->v2;
  
  for(j = 0; j < n[1]; ++j)
    for (k = 0; k < n[2]; ++k)
      interp_v += C[i_j_k_to_ijk(n,I,j,k)]*Ty(j,Y)*Tz(k,Z);
      
  return interp_v;
}
/* interpolation in X&Y&Z directions.
// ->return value: interpolation value.
*/
static double interpolation_Chebyshev_Tn_XYZ(Interpolation_T *const interp_s)
{
  Field_T *const field = interp_s->field;/* interesting field */
  const double X = General2ChebyshevExtrema(interp_s->X,0,field->patch);/* X coord of the interesting point */
  const double Y = General2ChebyshevExtrema(interp_s->Y,1,field->patch);/* Y coord of the interesting point */
  const double Z = General2ChebyshevExtrema(interp_s->Z,2,field->patch);/* Z coord of the interesting point */
  const Uint *const n = interp_s->field->patch->n;
  const double *C = 0;/* coeffs of expansion of field in Tn */
  double interp_v = 0;
  Uint i,j,k;
  
  make_coeffs_3d(field);
  C = field->v2;
  
  for (i = 0; i < n[0]; ++i)
    for (j = 0; j < n[1]; ++j)
      for (k = 0; k < n[2]; ++k)
        interp_v += C[i_j_k_to_ijk(n,i,j,k)]*Tx(i,X)*Ty(j,Y)*Tz(k,Z);
  
  return interp_v;
}

/* basis of expansion of interpolation with 
// "first kind Chebyshev" with extrema collocation
// which combined with some coefficients and simplifications
// for ease of using at interpolation implementation.
// ->return value: value of this compound basis at specific point.
*/
static double T(const Uint n,const Uint i,const double x)
{
  double t = 0;
  
  if (i == 0)
    t = 1;
  else if (i == n-1)
    t = Cheb_Tn((int)i,x);
  else
    t = 2*Cheb_Tn((int)i,x);
  
  return t;
}

/* freeing interpolation structure. */
void free_interpolation(Interpolation_T *interp_s)
{
  if (!interp_s)
    return;
  if (strstr_i(interp_s->method,"Natural_Cubic_Spline_1D"))
  {
    if (interp_s->N_cubic_spline_1d->b)
      free(interp_s->N_cubic_spline_1d->b);
    if (interp_s->N_cubic_spline_1d->c)
      free(interp_s->N_cubic_spline_1d->c);
    if (interp_s->N_cubic_spline_1d->d)
      free(interp_s->N_cubic_spline_1d->d);
    if (interp_s->N_cubic_spline_1d->Alloc_Mem)
    {
      free(interp_s->N_cubic_spline_1d->x);
      free(interp_s->N_cubic_spline_1d->f);
    }
  }
  else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline")
           || strstr_i(interp_s->method, "Hermite"))
  {
    free(interp_s->Hermite_spline_1d->fp);
    free(interp_s->Hermite_spline_1d->a);
    if (interp_s->Hermite_spline_1d->Alloc_Mem)
    {
      free(interp_s->Hermite_spline_1d->x);
      free(interp_s->Hermite_spline_1d->f);
      interp_s->Hermite_spline_1d->f = 0;
      interp_s->Hermite_spline_1d->x = 0;
      *interp_s->f = 0;
      *interp_s->x = 0;
      //free(interp_s->Hermite_spline_1d->fp);
      //free(interp_s->Hermite_spline_1d->a);
    }
  }
  
  free(interp_s);
}

// Getter and setter methods for bit-field flags in 1D interpolation
// sub-structures. These are necessary since the bit-fields cannot
// be accessed by reference.

// Note: getter methods currently aren't being used anywhere.
// Might be useful in future applications.
/*
// Gets array order flag
static Uint get_order_flag(Interpolation_T *const interp_s)
{
  if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
  { return (int)interp_s->N_cubic_spline_1d->Order; }
  else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
  { return (int)interp_s->Hermite_spline_1d->Order; }
  
  return 0;
}

static Uint get_warn_flag(Interpolation_T *const interp_s)
{
  if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
  { return (int)interp_s->N_cubic_spline_1d->No_Warn; }
  else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
  { return (int)interp_s->Hermite_spline_1d->No_Warn; }
  
  return 0;
}

static Uint get_alloc_mem_flag(Interpolation_T *const interp_s)
{
  if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
  { return (int)interp_s->N_cubic_spline_1d->Alloc_Mem; }
  else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
  { return (int)interp_s->Hermite_spline_1d->Alloc_Mem; }
  
  return 0;
}
*/

// Sets array order flag
static void set_interp_order_flag(Interpolation_T *const interp_s, Uint flag)
{
  if (flag != 0 && flag != 1)
  { Error0("set_order_flag: invalid flag value."); }
  
  // The if statements are necessary because the Order flag
  // is a bit field and direct conversion of the flag value
  // could be machine-dependent.
  if (flag)
  {
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->Order = 1; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
    { interp_s->Hermite_spline_1d->Order = 1; }
  }
  else
  {
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->Order = 0; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
    { interp_s->Hermite_spline_1d->Order = 0; }
  }
}

// Sets interpolation warning flag
void set_interp_warn_flag(Interpolation_T *const interp_s, Uint flag)
{
  if (flag != 0 && flag != 1)
  { Error0("set_order_flag: invalid flag value."); }
  
  // The if statements are necessary because the No_Warn flag
  // is a bit field and direct conversion of the flag value
  // could be machine-dependent.
  if (flag)
  {
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->No_Warn = 1; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
    { interp_s->Hermite_spline_1d->No_Warn = 1; }
  }
  else
  {
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->No_Warn = 0; }
    else if (strstr_i(interp_s->method, "Hermite_Cubic_Spline"))
    { interp_s->Hermite_spline_1d->No_Warn = 0; }
  }
}

// Sets interpolation alloc_mem flag
void set_interp_alloc_mem_flag(Interpolation_T *const interp_s, Uint flag)
{
  if (flag != 0 && flag != 1)
  { Error0("set_order_flag: invalid flag value."); }
  
  // The if statements are necessary because the No_Warn flag
  // is a bit field and direct conversion of the flag value
  // could be machine-dependent.
  if (flag)
  {
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->Alloc_Mem = 1; }
    else if (strstr_i(interp_s->method, "Hermite"))
    { interp_s->Hermite_spline_1d->Alloc_Mem = 1; }
  }
  else
  {
    if (strstr_i(interp_s->method, "Natural_Cubic_Spline_1D"))
    { interp_s->N_cubic_spline_1d->Alloc_Mem = 0; }
    else if (strstr_i(interp_s->method, "Hermite"))
    { interp_s->Hermite_spline_1d->Alloc_Mem = 0; }
  }
}

