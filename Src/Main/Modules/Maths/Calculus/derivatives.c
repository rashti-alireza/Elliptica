/*
// Alireza Rashti
// July 2018
*/

#include "derivatives.h"
#define DELIMIT '|'
#define COMMA ','
#define MAX_STR_LEN 100

/* Covariant Derivative function.
// it takes covariant derivative of a given field.
//
// usage examples: 
// ==============
// double *Dphi     = Covariant_Derivative(phi,"x");// => d(phi)/dx
// double *DxA_U2   = Covariant_Derivative(A_U2,"x");// => D(A_U2)/Dx = d(A(2))/dx + Gamma(2,0,k)*A(-k)
// double *DzB_U1D0 = Covariant_Derivative(B_U1D2,"z");// => D(B_U1D2)/Dz = d(B(1,2))/dz + Gamma(1,2,k)*B(-k,2) - Gamma(k,2,2)*B(1,-k)
//
// NOTE1 : each field has been defined on just one patch.
// NOTE2 : this function allocate memory for the result.
// NOTE3 : "x","y","z" are Cartesian coords
*/
double *Covariant_Derivative(Field_T *const f,const char *task)
{
  /* The idea is to take the name of the field and then make the possible
  // strings from its different components, e.g, in case f->name = "A_U1D2" 
  // we need all component "A_U?D?", having made strings then
  // a function can calculate the partial derivativce part and a loop can calculate
  // Gamma*A_U?D? parts for each index. */
  UNUSED(f);
  UNUSED(task);
  abortEr(NO_JOB);
  return 0;
}

/* partial derivative function. it takes partial derivative for each field
// according to the task and return the result.
//
// usage example:
// ==============
// double *dphi = Partial_Derivative(phi,"x"); //=> taking d(phi)/dx
//
//
// NOTE1 : each field has been defined on just one patch.
// NOTE2 : this function also in case of spectral derivative makes the 
// spectral coefficients and save them in f->v2;
// NOTE3 : this function allocate memory for the result.
//
// list of available derivative methods:
// =====================================
// o. spectral = for spectral method
//
// task format example:
// ===================
//
// Cartesian 'x' derivative => *task = "x"
// Cartesian 'z' followed by Cartesian 'y' derivative => *task = "z,y"
// Cartesian xx derivative => *task = "x,x"
// and so forth. furthermore for derivative in curvilinear coords we have:
// Curvilinear 'a' derivative followed by 'b' followed by 'c' => *task ="a,b,c"
// if one wants to override the default derivative method defined in 
// the input file, they can append the task by " DELIMIT derivative type"; e.g.
// Cartesian x derivative with finite difference:
// => *task = "x DELIMIT Finite_Difference", DELIMIT is a macro defined above.
// so for example if DELIMIT is | then *task = "x | Finite_Difference".
// Note: curvilinear coordinate must be identified by a,b,c not X,Y,Z.
//
// ->return value: partial derivative of Field_T f accordingly, null for error.
*/
double *Partial_Derivative(Field_T *const f,const char *task)
{
  /* check up */
  if (!f)
    abortEr("The field is empty!\n");
  if (!f->v && !f->v2)  
    abortEr("The field is empty!\n");
  if (!task)
    abortEr("The task is empty!\n");

  double *r = 0;
  const char *der_par = PgetsEZ("Derivative_Method");
  unsigned Ndir;
  Method_T method_e = derivative_method(der_par,task);
  Dd_T  *dir_e = derivative_direction(task,&Ndir);
  
  if (method_e == SPECTRAL)
    r = take_spectral_derivative(f,dir_e,Ndir);
  else
    abortEr("There is no such derivative method defined for this function.\n");
    
  free(dir_e);
  return r;
}

/* taking spectral derivative of the given field 
// based on direction of derivative.
// ->return value: values of the resultant.
*/
static double *take_spectral_derivative(Field_T *const f,const Dd_T  *const dir_e,const unsigned Ndir)
{
  double *deriv = 0;
  Field_T *ff[2];
  unsigned bck,frd;
  Patch_T p_tmp1 = make_temp_patch(f->patch);
  Patch_T p_tmp2 = make_temp_patch(f->patch);
  unsigned i;
  
  assert(Ndir);
  
  /* 3-D fields */
  if (strstr(f->attr,"(3dim)"))
  {
    /* check if it is a second order derivative and it is
    // possible to apply second oreder derivative formulas.
    */
    if (IsSecondOrderFormula(f,dir_e,Ndir))
      deriv = spectral_derivative_2ndOrder(f,dir_e[0]);
    else
    {
      ff[0] = add_field("tmp1","(3dim)",&p_tmp1,NO);
      ff[1] = add_field("tmp2","(3dim)",&p_tmp2,NO);
      
      deriv = spectral_derivative_1stOrder(f,dir_e[0]);
      ff[0]->v = deriv;
      frd = 0;
      for (i = 1; i < Ndir; ++i)
      {
        frd = i%2;
        bck = i/2;
        deriv = spectral_derivative_1stOrder(ff[bck],dir_e[i]);
        
        /* make next ff ready */
        ff[frd]->v = deriv;
        free_v(ff[bck]);
        free_v2(ff[bck]);
      }
      
      ff[frd]->v = 0;
      /* free leftovers */
      remove_field(ff[0]);
      remove_field(ff[1]);
      
      free_temp_patch(&p_tmp1);
      free_temp_patch(&p_tmp2);
    }
  }/* end of if (strstr(f->attr,"(3dim)")) */
  else
    abortEr("No such Dimension is defined for this function.\n");
  
  return deriv;
}

/* find out if one can apply second order derivative formula for
// this derivative.
// ->return value: 1 if applicable , 0 otherwise.
*/
static unsigned IsSecondOrderFormula(Field_T *const f,const Dd_T *const dir_e,const unsigned Ndir)
{
  unsigned r = 0;
  SpecDerivative_Func_T *df[3];
  Dd_T dp[3];
  int c;
  unsigned d;
  
  if (Ndir != 2)/* if it doesn't have 2 derivatives */
    return 0;
  else if (dir_e[0] != dir_e[1])/* if directions are not matched */
    return 0;
  else
  {
    get_SpecDerivative_func_2ndOrder(f->patch,df);
    c = get_dp_2ndOrder(f->patch,df,dir_e[0],dp);
    if (c > 1)/* if it depends more than 1 variable */
      return 0;
    
    c = -1;
    for (d = 0; d < 3; ++d)
    {
      if (dp[d] != UNDEFINED_DIR)
        c = (int)d;
    }

    /* since it is only depends on one directions
    // so second order formula is likely applicable.
    // further study required.
    */
    r = 0;
    if (c != -1)
    {
      r = JacobianFormat_2ndOrder(f->patch,dir_e[0],dp[c]);
    }
  }
  
  return r;
}

/* study if the format of Jacobian is matched with second order derivative.
// it must be like that: d^2f(X)/dX^2 = (dN/dX)^2 d^2f/dN^2. if the transformation
// is not like that, Jacobian format is not matched.
// ->return value: 1 if matches, 0 otherwise.
*/
static unsigned JacobianFormat_2ndOrder(const Patch_T *const patch,const Dd_T dir,Dd_T dp)
{
  unsigned r = 0;
  
  if (patch->coordsys == Cartesian)
  {
    if (dir%3 == dp%3)
      r = 1;
  }
  else
     abortEr(NO_JOB);
  
  return r;
}

/* based on basis and collocation,
// get the pertinent function for second order spectral derivative. 
*/
static void get_SpecDerivative_func_2ndOrder(const Patch_T *const patch,SpecDerivative_Func_T **func)
{
  unsigned i;
  
  for (i = 0; i < 3; ++i)
  {
    if (patch->basis[i] == Chebyshev_Tn_BASIS 
        && patch->collocation[i] == Chebyshev_Extrema)
    {
      func[i] = derivative_ChebyshevExtrema_Tn_2ndOrder;
    }
    else if (patch->basis[i] == Chebyshev_Tn_BASIS 
        && patch->collocation[i] == Chebyshev_Nodes)
    {
      func[i] = derivative_ChebyshevNodes_Tn_2ndOrder;
    }
    else
      abortEr("There is no such basis or collocation defined for this function.\n");
  }
}

/* based on second order spectral derivative function and dependencies 
// it finds the Dd_T dp used for direction of partial derivative.
// also it returns the number of dependecy of dir on a,b,c, which obviously
// if it greaters than 1, it means that 2nd order formula is not applicable.
// ->return value: number of dependency
*/
static int get_dp_2ndOrder(const Patch_T *const patch,SpecDerivative_Func_T **func,const Dd_T dir,Dd_T *dp)
{
  unsigned depend[3];
  unsigned i;
  int cnt;
  
  get_dependency(patch,dir,depend);
  
  cnt = 0;
  for (i = 0; i < 3; ++i)
  {
    dp[i] = UNDEFINED_DIR;
    if (depend[i])/* if this direction depends on _a_,_b_,_c_ */
    {
      if (func[i] == derivative_ChebyshevExtrema_Tn_2ndOrder)
      {
        dp[i] = (Dd_T)i;/* means _N0_or _N1_or _N2_ */
        cnt++;
      }
      else if (func[i] == derivative_ChebyshevNodes_Tn_2ndOrder)
      {
        dp[i] = (Dd_T)i;/* means _N0_or _N1_or _N2_ */
        cnt++;
      }
      else
        abortEr("There is no such derivative function defined for this function.\n");

    }
  }
  
  return cnt;
}

/* taking second order 3-D spectral derivative in 
// the specified direction dir_e on the patch determined by field.
// it is worth explaining that some bases are expanded in 
// different coordinate than x,y and z or a,b and c; for example Chebyshev Tn
// which expanded in normal coords N0,N1 and N2. so in taking
// derivative one needs to consider this. the variable being used
// for this is "Dd_T dp[3]".
// note: it allocates memory for resultant in each invoking.
// note: second order means that: d^2f(x)/dx^2 = (dN/dx)^2 *d^2f(N)/dN^2
// and one can take d^2f(N)/dN^2 by some special formula. for example:
// in Chebyshev Tn we have d^2Tn(N)/d^2N = some known formula.
// so one can take advantage of this to calculate the derivative.
// ->return value: second order spectral derivative.
*/
static double *spectral_derivative_2ndOrder(Field_T *const f,const Dd_T dir_e)
{
  Patch_T *const patch = f->patch;
  double *der = alloc_double(patch->nn);
  SpecDerivative_Func_T *df[3];/* spectral derivative function in each direction */
  Dd_T dp[3];/* see above explanation */
  double *df_dp[3];
  unsigned nn = patch->nn;
  unsigned i ,c;
  Dd_T d;
  Flag_T flg[3];
  
  get_SpecDerivative_func_2ndOrder(patch,df);
  get_dp_2ndOrder(patch,df,dir_e,dp);
  
  c = 0;
  for (d = 0; d < 3; ++d)
  {
    flg[d] = NO;
    df_dp[d] = 0;
    if (dp[d] != UNDEFINED_DIR)
    {
      df_dp[d] = df[d](f,dp[d]);/* taking derivative with respect to dp[d] */
      flg[d] = YES;
      c++;
    }
  }
  
  /* since this is second order, we have to have only one direction */
  assert(c == 1);
  
  if (flg[0] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
    {
      double j = dq2_dq1(patch,dp[0],dir_e,i);
      der[i] = df_dp[0][i]*SQR(j);
    }
  }
  else if (flg[1] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
    {
      double j = dq2_dq1(patch,dp[1],dir_e,i);
      der[i] = df_dp[1][i]*SQR(j);
    }
  }
  else if (flg[2] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
    {
      double j = dq2_dq1(patch,dp[2],dir_e,i);
      der[i] = df_dp[2][i]*SQR(j); 
    }
  }


  for (d = 0; d < 3; ++d)
    if (df_dp[d])
      free(df_dp[d]);
  
  return der;
}

/* finding all of types of derivatives and put them into 
// an array of Dd_T with size *n in order that they have been written.
// note: this function allocate memory.
// ->return value: array of Dd_T and number of this arrays
*/
static Dd_T *derivative_direction(const char *const task,unsigned *const n)
{
  Dd_T *e = 0;
  char *savestr,*str = dup_s(task);
  char *tok = tok_s(str,DELIMIT,&savestr);
  
  if (!tok)
    abortEr_s("There is No direction in %s.\n",task);
  
  *n = 0;
  tok = tok_s(tok,COMMA,&savestr);  
  while (tok)
  {
    e = realloc(e,(*n+1)*sizeof(*e));
    pointerEr(e);
    e[*n] = str2enum_direction(tok);
    tok = tok_s(0,COMMA,&savestr);
    ++(*n);
  }
  
  free(str);
  
  return e;
}

/* getting both task and par, it finds the derivative method.
// if task has a derivative method, it is preferred over par.
// ->return value: derivative method in enum format.
*/
static Method_T derivative_method(const char *const par,const char *const task)
{
  Method_T type = UNDEFINED_METHOD;
  char *s = dup_s(task);
  char *rs = 0;
  
  tok_s(s,DELIMIT,&rs);
  type = str2enum_method(rs);
  free(s);
  
  /* check parameter if no info is in task */
  if (type == UNDEFINED_METHOD)
  {
    type = str2enum_method(par);
  } 
  
  /* if still no type has been found => no info in parameter and task */
  if (type == UNDEFINED_METHOD)
    abortEr("No Derivative Method is defined "
      "in parameter file or in function task.\n");
  
  return type;
}

/* getting a derivative method in string format and 
// returing the Method_T.
// ->return value: Method_T and if not found UNDEFINED_Method.
*/
static Method_T str2enum_method(const char *const str)
{
  Method_T type = UNDEFINED_METHOD;
  
  if (strstr_i(str,"Spectral"))
    type = SPECTRAL;
  //else if (strcmp_i(str,"Finite_Difference"))
    //type = FINITE_DIFF;
  
  return type;
}

/* getting a derivative direction in string format and 
// returning the Dd_T.
// ->return value: Dd_T and if not found error.
*/
static Dd_T str2enum_direction(const char *const str)
{
  if (strcmp_i(str,"x"))
    return _x_;
  else if (strcmp_i(str,"y"))
    return _y_;
  else if (strcmp_i(str,"z"))
    return _z_;
  else if (strcmp_i(str,"a"))
    return _a_;
  else if (strcmp_i(str,"b"))
    return _b_;
  else if (strcmp_i(str,"c"))
    return _c_;
  else
    abortEr_s("There is no such %s derivative defined!\n",str);
  
  return UNDEFINED_DIR;
}

/* taking 3-D spectral derivative in the specified direction dir_e 
// on a patch determined by field.
// it is worth explaining that some bases are expanded in 
// different coordinate than x,y and z or a,b and c; for example Chebyshev Tn
// which expanded in normal coords N0,N1 and N2. so in taking
// derivative one needs to consider this. the variable being used
// for this is "Dd_T dp[3]".
// note: it allocates memory for resultant in each invoking.
// ->return value: sepctral derivative.
*/
static double *spectral_derivative_1stOrder(Field_T *const f,const Dd_T dir_e)
{
  Patch_T *const patch = f->patch;
  double *der = alloc_double(patch->nn);
  SpecDerivative_Func_T *df[3];/* spectral derivative function in each direction */
  Dd_T dp[3];/* see above explanation */
  double *df_dp[3];
  unsigned nn = patch->nn;
  unsigned i;
  Dd_T d;
  Flag_T flg[3];
  
  get_SpecDerivative_func_1stOrder(patch,df);
  get_dp_1stOrder(patch,df,dir_e,dp);
  
  for (d = 0; d < 3; ++d)
  {
    flg[d] = NO;
    df_dp[d] = 0;
    if (dp[d] != UNDEFINED_DIR)
    {
      df_dp[d] = df[d](f,dp[d]);/* taking derivative with respect to dp[d] */
      flg[d] = YES;
    }
  }
  
  if (flg[0] == YES && flg[1] == YES && flg[2] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
      der[i] = df_dp[0][i]*dq2_dq1(patch,dp[0],dir_e,i) + 
               df_dp[1][i]*dq2_dq1(patch,dp[1],dir_e,i) +
               df_dp[2][i]*dq2_dq1(patch,dp[2],dir_e,i);
  }
  else if (flg[0] == YES && flg[1] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
      der[i] = df_dp[0][i]*dq2_dq1(patch,dp[0],dir_e,i) + 
               df_dp[1][i]*dq2_dq1(patch,dp[1],dir_e,i);
  }
  else if (flg[1] == YES && flg[2] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
      der[i] = df_dp[1][i]*dq2_dq1(patch,dp[1],dir_e,i) + 
               df_dp[2][i]*dq2_dq1(patch,dp[2],dir_e,i);
  }
  else if (flg[0] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
      der[i] = df_dp[0][i]*dq2_dq1(patch,dp[0],dir_e,i);
  }
  else if (flg[1] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
      der[i] = df_dp[1][i]*dq2_dq1(patch,dp[1],dir_e,i); 
  }
  else if (flg[2] == YES)
  {
    /* OpenMP_1d_Pragma(omp parallel for) */
    for (i = 0; i < nn; ++i)
      der[i] = df_dp[2][i]*dq2_dq1(patch,dp[2],dir_e,i); 
  }


  for (d = 0; d < 3; ++d)
    if (df_dp[d])
      free(df_dp[d]);
  
  return der;
}

/* taking derivative for f(N) = sum_{0}^{n-1}a_n*T_n(N), in the specified
// patch.
// note: N is in [-1,1].
// ->return value: df(N)/dN.
*/
static double *derivative_ChebyshevExtrema_Tn_1stOrder(Field_T *const f,const Dd_T dir)
{
  assert(dir <= _N2_);
  make_coeffs_1d(f,dir);
  
  Patch_T *const patch = f->patch;
  const unsigned *const n = patch->n;
  const unsigned nn = total_nodes_patch(patch);
  const unsigned B = n[dir]-1;
  double *der = alloc_double(nn);
  double *x = make_1Dcollocation_ChebExtrema(n[dir]);
  const double *const coeffs = f->v2;
  unsigned l;
  
  if (dir == 0)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 1; c < B; ++c)
      {
        unsigned C = L(n,c,j,k);//coeff_ind(i,j,k,c,n,dir);
        der[l] += c*coeffs[C]*Cheb_Un((int)c-1,x[i]);
      }
      der[l] *= 2;
      der[l] += B*coeffs[L(n,B,j,k)]*Cheb_Un((int)B-1,x[i]);
    }
  }
  else if (dir == 1)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 1; c < B; ++c)
      {
        unsigned C = L(n,i,c,k);//coeff_ind(i,j,k,c,n,dir);
        der[l] += c*coeffs[C]*Cheb_Un((int)c-1,x[j]);
      }
      der[l] *= 2;
      der[l] += B*coeffs[L(n,i,B,k)]*Cheb_Un((int)B-1,x[j]);
    }
  }
  else /* (dir == 2) */
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 1; c < B; ++c)
      {
        unsigned C = L(n,i,j,c);//coeff_ind(i,j,k,c,n,dir);
        der[l] += c*coeffs[C]*Cheb_Un((int)c-1,x[k]);
      }
      der[l] *= 2;
      der[l] += B*coeffs[L(n,i,j,B)]*Cheb_Un((int)B-1,x[k]);
    }
  }
  free(x);
  
  return der;
}

/* taking derivative for f(N) = c(0)+2*\sum_{j=1}^{n-1} c(j)Tj(N), in the specified
// patch for Chebyshev basis with Chebyshev nodes collocation.
// note: N is in (-1,1).
// ->return value: df(N)/dN.
*/
static double *derivative_ChebyshevNodes_Tn_1stOrder(Field_T *const f,const Dd_T dir)
{
  assert(dir <= _N2_);
  make_coeffs_1d(f,dir);
  
  Patch_T *const patch = f->patch;
  const unsigned *const n = patch->n;
  const unsigned nn = patch->nn;
  const unsigned B = n[dir];
  double *der = alloc_double(nn);
  double *x = make_1Dcollocation_ChebNodes(n[dir]);
  const double *const coeffs = f->v2;
  unsigned l;
  
  if (dir == 0)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 1; c < B; ++c)
      {
        unsigned C = L(n,c,j,k);//coeff_ind(i,j,k,c,n,dir);
        der[l] += c*coeffs[C]*Cheb_Un((int)c-1,x[i]);
      }
      der[l] *= 2;
    }
  }
  else if (dir == 1)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 1; c < B; ++c)
      {
        unsigned C = L(n,i,c,k);//coeff_ind(i,j,k,c,n,dir);
        der[l] += c*coeffs[C]*Cheb_Un((int)c-1,x[j]);
      }
      der[l] *= 2;
    }
  }
  else /* (dir == 2) */
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 1; c < B; ++c)
      {
        unsigned C = L(n,i,j,c);//coeff_ind(i,j,k,c,n,dir);
        der[l] += c*coeffs[C]*Cheb_Un((int)c-1,x[k]);
      }
      der[l] *= 2;
    }
  }
  free(x);
  
  return der;
}

/* taking second order derivative for f(N) = sum_{0}^{n-1}a_n*T_n(N), in the specified
// patch.
// note: N is in [-1,1].
// ->return value: d^2f(N)/d^2N.
*/
static double *derivative_ChebyshevExtrema_Tn_2ndOrder(Field_T *const f,const Dd_T dir)
{
  assert(dir <= _N2_);
  make_coeffs_1d(f,dir);
  
  Patch_T *const patch = f->patch;
  const unsigned *const n = patch->n;
  const unsigned nn = total_nodes_patch(patch);
  const unsigned B = n[dir]-1;
  double *der = alloc_double(nn);
  double *x = make_1Dcollocation_ChebExtrema(n[dir]);
  const double *const coeffs = f->v2;
  unsigned l;
  
  if (dir == 0)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 2; c < B; ++c)
      {
        unsigned C = L(n,c,j,k);
        der[l] += coeffs[C]*d2T_dx2((int)c,x[i]);
      }
      der[l] *= 2;
      der[l] += coeffs[L(n,B,j,k)]*d2T_dx2((int)B,x[i]);
    }
  }
  else if (dir == 1)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 2; c < B; ++c)
      {
        unsigned C = L(n,i,c,k);
        der[l] += coeffs[C]*d2T_dx2((int)c,x[j]);
      }
      der[l] *= 2;
      der[l] += coeffs[L(n,i,B,k)]*d2T_dx2((int)B,x[j]);
    }
  }
  else /* (dir == 2) */
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 2; c < B; ++c)
      {
        unsigned C = L(n,i,j,c);
        der[l] += coeffs[C]*d2T_dx2((int)c,x[k]);
      }
      der[l] *= 2;
      der[l] += coeffs[L(n,i,j,B)]*d2T_dx2((int)B,x[k]);
    }
  }
  free(x);
  
  return der;
}

/* taking second order derivative for f(N) = sum_{0}^{n-1}a_n*T_n(N), in the specified
// patch for chebyshev basis with chebyshev nodes collocation.
// note: N is in (-1,1).
// ->return value: d^2f(N)/d^2N.
*/
static double *derivative_ChebyshevNodes_Tn_2ndOrder(Field_T *const f,const Dd_T dir)
{
  assert(dir <= _N2_);
  make_coeffs_1d(f,dir);
  
  Patch_T *const patch = f->patch;
  const unsigned *const n = patch->n;
  const unsigned nn = total_nodes_patch(patch);
  const unsigned B = n[dir];
  double *der = alloc_double(nn);
  double *x = make_1Dcollocation_ChebNodes(n[dir]);
  const double *const coeffs = f->v2;
  unsigned l;
  
  if (dir == 0)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 2; c < B; ++c)
      {
        unsigned C = L(n,c,j,k);
        der[l] += coeffs[C]*d2T_dx2((int)c,x[i]);
      }
      der[l] *= 2;
    }
  }
  else if (dir == 1)
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 2; c < B; ++c)
      {
        unsigned C = L(n,i,c,k);
        der[l] += coeffs[C]*d2T_dx2((int)c,x[j]);
      }
      der[l] *= 2;
    }
  }
  else /* (dir == 2) */
  {
    /* OpenMP_2d_Pragma(omp parallel for) */
    for (l = 0; l < nn; ++l)
    {
      unsigned i,j,k;
      unsigned c;
      
      IJK(l,n,&i,&j,&k);
      
      for (c = 2; c < B; ++c)
      {
        unsigned C = L(n,i,j,c);
        der[l] += coeffs[C]*d2T_dx2((int)c,x[k]);
      }
      der[l] *= 2;
    }
  }
  free(x);
  
  return der;
}


/* based on basis and collocation,
// get the pertinent function for spectral derivative. 
*/
static void get_SpecDerivative_func_1stOrder(const Patch_T *const patch,SpecDerivative_Func_T **func)
{
  unsigned i;
  
  for (i = 0; i < 3; ++i)
  {
    if (patch->basis[i] == Chebyshev_Tn_BASIS 
        && patch->collocation[i] == Chebyshev_Extrema)
    {
      func[i] = derivative_ChebyshevExtrema_Tn_1stOrder;
    }
    else if (patch->basis[i] == Chebyshev_Tn_BASIS 
        && patch->collocation[i] == Chebyshev_Nodes)
    {
      func[i] = derivative_ChebyshevNodes_Tn_1stOrder;
    }
    else
      abortEr("There is no such basis or collocation defined for this function.\n");
  }
}

/* based on coordinate system and given direction,
// it find this direction depends on what "a,b,c".
// note: it is ONLY for a,b,c; since it is assumed we want
// to expand or fields on a,b,c.
// for example:
// ============
//
// in spherical coords, if (dir == r) => r = r(r)= r(a)
// in spherical coords, if (dir == x) => x = x(r,theta,phi) = x(a,b,c)
// Notation: dep[?] = 1 means dir depends on ?. 
// and if dep[?] = 0 it means it is not depended.
*/
static void get_dependency(const Patch_T *const patch,const Dd_T dir, unsigned *dep)
{
  unsigned i;
  
  for (i = 0; i < 3; i++)
  {
    dep[i] = 0;
  }
  
  if (patch->coordsys == Cartesian)
  {
    dep[dir%3] = 1;/* means that for example _y_%3 = 1 = _b_ */
  }
  /* x(a,b,c), y(a,b,c), z(a,b,c) */
  else if (patch->coordsys == CubedSpherical) 
  {
    dep[0] = 1;
    dep[1] = 1;
    dep[2] = 1;
  }
  else
     abortEr("There is no coordinate defined for this function.\n");
}

/* based on spectral derivative function and dependencies 
// it finds the Dd_T dp used for direction of partial derivative.
*/
static void get_dp_1stOrder(const Patch_T *const patch,SpecDerivative_Func_T **func,const Dd_T dir,Dd_T *dp)
{
  unsigned depend[3];
  unsigned i;
  
  get_dependency(patch,dir,depend);
  
  for (i = 0; i < 3; ++i)
  {
    dp[i] = UNDEFINED_DIR;
    if (depend[i])/* if this direction depends on _a_,_b_,_c_ */
    {
      if (func[i] == derivative_ChebyshevExtrema_Tn_1stOrder)
        dp[i] = i;/* means _N0_or _N1_or _N2_ */
      else if (func[i] == derivative_ChebyshevNodes_Tn_1stOrder)
        dp[i] = i;/* means _N0_or _N1_or _N2_ */
      else
        abortEr("There is no such derivative function defined for this function.\n");

    }
  }
}

/* making collocation point for Chebyshev Extrema in [-1,1]
// Note: it allocates memory.
//-> return value: Chebyshev Extrema points
*/
static double *make_1Dcollocation_ChebExtrema(const unsigned N)
{
    const double t0 = M_PI/(N-1);
    double *x = alloc_double(N);
    unsigned i;
    
    for (i = 0; i < N; i++)
      x[i] = cos(i*t0);

  return x;
}

/* making collocation point for Chebyshev Nodes in [-1,1]
// Note: it allocates memory.
//-> return value: Chebyshev Nodes points
*/
static double *make_1Dcollocation_ChebNodes(const unsigned N)
{
    const double t0 = M_PI/N;
    double *x = alloc_double(N);
    unsigned i;
    
    for (i = 0; i < N; i++)
      x[i] = cos((i+0.5)*t0);

  return x;
}

